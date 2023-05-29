% ======================================================================
%>@file nonLinearFEM.m
%>@brief Nonlinear Finite Element Analysis 
%>
%>@param nelx (@b integer) Number of elements in the horizontal direction
%>@param nely (@b integer)  Number of elements in the vertical direction
%>@param def (@b integer) Deformation type - Minimizaton (-1) or maximization
%>(1)
%>@param loc (@b integer) Location of the output displacement (1,2,3)
%>@param multi (@b integer) Number of output displacements
%>@param vol (@b integer) Volume fraction for the specimen
%>@param rmin (@b integer) Minimum radius for the filter
%>@param nG (@b integer) Number of grid points
%>@param constraints (@b cell) constraints to be used in the optimization
%>problem (volume, damage and buckling)
%>
%>@retval Test topology
%>
%>@details
% ======================================================================
function [KE,D,U,STRESSES,STRAINS, d_STRAINS, d_stresses, B_matrix, Kt, Plastic, TOTAL, Kst,KL,Klt, B_g,CT] = nonLinearFEM(Data, Filters, FeaData, OptData,i)

    F = FeaData.F(:,3*(i-1)+1);
    U=FeaData.U(:,3*(i-1)+1);


    % Initial Global Stiffness Matrix nd Displacement Vectors 
    inc = 0.1;
    A = zeros(3,4);
    M = zeros(4,4);
    R = sparse(2*(Data.nely+1)*(Data.nelx+1),1);
    strains = zeros(Data.nelx*Data.nely,3);
    stresses = zeros(Data.nelx*Data.nely,3);
    STRAINS=zeros(Data.nelx*Data.nely,1/inc,3);
    STRESSES=zeros(Data.nelx*Data.nely,1/inc,3);
    Sstrains = zeros(Data.nelx*Data.nely,3);
    Sstresses = zeros(Data.nelx*Data.nely,3);
    Fs = zeros(2*(Data.nely+1)*(Data.nelx+1),1);
    MAT=ones(Data.nelx*Data.nely,1);
    TOTAL = zeros(Data.nelx*Data.nely,3);
    TOTALt=zeros(Data.nelx*Data.nely,1/inc,3);

    
    %Initial parameters
    it = 1;
    plastic=zeros(Data.nelx*Data.nely,1);
    Plastic = zeros(Data.nelx*Data.nely,1/inc);
    AC=zeros(Data.nelx*Data.nely,3);

    % Elasticity tensor
    [~, D, ~]=lk(Data);

    % Shape functions derivative matrix (BL+BN)
    [B, ~] = computeB(1, Data, zeros(3,4));

    deltaF = inc * F(:,1);
    CT=D.*ones(3,3,Data.nelx*Data.nely);

    for i=1:1/inc

        detaU=1; 
        loopFEM=0;
        Fs = Fs + deltaF;
        Ur = zeros(2*(Data.nely+1)*(Data.nelx+1),1);
        Uk = zeros(2*(Data.nely+1)*(Data.nelx+1),1);

        % Iteratively solve FEM - Newton Raphson
        while detaU>0.05
            loopFEM=loopFEM+1; hv1=1;
            Fint=zeros(2*(Data.nelx+1)*(Data.nely+1),1);
            Kt=zeros(2*(Data.nelx+1)*(Data.nely+1),2*(Data.nelx+1)*(Data.nely+1));
            Kst=zeros(2*(Data.nelx+1)*(Data.nely+1),2*(Data.nelx+1)*(Data.nely+1));
            Klt=zeros(2*(Data.nelx+1)*(Data.nely+1),2*(Data.nelx+1)*(Data.nely+1));

            for elx=1:Data.nelx
                for ely=1:Data.nely
                    
                    %Degrees of freedom for assembly (3-2-1-4)
                    edof = Data.edofMat(hv1,:);

                    %Element displacements (displacements of the increment)
                    Ue=Ur(edof,1)-Uk(edof,1);

                    %Changes in the displacement (iteration)
                    Uek = Uk(edof,1);

                    %Global displacements
                    UeG=U(edof,1)+Ur(edof,1)-Uk(edof,1);

                    if it >1
                        B= B_matrix(hv1);
                    end
                    
                    %Elastoplastic matrix
                    [KE, D, ~]=lk(Data);
                    D = ((1-OptData.vxPhys(ely,elx))*Data.Emin.^(Data.penal)+OptData.vxPhys(ely,elx)) * D;

                    
                    %Stress and strain computation
                    ESTRAIN = Sstrains(hv1,:)';
                    plastic_exist=plastic(hv1);
                    total = TOTAL(hv1,:)';
                    [ESTRAIN,estrain,d_estrain_e,ESTRESS,estress,A, M, Ct,plastic_inc, ac,total] = computeStrainStress(hv1, ESTRAIN, B,D,UeG,Uek,Data.matbehav,Data.t,OptData.vxPhys(ely,elx), Data.penal, plastic_exist,total,Data.Emin);
                    strains(hv1,:) = estrain';
                    d_strains(hv1,:) = d_estrain_e';
                    stresses(hv1,:) = estress';
                    d_stresses = zeros(3,8);
                    

                    Sstrains(hv1,:) = ESTRAIN';
                    Sstresses(hv1,:) =ESTRESS';
                    plastic(hv1)=  plastic(hv1)+ plastic_inc;
                    AC(hv1,:)=ac';
                    CT(:,:,hv1)=Ct;
                    TOTAL(hv1,:) = total';
                    %Updated shape functions derivative matrix (BL+BN)
                    [B, G] = computeB(it, Data, A);

                    %Element stiffness matrix
                    [K] = computeKE(B,M,G,D,Data); 

                    KL = K(:,:,1); 
                    KS = K(:,:,2); 
                    KN = K(:,:,3);
                    KEt = (KL+KS+KN);

                    %Stiffness matrix assembly
                    Kt(edof,edof) = Kt(edof,edof) + KEt; 
                    Kst(edof,edof) = Kst(edof,edof) + KS;
                    Klt(edof,edof) = Klt(edof,edof) + KL;

                    %Element and global internal forces
                    fe = computeInternalForce(B,ESTRESS',Data.A0,Data.integration,Data.t);
                    Fint(edof) = Fint(edof) + fe(:,1);

                    %Save element B matrices 
                    B_matrix(hv1)=B;
                    B_g(:,:,:,hv1)= G;

                    hv1=hv1+1;

                end
            end

            %Residual vector
            R=Fs(FeaData.freedofs)-Fint(FeaData.freedofs,1);

            %Equilibrium equation
            Kt(FeaData.din,FeaData.din)=Kt(FeaData.din,FeaData.din)+ Data.Kin; 
            Kt(FeaData.dout,FeaData.dout)=Kt(FeaData.dout,FeaData.dout)+ Data.Kout;
            Kt=(Kt+Kt')/2;
            Uk(FeaData.freedofs)=Kt(FeaData.freedofs,FeaData.freedofs)\R;
            Uk(FeaData.fixeddofs)=0;
            Ur = Ur + Uk;
            
            Nconv = reshape((Filters.Hn*OptData.vxPhys(:)./Filters.Hns),Data.nely+1,Data.nelx+1);
            
            if sum(sum(Nconv))>1
                Uconv1=reshape(Uk(1:2:end),Data.nely+1,Data.nelx+1);
                Uconv2=reshape(Uk(2:2:end),Data.nely+1,Data.nelx+1);
                Uconv1(Nconv<0.002)=0;
                Uconv2(Nconv<0.002)=0;
                Uconv=zeros((Data.nelx+1)*(Data.nely+1)*2,1);
                Uconv(1:2:end) = Uconv1(:);
                Uconv(2:2:end) = Uconv2(:);
            end
            
            it = it+1;
            detaU = norm(Uconv)/norm(Ur);

        end
        Plastic(:,i) =  plastic;
        
        %Displacements update
        U(:,1) = U(:,1)+ Ur;
        fprintf('Load: %d U: %.3f\n',i,U(end,1));


        
        STRESSES(:,i,1)=stresses(:,1);
        STRESSES(:,i,2)=stresses(:,2);
        STRESSES(:,i,3)=stresses(:,3);
        
        STRAINS(:,i,1)=strains(:,1);
        STRAINS(:,i,2)=strains(:,2);
        STRAINS(:,i,3)=strains(:,3);
        
        d_STRAINS(:,i,1)=d_strains(:,1);
        d_STRAINS(:,i,2)=d_strains(:,2);
        d_STRAINS(:,i,3)=d_strains(:,3);
        
        
        TOTALt(:,i,1)=TOTAL(:,1);
        TOTALt(:,i,2)=TOTAL(:,2);
        TOTALt(:,i,3)=TOTAL(:,3);
        

    end
    fprintf('Acabou o %d deslocamento \n',j);
   
end
%======================================================================
%>@file nonLinearFEM.m
%>@brief Nonlinear Finite Element Analysis 
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================