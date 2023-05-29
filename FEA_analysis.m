% ======================================================================
%>@file FEA_analysis.m
%>@brief Establishes the Finite Element analysis
%>
%>@param Filters (@b struct) Structure with the data to filter
%>@param Data (@b struct) Initial data that defines the problem
%>(nelx,nely,def,loc,vol,penal,Kin,Kout,sigma_vm_max,t,A0,a,L,Fin,Vh,rnmin,nu,E0,Emin)
%>@param loop (@b integer) Number of iteration
%>@retval OptData (@b struct) Objective function value and FO sensitivities
%>@retval FeaData (@b struct) Displacements, Stresses, Strains, BC's
%>
%>@details
% =====================================================================
function [OptData, FeaData]=FEA_analysis(Filters, Data, FeaData, OptData, loop)

% dout(1) = dout(1)-1;
%Mechanism design - Displacement inverter BC
% F = sparse(2*(Data.nely+1)*(Data.nelx+1),2);
% U = zeros(2*(Data.nely+1)*(Data.nelx+1),2);
% din = 2*(Data.nely+1)-1;
% dout = ((Data.nely+1)*(Data.nelx+1)*2)-1;
% F(din,1) = Fin;
% F(dout,2) = -1;

% fixeddofs = union(1:1:22,2*(Data.nely+1):2*(Data.nely+1):2*(Data.nelx+1)*(Data.nely+1));
% alldofs = 1:2*(Data.nely+1)*(Data.nelx+1);
% freedofs = setdiff(alldofs,fixeddofs);
%     J(1,1) = 0.25; J(1,end) = 0.25; J(end,1) =0.25; J(end,end) =0.25; 
%     J(2:end-1,end)=0.5; J(2:end-1,1)=0.5; J(end,2:end-1)=0.5; J(1,2:end-1)=0.5; 


%% Linear Finite Element Analysis
if Data.analysis == 'l'
    
    % Initial Global Stiffness Matrix, Load and Displacement Vectors 
    [~,D,]=lk(Data);
    [B, ~] = computeB(1, Data, zeros(3,4));

    [KE] = computeKE(B,zeros(4,4),zeros(4,8),D,Data);

    [Kt] = stiffAssemb(Data,OptData,KE);

    J = ones(size(FeaData.din,2),size(FeaData.din,2));
    Kt(FeaData.din,FeaData.din)=Kt(FeaData.din,FeaData.din) + Data.Kin*J; 
    Kt(FeaData.dout,FeaData.dout)=Kt(FeaData.dout,FeaData.dout)+ Data.Kout;
    Kt=(Kt+Kt')./2;
    U(FeaData.freedofs,:)=Kt(FeaData.freedofs,FeaData.freedofs)\FeaData.F(FeaData.freedofs,:);
    U(FeaData.fixeddofs,:)=0;
    U1=U(:,1); 
    U2=U(:,2); 
    U3=U(:,3);
    if Data.ndis >1
        U4=U(:,4); 
        U5=U(:,5); 
        U6=U(:,6);
    end


%     %% Compute von Mises stress
%     stresses = zeros(Data.nelx*Data.nely,3);
%     strains = zeros(Data.nelx*Data.nely,3);
%     Data.sigma_vm = zeros(Data.nely,Data.nelx);
%     hv1=1;
%     for elx=1:Data.nelx
%         for ely=1:Data.nely
%             edof = Data.edofMat(hv1,:);
%             Ue = U(edof,1);
%             [estrains, estresses, d_estress] = computeStrainStress_linear(B,D,Ue, Data.t,OptData.vxPhys(ely,elx), Data.penal);
%             stresses(hv1,:) = estresses';
%             strains(hv1,:)=estrains';
%             d_stresses{hv1} = d_estress;
%             hv1 = hv1 + 1;
%         end
%     end
%     sx = stresses(:,1); sy = stresses(:,2); sxy = stresses(:,3);
%     sigma_vm = reshape((sx.^2 + sy.^2 -sx.*sy + 3*sxy.^2).^0.5,Data.nely,Data.nelx);
%     
%     
    
    %% Objective-function and sensitivity analysis
    tout = U1(dout(1)); 
    te = reshape(sum((U2(Data.edofMat)*KE).*U1(Data.edofMat),2),Data.nely,Data.nelx);
    dtout = Data.penal*((1-OptData.vxPhys)*Data.Emin^(Data.penal-1)+OptData.vxPhys).*te;
    tin  = sum(U1(din))/length(din);
    te = reshape(sum((U3(Data.edofMat)*KE).*U1(Data.edofMat),2),Data.nely,Data.nelx)/length(din);
    dtin = Data.penal*((1-OptData.vxPhys)*Data.Emin^(Data.penal-1)+OptData.vxPhys).*te;
    dc1 = (dtout./tin - tout./tin^2*dtin); 
    obj1=tout./tin; 
    
    OptData.df0dx = dc1(:);
    OptData.f0val = obj1;
    
    if Data.ndis > 1
        tout = U4(dout(2)); 
        te = reshape(sum((U5(Data.edofMat)*KE).*U4(Data.edofMat),2),Data.nely,Data.nelx);
        dtout = Data.penal*((1-OptData.vxPhys)*Data.Emin^(Data.penal-1)+OptData.vxPhys).*te;
        tin  = sum(U4(din))/length(din);
        te = reshape(sum((U6(Data.edofMat)*KE).*U4(Data.edofMat),2),Data.nely,Data.nelx)/length(din);
        dtin = Data.penal*((1-OptData.vxPhys)*Data.Emin^(Data.penal-1)+OptData.vxPhys).*te;
        dc2 = (dtout./tin - tout./tin^2*dtin); 
        obj2=tout./tin; 
        
        OptData.df0dx = OptData.df0dx + dc2(:);
        OptData.f0val = OptData.f0val + obj2;
    end
    


elseif Data.analysis == 'n'
    
    dc = zeros(Data.nely,Data.nelx);
    obj=0;
    for i = 1:Data.ndis
        [KE,D,u,stresses,strains,d_strains, d_stresses, B_matrix, Kt, Plastic, TOTAL,Kst,KL,Klt,B_g,CT]=nonLinearFEM(Data,Filters,FeaData, OptData,i);
        STRESSES(:,:,:,i) = stresses;
        STRAINS(:,:,:,i) = strains;
        FeaData.U(:,3*(i-1)+1) = u;
        FeaData.Stresses(:,i) = stresses(:,end);
        FeaData.Strains(:,i) = strains(:,end);
        B_M(i,:)=B_matrix;
        KT(:,:,i) = Kt;
        FeaData.SigmaVM(:,i) = (STRESSES(:,end,1).^2 + STRESSES(:,end,2).^2-STRESSES(:,end,1).*STRESSES(:,end,2)+3*STRESSES(:,end,3).^2).^(1/2);
    end
    

    %% Compliant mechanisms 
    %Solve the linear system for the adjoint loads (input and output
    %displacements)
    if Data.ndis == 1
        FeaData.U(FeaData.freedofs,2:3)=Kt(FeaData.freedofs,FeaData.freedofs)\FeaData.F(FeaData.freedofs,2:3);
    elseif Data.ndis > 1
        for i = 1:Data.ndis
            Kt1 = KT(:,:,i);
            FeaData.U(FeaData.freedofs,(3*(i-1)+1:3*(i-1)+3))=Kt1(FeaData.freedofs,FeaData.freedofs)\FeaData.F(FeaData.freedofs,(3*(i-1)+1:3*(i-1)+3));
        end
    end
    
    
    %% Sensitivity analysis
    for i = 1:Data.ndis
        %Objective-funtion and sensitivity analysis
        hv2 = 1;
        for elx=1:Data.nelx
            for ely=1:Data.nely 

                %Degrees of freedom for assembly (3-2-1-4)
                edof = Data.edofMat(hv2,:);
                Ud = FeaData.U(:,3*(i-1)+1:3*(i-1)+3);
                %Compute elements internal forces
    %             Ue=Ud(edof,:);
                B=B_M(i,hv2);
                estress = [STRESSES(hv2,end,1,i);STRESSES(hv2,end,2,i);STRESSES(hv2,end,3,i)];
                Fe = computeInternalForce(B,estress',Data.A0,Data.integration,Data.t);

                %Input and output locations displacements
                tout=Ud(FeaData.dout(i),1);
                tin=sum(Ud(FeaData.din,1))/length(FeaData.din); 

                %Sensitivity computation
                dtout(ely,elx)= Data.penal*((1-OptData.vxPhys(ely,elx))*Data.Emin.^(Data.penal-1)+OptData.vxPhys(ely,elx)).*Ud(edof,2)'*Fe;
                dtin(ely,elx)= Data.penal*((1-OptData.vxPhys(ely,elx))*Data.Emin.^(Data.penal-1)+OptData.vxPhys(ely,elx)).*Ud(edof,3)'*Fe;
                dt(ely,elx) = dtout(ely,elx)./tin - tout./tin^2*dtin(ely,elx); 
                hv2=hv2+1;

            end
        end

        dc = dc + dt;         
        obj = obj + tout./tin;
    end
    
    OptData.f0val = obj;
    OptData.df0dx = dc(:);
    
   
end
%======================================================================
%>@file FEA_analysis.m
%>@brief Establishes the Finite Element analysis
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================
