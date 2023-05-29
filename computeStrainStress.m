% ======================================================================
%>@file computeStrainStress.m
%>@brief Compute stress and strain of an element considering linear
%>
%>@param B (@b integer) Number of elements in the horizontal direction
%>@param D (@b integer)  Number of elements in the vertical direction
%>@param Ue (@b integer) Deformation type - Minimizaton (-1) or maximization
%>(1)
%>@param t (@b integer) Location of the output displacement (1,2,3)
%>@param vxPhys (@b integer) Number of output displacements
%>@param penal (@b integer) Volume fraction for the specimen
%>
%>@retval ESTRAIN (@b array) Element strain components (Green-Lagrange)
%>@retval ESTRESS (@b array) Element stress components (Green-Lagrange)
%>@retval d_estrain (@b array) Element strain increment components
%>@retval estrain (@b array) Element strain increment components
%>@retval estress (@b array) Element stress increment components
%>@retval A (@b matrix) Element stress increment components
%>@retval M (@b matrix) Element stress increment components
%>@retval Ct (@b matrix) Consistent matrix
%>@retval plastic (@b array) Element equivalent plastic strain
%>@retval ac (@b array) Element plastic direction components
%>@retval total (@b array) Element total strain components
%>
%>@details
% ======================================================================
function [ESTRAIN,estrain,d_estrain,ESTRESS,estress,A, M, Ct,plastic, ac,total] = computeStrainStress(hv1, ESTRAIN, B,D,Ue,Uek, material_behavior,t,vxPhys, penal, plastic_exist,total,Emin)

        S = zeros(2,2);
        d_estrain = zeros(3,1);
        d_ESTRAIN = zeros(3,1);
        
        
        if material_behavior=='l'
            ESTRAIN=zeros(3,1);
            ESTRESS = zeros(3,1);
        end
        
        A = zeros(3,4);
        M = zeros(4,4);
        
        a1=[-1 1 1 -1];
        a2=[-1 -1 1 1];

        csi=[-11/19 11/19];
        wcsi=[1 1];

        eta=[-11/19 11/19];
        weta=[1 1];

        npi = 0; nnode = 4; coords = [0 1 1 0; 0 0 1 1]; 
        
        for k=1:2 %eta
          for i=1:2 %csi

            npi=npi+1;
            
            if length(fieldnames(B)) == 1
                b = B.BL(:,:,npi);    
            elseif length(fieldnames(B)) == 2
                b = B.BL(:,:,npi) + B.BN(:,:,npi);  
            end

            %Derivatives of Shape Functions local referencial
            dNdcsi=zeros(nnode,1);
            dNdeta=zeros(nnode,1);
            for j=1:nnode
                dNdcsi(j)=(1/4).*a1(j).*(1+a2(j).*eta(k));
                dNdeta(j)=(1/4).*a2(j).*(1+a1(j).*csi(i));
            end

            %Determine Jacobian Matrix --------------------
            jac=zeros(2);
            for j=1:nnode
                jac(1,1)=jac(1,1)+dNdcsi(j).*coords(1,j);
                jac(1,2)=jac(1,2)+dNdcsi(j).*coords(2,j);
                jac(2,1)=jac(2,1)+dNdeta(j).*coords(1,j);
                jac(2,2)=jac(2,2)+dNdeta(j).*coords(2,j);
            end

            %Determine Inverse Jacobian Matrix -------------------- 
            detj =jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1);

            if material_behavior == 'l'
                ESTRAIN = ESTRAIN + b*(Ue+Uek).*detj.*wcsi(i).*weta(k).*t;
                ESTRESS = ESTRESS + D*b*(Ue+Uek).*detj.*wcsi(i).*weta(k).*t; 
                Ct = D;
                plastic_exist=plastic_exist+ 0;
                ac = zeros(3,1);
                
            elseif material_behavior == 'n'
                    d_ESTRAIN = d_ESTRAIN + b*Uek.*detj.*wcsi(i).*weta(k).*t;
 
            end
           
          end
        end

        if material_behavior == 'n'
            Ex=ESTRAIN(1); Ey=ESTRAIN(2); Exy=ESTRAIN(3); 
            F = [1 + Ex, Exy; Exy, 1+Ey];
            J=det(F);
            
            %Tranform delta_E em delta_e
            d_Estrain_tensor = [d_ESTRAIN(1) d_ESTRAIN(3); d_ESTRAIN(3) d_ESTRAIN(2)];
            d_estrain_tensor = F'\d_Estrain_tensor/F;
            d_estrain = [d_estrain_tensor(1,1); d_estrain_tensor(2,2); d_estrain_tensor(1,2)];

            total = total + d_estrain;
            
            %Transform E em e
            Estrain_tensor = [Ex Exy; Exy Ey];
            estrain_tensor = F'\Estrain_tensor/F;
            estrain = [estrain_tensor(1,1); estrain_tensor(2,2); estrain_tensor(1,2)];
            
            %Compute cauchy stress and delta elastic strain
            [estress, Ct, plastic,ac,d_estrain_e] = materialBehavior(estrain, d_estrain, D, plastic_exist, vxPhys, penal,Emin, hv1);
            
%             plastic = plastic_inc+plastic;
            %Add delta elastic strain to elastic strain
            estrain = estrain + d_estrain_e;
            
            %Tranform normal elastic strain into Green-Lagrange
            estrain_tensor = [estrain(1), estrain(3); estrain(3), estrain(2)];
            ESTRAIN_tensor = F'*estrain_tensor*F;
            ESTRAIN =[ESTRAIN_tensor(1,1); ESTRAIN_tensor(2,2); ESTRAIN_tensor(1,2)];
            Ex=ESTRAIN(1); Ey=ESTRAIN(2); Exy=ESTRAIN(3); 
            
            %Transform cauchy stress into 2nd Piola-Kirchoff
            sx=estress(1); sy=estress(2); sxy=estress(3); 
            sigma_tensor = [sx sxy; sxy, sy];
            S = J * F\sigma_tensor/F';
            Sx=S(1,1); Sy=S(2,2); Sxy=S(1,2);
            ESTRESS=[Sx; Sy; Sxy];
            
            div=1;
            A = [Ex 0 0.5*Exy 0; 0 0.5*Exy 0 Ey; 0.5*Exy Ex Ey 0.5*Exy];
            M = [Sx/div 0 Sxy/div 0; 0 Sx/div 0 Sxy/div; Sxy/div 0 Sy/div 0; 0 Sxy  0 Sy/div]; 

        else
            
            plastic = plastic_exist;
            
            
            Ex=ESTRAIN(1); Ey=ESTRAIN(2); Exy=ESTRAIN(3); 
            F = [1 + Ex, Exy; Exy, 1+Ey];
            J=det(F);
            
            Estrain_tensor = [Ex Exy; Exy Ey];
            estrain_tensor = F'\Estrain_tensor/F;
            estrain = [estrain_tensor(1,1); estrain_tensor(2,2); estrain_tensor(1,2)];
            total = estrain;
            
            Sx=ESTRESS(1); Sy=ESTRESS(2); Sxy=ESTRESS(3); 
            A = [Ex 0 0.5*Exy 0; 0 0.5*Exy 0 Ey; 0.5*Exy Ex Ey 0.5*Exy];
            M = [Sx 0 Sxy 0; 0 Sx 0 Sxy; Sxy 0 Sy 0; 0 Sxy  0 Sy]; 
        
            ESTRESS = ((1-vxPhys)*Emin.^(penal)+vxPhys) * ESTRESS; 
            Sx=ESTRESS(1); Sy=ESTRESS(2); Sxy=ESTRESS(3); 

            S = [Sx Sxy; Sxy, Sy];
            sigma_tensor = 1/J * F*S*F';
            sx = sigma_tensor(1,1); sy = sigma_tensor(2,2); sxy = sigma_tensor(1,2);
            estress=[sx;sy;sxy];
            
        end
   
end
%======================================================================
%>@file computeStrainStress.m
%>@brief Compute stress and strain of an element considering nonlinear
%>behavior
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================