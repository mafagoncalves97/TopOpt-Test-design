% ======================================================================
%>@file computeStrainStress_linear.m
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
%>@retval ESTRAIN (@b array) Element strain components
%>@retval ESTRESS (@b array) Element stress components
%>@retval d_estress (@b array) Element stress derivative components
%>(without penalization)
%>
%>@details
% ======================================================================
function [ESTRAIN, ESTRESS, d_estress] = computeStrainStress_linear(B,D,Ue, t,vxPhys, penal)

        d_estress = zeros(3,8);
        ESTRAIN=zeros(3,1);
        ESTRESS = zeros(3,1);
        
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
            
            b = B.BL(:,:,npi);    


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

            ESTRAIN = ESTRAIN + b*Ue.*detj.*wcsi(i).*weta(k).*t;
            ESTRESS = ESTRESS + D*b*Ue.*detj.*wcsi(i).*weta(k).*t; 
            
            d_estress = d_estress + D*b.*detj.*wcsi(i).*weta(k).*t;

          end
        end

        d_estress = vxPhys^penal * d_estress;
        ESTRESS = vxPhys^penal * ESTRESS; 
            
end
%======================================================================
%>@file computeStrainStress_linear.m
%>@brief Compute stress and strain of an element considering linear
%>behavior
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================