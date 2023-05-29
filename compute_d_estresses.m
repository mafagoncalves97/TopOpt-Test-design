% ======================================================================
%>@file compute_d_estresses.m
%>@brief Compute derivative of stresses related to U
%>
%>@param B (@b integer) 
%>@param Ct (@b integer)  
%>@param t (@b integer) 
%>
%>@retval d_estress (@b array) Element internal force components
%>
%>@details
% ======================================================================
function [d_estress] = compute_d_estresses(B, Ct, t)

        S = zeros(2,2);
        d_estress = zeros(3,1);
        d_ESTRAIN = zeros(3,1);
        
       
        
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

                
            d_estress = d_estress + Ct*b.*detj.*wcsi(i).*weta(k).*t;

          end
        end
%======================================================================
%>@file compute_d_estresses.m
%>@brief Compute derivative of stresses related to U
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================