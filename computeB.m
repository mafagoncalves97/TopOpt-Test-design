% ======================================================================
%>@file computeB.m
%>@brief Compute shape functions derivative matrices (linear and nonlinear)
%>
%>@param it (@b integer) 
%>@param integration_type (@b integer)  
%>@param A (@b integer) 
%>@param a (@b integer) 
%>@param L (@b integer) Number of output displacements
%>
%>@retval B
%>@retval G
%>
%>@details
% ======================================================================
function [B, G] = computeB(it, Data,A)

    L = Data.L; a = Data.a;

    if Data.integration == 'f'

        a1=[-1 1 1 -1];
        a2=[-1 -1 1 1];

        % Aproximando 1/3**0.5 -> 11/19

        csi=[-11/19 11/19];

        eta=[-11/19 11/19];
        
       
        npi = 0; nnode = 4; coords = [0 1 1 0; 0 0 1 1]; 

        %Scan all integration points
        for k=1:2 %eta
          for i=1:2 %csi

            npi=npi+1;

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
            COFACTOR = zeros(2);
            detj =jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1);
            COFACTOR(1,1) = +jac(2,2);
            COFACTOR(1,2) = -jac(2,1);
            COFACTOR(2,1) = -jac(1,2);
            COFACTOR(2,2) = +jac(1,1);
            jacinv = COFACTOR' ./ detj;

            %Construct B Matrix and Change to Golval Ref  
            Bl=zeros(3,8);
            g = zeros(4,8);
            for j=1:4
                dNdce(1,1)=dNdcsi(j);
                dNdce(2,1)=dNdeta(j);
                dNdxy=jacinv*dNdce;
                Bl(1,j*2-1)=dNdxy(1,1);
                Bl(2,j*2)=dNdxy(2,1);
                Bl(3,j*2-1)=dNdxy(2,1);
                Bl(3,j*2)=dNdxy(1,1);

                g(1,j*2-1)=dNdxy(1,1);
                g(2,j*2-1)=dNdxy(2,1);
                g(3,j*2)=dNdxy(1,1);
                g(4,j*2)=dNdxy(2,1);
            end

            BL(:,:,npi)=Bl/L;
            G(:,:,npi) = g;
            Bnl = A * g;
            BN(:,:,npi)=Bnl;
            
          end
        end
        
            B.BL = BL;         
            B.BN = BN;
     

    elseif Data.integration == 'r'

        if it == 1
            BL=1/(4*a)*[-1 0 1 0 1 0 -1 0; 0 -1 0 -1 0 1 0 1; -1 -1 -1 1 1 1 1 -1];
            B.BL = BL;
            
            G = 1/(4*a)* [-1 0 1 0 1 0 -1 0; -1 0 -1 0 1 0 1 0;
                          0 -1 0 1 0 1 0 -1; 0 -1 0 -1 0 1 0 1];
            BN = zeros(3,8);          
            B.BN = BN;
        else
            BL = 1/(4*a)*[-1 0 1 0 1 0 -1 0; 0 -1 0 -1 0 1 0 1; -1 -1 -1 1 1 1 1 -1];
            B.BL = BL;
            
            G = 1/(4*a)* [-1 0 1 0 1 0 -1 0; -1 0 -1 0 1 0 1 0;
                          0 -1 0 1 0 1 0 -1; 0 -1 0 -1 0 1 0 1];
            BN = A*G;
            B.BN = BN;
        end

    end
    
end
%======================================================================
%>@file computeB.m
%>@brief Compute shape functions derivative matrices (linear and nonlinear)
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================