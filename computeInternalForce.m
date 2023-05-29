% ======================================================================
%>@file computeInternalForce.m
%>@brief Compute element internal force
%>
%>@param B (@b integer) 
%>@param estress (@b integer)  
%>@param A0 (@b integer) 
%>(1)
%>@param integration_type (@b integer) 
%>@param t (@b integer) 
%>
%>@retval fe (@b array) Element internal force components
%>
%>@details
% ======================================================================
function fe = computeInternalForce(B,estress,A0,integration_type,t)
% Returns the  element internal load for reduced and full integration

    fe=zeros(8,1);

    if integration_type == 'f'
        
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
            else
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
            
            fe = fe + b'*estress'.*detj.*wcsi(i).*weta(k).*t*A0;         
            
            
          end
        end
        

    elseif integration_type == 'r'
        
        B = B.BL + B.BN;
        fe=B'*D*B*Ue*A0;

    end

end
%======================================================================
%>@file computeInternalForce.m
%>@brief Compute element internal force
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================