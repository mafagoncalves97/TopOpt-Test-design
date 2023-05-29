% ======================================================================
%>@file stressStiffU_derivative.m
%>@brief Compute element stress stiffness matrix U-derivative
%>
%>@param B (@b integer) Number of elements in the horizontal direction
%>@param Bg (@b integer)  Number of elements in the vertical direction
%>@param D (@b integer) Deformation type - Minimizaton (-1) or maximization
%>(1)
%>@param t (@b integer) Location of the output displacement (1,2,3)
%>
%>@retval dZdu Stress stiffness matrix U-derivative
%>
%>@details
% ======================================================================
function [dZdu] = stressStiffU_derivative(B,Bg,D,t)

    d_estress = zeros(3,1);
    for u=1:8
    
        a1=[-1 1 1 -1];
        a2=[-1 -1 1 1];
    
        csi=[-11/19 11/19];
        wcsi=[1 1];
    
        eta=[-11/19 11/19];
        weta=[1 1];
    
        npi = 0; nnode = 4; coords = [0 1 1 0; 0 0 1 1]; 

        tt= 0; Uvec = zeros (8 ,1); Uvec (u ,1) = 1;
        
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
            detj =jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1);
    
            Bl = Bg(:,:,npi);
            Bi = B(1).BL(:,:,npi);
            d_estress = d_estress + D(:,:,1)*Bi*Uvec;
            
            %Compute tt
            tt = tt +(Bl'* kron ( eye(2) ,[d_estress(1) ,d_estress(3) ;d_estress(3) ,d_estress(2) ])*Bl).*detj.*wcsi(i).*weta(k).*t;
           
          end
        end
    
        
        dZdu(:,:,u) = tt;
    
    end
end
%======================================================================
%>@file stressStiffU_derivative.m
%>@brief Compute element stress stiffness matrix U-derivative
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================