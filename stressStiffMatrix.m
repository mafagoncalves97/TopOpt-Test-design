% ======================================================================
%>@file stressStiffMatrix.m
%>@brief Compute element stress stiffness matrix
%>
%>@param stresses (@b integer)
%>@param Bg (@b integer)  
%>@param nelx (@b integer) 
%>@param nely (@b integer) 
%>@param t (@b integer) 
%>@param A0 (@b integer) 
%>
%>@retval GE (@b matrix) Element stress stiffness matrix
%>
%>@details
% ======================================================================
function [GE] = stressStiffMatrix(stresses,Bg,nelx,nely,t,A0,vxPhys,penal)

    X = vxPhys(:);
    for nEle=1:nelx*nely
        xx = X(nEle);
        a1=[-1 1 1 -1];
        a2=[-1 -1 1 1];
    
        csi=[-11/19 11/19];
        wcsi=[1 1];
    
        eta=[-11/19 11/19];
        weta=[1 1];
    
        npi = 0; nnode = 4; coords = [0 1 1 0; 0 0 1 1]; 
        kg = zeros(8,8); 

        sx = stresses(nEle,end,1); sy = stresses(nEle,end,2); sxy = stresses(nEle,end,3); 
        S = [sx, sxy, 0, 0; sxy, sy, 0,0; 0,0,sx,sxy; 0,0,sxy,sy];
        S = S./((1-xx).*0.001^penal+xx);
        
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
            
            %Compute KL 
            kgi=Bl'*S*Bl;
            kg=kg + kgi.*detj.*wcsi(i).*weta(k).*t;
    
           
          end
        end
    
        
        GE(:,:,nEle) = kg*A0;
    
    end
end
%======================================================================
%>@file stressStiffMatrix.m
%>@brief Compute element stress stiffness matrix
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================