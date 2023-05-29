% ======================================================================
%>@file computeKE_ep.m
%>@brief Compute elemental stiffness matrix in elastoplasticity
%>
%>@param B (@b integer) Number of elements in the horizontal direction
%>@param t (@b integer)  Number of elements in the vertical direction
%>@param Ct (@b integer) Minimum radius for the element filter
%>@param A0 (@b integer) Minimum radius for the node filter
%>
%>@retval KE_ep 
%>
%>@details
% ======================================================================
function [KE_ep] = computeKE_ep(B, t, Ct, A0)


a1=[-1 1 1 -1];
    a2=[-1 -1 1 1];

    csi=[-11/19 11/19];
    wcsi=[1 1];

    eta=[-11/19 11/19];
    weta=[1 1];

    npi = 0; nnode = 4; coords = [0 1 1 0; 0 0 1 1]; 
    kl = zeros(8,8);
    
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

        Bl = B.BL;
        Bnl = B.BN;
       
        kli=(Bl(:,:,npi)+Bnl(:,:,npi))'*Ct*(Bl(:,:,npi)+Bnl(:,:,npi));
        kl=kl + kli.*detj.*wcsi(i).*weta(k).*t;
        
      end
    end
    
    KE_ep = kl.*A0;
%======================================================================
%>@file computeKE_ep.m
%>@brief Compute elemental stiffness matrix in elastoplasticity
%>@details
%>
%>@author Mafalda Gonçalves
%>@date since 2020
%======================================================================