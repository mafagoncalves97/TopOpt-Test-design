% ======================================================================
%>@file computeKE.m
%>@brief Compute element stiffness matrix
%>
%>@param B (@b integer)
%>@param t (@b integer)
%>@param M (@b integer) 
%>@param G (@b integer) 
%>@param D (@b integer) 
%>@param integration_type (@b integer) 
%>@param A0 (@b integer) 
%>@param analysis_type (@b integer) 
%>
%>@retval Kt (@b matrix) Global stiffness matrix
%>
%>@details
% ======================================================================
function [Kt] = computeKE(B,M,G,D,Data)
% Element stiffness matrix computation
% Returns the element matrix for two types of integration (reduced and
% full)

t = Data.t;

if Data.integration == 'f'
    
    a1=[-1 1 1 -1];
    a2=[-1 -1 1 1];

    csi=[-11/19 11/19];
    wcsi=[1 1];

    eta=[-11/19 11/19];
    weta=[1 1];

    npi = 0; nnode = 4; coords = [0 1 1 0; 0 0 1 1]; 
    kl = zeros(8,8); kn = zeros(8,8); ks = zeros(8,8); kns = zeros(8,8);
    knI = zeros(8,8); knJ = zeros(8,8); knK = zeros(8,8);
    
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
        
        %Compute KL 
        kli=Bl(:,:,npi)'*D*Bl(:,:,npi);
        kl=kl + kli.*detj.*wcsi(i).*weta(k).*t;
        
        if Data.analysis == 'n'
            %Compute KS
            ksi = G(:,:,npi)'*M*G(:,:,npi);
            ks = ks + ksi.*detj.*wcsi(i).*weta(k).*t;

            %Compute KN
            kni=Bl(:,:,npi)'*D*Bnl(:,:,npi);
            knI=knI + kni.*detj.*wcsi(i).*weta(k).*t;

            knj=Bnl(:,:,npi)'*D*Bl(:,:,npi);
            knJ=knJ + knj.*detj.*wcsi(i).*weta(k).*t;

            knk=Bnl(:,:,npi)'*D*Bnl(:,:,npi);
            knK=knK + knk.*detj.*wcsi(i).*weta(k).*t;

            %Compute KNs
            kns = kns + (0.5 *Bl(:,:,npi)'*D*Bnl(:,:,npi)+Bnl(:,:,npi)'*D*Bl(:,:,npi)+0.5*Bnl(:,:,npi)'*D*Bnl(:,:,npi)).*detj.*wcsi(i).*weta(k).*t;
        end
      end
    end
    
    kn = knI + knJ + knK;
    if Data.analysis == 'n'
        
        Kt(:,:,1) = kl*Data.A0;
        Kt(:,:,2) = ks*Data.A0;
        Kt(:,:,3) = kn*Data.A0;

        Ks(:,:,1)=kl;
        Ks(:,:,2)=kns;
        
    else
        Kt = kl * Data.A0;
    end

elseif Data.integration == 'r'
    
    KL=B.BL'*D*B.BL*Data.A0;
    KS=G'*M*G*Data.A0;
    KN=(B.BL'*D*B.BN+B.BN'*D*B.BL+B.BN'*D*B.BN)*Data.A0;
    
    Kt(:,:,1) = KL;
    Kt(:,:,2) = KS;
    Kt(:,:,3) = KN;
    
end
%======================================================================
%>@file computeKE.m
%>@brief Compute element stiffness matrix
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================