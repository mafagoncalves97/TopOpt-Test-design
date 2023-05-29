% ======================================================================
%>@file DamageConstraint.m
%>@brief Runs the damage constraint
%>
%>@param Data (@b integer) Number of elements in the horizontal direction
%>@param FeaData (@b integer)  Number of elements in the vertical direction
%>@param vxPhys (@b integer) Deformation type - Minimizaton (-1) or maximization
%>
%>@retval g, dg Damage constraint and its sensitivity
%>
%>@details
% ======================================================================
function [dg,g_constraint] = DamageConstraint(Data, FeaData, OptData,Filters)
vxPhys = OptData.vxPhys;

%Import data
nelx = Data.nelx; nely = Data.nely; 
freedofs = FeaData.freedofs; fixeddofs = FeaData.fixeddofs;
edofMat = Data.edofMat;
U = FeaData.U(:,1); Emin = Data.Emin; 
KE = FeaData.KE; t = Data.t;
penal = Data.penal; Kt = FeaData.Kglobal;
sigma_vm_max_0 = Data.sigma_vm_max_0;
sigma_vm = reshape(FeaData.sigmaVM,nely,nelx);
d_stresses = FeaData.d_stresses;
material_behavior = FeaData.material_behavior;
stresses = FeaData.stresses;


%% Define damaged model

BetaT = zeros(nely,nelx);
dg = zeros(nely,nelx);
h = zeros(nely,nelx);
ratio = zeros(nely,nelx);
dXdam_dX = zeros(nely,nelx);
dBeta_dX = zeros(nely,nelx);
dfi_dX = zeros(nely,nelx);
dh_dX = zeros(nely, nelx);
dVM_dX = zeros(nely,nelx);
% FF=zeros(242,1);
% 
phi = 0.1; delta=50; tau = 0.01; ep=0.01; q=0.5;
X_damaged = vxPhys;

%% Prepare analysis
if material_behavior == 'n'
    [dU_dx, KE_Elastoplastic] = Uderivative_dX(FeaData,Data, vxPhys, penal, Emin);
    hv1 = 1;
    for elx = 1:nelx
        for ely = 1:nely
            d_stresses(:,:,hv1) = compute_d_estresses(FeaData.B(hv1), FeaData.D(:,:,hv1), t);
            hv1 = hv1 +1;
        end
    end
end
%% Check stresses

%Update damaged model
for elx=1:nelx
    for ely=1:nely
        fi = 1-ep + ep/vxPhys(ely,elx);
        sigma_vm_max = sigma_vm_max_0;
        ratio(ely,elx) = sigma_vm(ely,elx)/(sigma_vm_max*fi);
        h(ely,elx) = sigma_vm(ely,elx)/(sigma_vm_max*fi);
        if ratio(ely,elx) <= (1-phi)
            BetaT(ely,elx) = 1;
        elseif ratio(ely,elx) > (1-phi) && ratio(ely,elx) <1
            BetaT(ely,elx) = exp((6353/828)*(vxPhys(ely,elx)^(q))*(h(ely,elx)-0.9)^20);
        elseif ratio(ely,elx) >=1
            BetaT(ely,elx) = exp(50*((vxPhys(ely,elx)^(q))*(h(ely,elx)-0.99))^2);
        end
  
        X_damaged(ely,elx) = Emin + BetaT(ely,elx) * (vxPhys(ely,elx) - Emin);
    end
end

%% Calculate damage constraint

W_original = sum(sum(vxPhys))/(nelx*nely);
W_damaged = sum(sum(X_damaged))/(nelx*nely);
g_constraint = (W_damaged/W_original)-1;
if g_constraint>0
    fprintf('Damage constraint active - %.3f\n', g_constraint)
end


FF = zeros(2*(nely+1)*(nelx+1),1);




%Volume constraint sensitivity
dv = ones(nely,nelx); 
sx = stresses(:,1); sy = stresses(:,2); sxy = stresses(:,3);
hv1 = 1;
%Damage constraint sensitivity
for elx = 1:nelx
    for ely = 1:nely
        
        edof = edofMat(hv1,:);
        Ue = U(edof,1);
        
        fi = 1-ep + ep/vxPhys(ely,elx);
        sigma_vm_max = sigma_vm_max_0;

        ratio(ely,elx) = sigma_vm(ely,elx)/(sigma_vm_max*fi);
        h(ely,elx) = sigma_vm(ely,elx)/(sigma_vm_max*fi);
        
        if ratio(ely,elx) > (1-phi)
            dVM_dsigma(hv1,:) = [(2*sx(hv1)-sy(hv1))/(2*sigma_vm(ely,elx));(2*sy(hv1)-sx(hv1))/(2*sigma_vm(ely,elx));(3*sxy(hv1))/(sigma_vm(ely,elx))]';
        else
            dVM_dsigma(hv1,:)=[0,0,0];
        end
        
        if material_behavior == 'l'
            dsigma_du = d_stresses{hv1};
            dKE_dX(:,:,hv1) = penal*((1-vxPhys(ely,elx))*Emin^(penal-1)+vxPhys(ely,elx))* KE;
            du_dX = ((vxPhys(ely,elx)^penal)*KE)\(-dKE_dX(:,:,hv1)*Ue);
            
        elseif material_behavior == 'n'
            dsigma_du = d_stresses(:,:,hv1); %D_ep * Bj
            du_dX = dU_dx(edof);
            dKE_dX(:,:,hv1) = penal * ((1-vxPhys(ely,elx))*Emin^(penal-1)+vxPhys(ely,elx))*KE_Elastoplastic(:,:,hv1);
        end
        

        %dB_dX
        if ratio(ely,elx) <= (1-phi)
            dBeta_dX(ely,elx) = 0;
        elseif ratio(ely,elx) > (1-phi) && ratio(ely,elx) <1
            dBeta_dX(ely,elx) = exp(((6353/828)*vxPhys(ely,elx)^(q)*(h(ely,elx)-0.9)).^20) * ((6353/828)*vxPhys(ely,elx)^(q))^20 * 20*(h(ely,elx)-0.9)^19;
        elseif ratio(ely,elx) >=1
            dBeta_dX(ely,elx) = exp(50*(vxPhys(ely,elx)^(q)*(h(ely,elx)-0.99)).^2)*50*vxPhys(ely,elx)^(2*q)*2*(h(ely,elx)-0.99);
        end
        
        %dB_Drho
        if ratio(ely,elx) <= (1-phi)
            dBeta_drho(ely,elx) = 0;
        elseif ratio(ely,elx) > (1-phi) && ratio(ely,elx) <1
            dBeta_drho(ely,elx) = exp(((6353/828)*vxPhys(ely,elx)^(q)*(h(ely,elx)-0.9)).^20) * ((6353/828)*(h(ely,elx)-0.9))^20 * 20 * q*vxPhys(ely,elx)^(20*q-1);
        elseif ratio(ely,elx) >=1
            dBeta_drho(ely,elx) = exp(50*(vxPhys(ely,elx)^(q)*(h(ely,elx)-0.99)).^2)*50*(h(ely,elx)-0.99)^2*2*q*vxPhys(ely,elx)^(2*q-1);
        end

        dfi_drho(ely,elx) =  -ep/(vxPhys(ely,elx)^2);
        
        dsigmavm_drho(ely,elx) = (-sigma_vm(ely,elx)/(sigma_vm_max*fi^2))*dfi_drho(ely,elx);


        dBeta(ely,elx) = dBeta_dX(ely,elx)*dsigmavm_drho(ely,elx) + dBeta_drho(ely,elx);

        dXdam_dX(ely,elx) = BetaT(ely,elx) + (vxPhys(ely,elx)-Emin)*(dBeta(ely,elx));

        dgT(ely,elx) =  (1/vxPhys(ely,elx))*dXdam_dX(ely,elx) - X_damaged(ely,elx)/(vxPhys(ely,elx)^2);
        
        FF(edof,1) = FF(edof,1)+((1/vxPhys(ely,elx)) *(vxPhys(ely,elx)-Emin)*dBeta_dX(ely,elx)*(1/(sigma_vm_max*fi))*dVM_dsigma(hv1,:)*dsigma_du)';
        
        hv1 = hv1 + 1;
    end
end 

lambda(freedofs,:)=Kt(freedofs,freedofs)\FF(freedofs,:);
lambda(fixeddofs,:)=0;

hv1=1;
for elx=1:nelx
    for ely=1:nely
        edof = edofMat(hv1,:);
%         xx(ely,elx) = lambda(edof,1)'*(-penal*vxPhys(ely,elx)^(penal-1)*KE*U(edof,1));
        xx(ely,elx) = lambda(edof,1)'*(-dKE_dX(:,:,hv1)*U(edof,1));
        dg(ely,elx) = xx(ely,elx)+dgT(ely,elx);
        hv1 = hv1+1;

    end
end

%%Output
dg(:) = Filters.H*(dg(:)./Filters.Hs);
OptData.vxPhys = vxPhys;
%======================================================================
%>@file DamageConstraint.m
%>@brief Runs the damage constraint
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================