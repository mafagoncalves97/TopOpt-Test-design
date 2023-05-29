% ======================================================================
%>@file BucklingConstraint.m
%>@brief Runs the buckling constraint
%>
%>@param Data (@b integer) Number of elements in the horizontal direction
%>@param FeaData (@b integer)  Number of elements in the vertical direction
%>@param vxPhys (@b integer) Deformation type - Minimizaton (-1) or maximization
%>
%>@retval b, db Damage constraint and its sensitivity
%>
%>@details
% ======================================================================
function     [b,dbdx, dbdx2] = BucklingConstraint(vxPhys,Data,FeaData)

%Import data
nelx = Data.nelx; nely = Data.nely; A0 = Data.A0; t = Data.t;
freedofs = FeaData.freedofs;edofMat = Data.edofMat;
stresses = FeaData.stresses; U = FeaData.U(:,1);
B = FeaData.B;
KE = FeaData.KE;
D = FeaData.D;
penal = Data.penal;
Bg = FeaData.G;



nEig = 2;
pAgg = 80;
lambda_min = 1;

%% Stress stiffness matrix - G

%Elemental stress stifness matrix - ZE
[ZE] = stressStiffMatrix(stresses,Bg, nelx, nely,t, A0,vxPhys,penal);

%Assembly global stress matrix - Z
Z = zeros(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
G_global = zeros(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
Kt = zeros(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
hv1 = 1;
for elx = 1:nelx
    for ely = 1:nely
        edof = edofMat(hv1,:);
        G_global(edof,edof) = G_global(edof,edof) + ((1-vxPhys(ely,elx))*Data.Emin^penal+vxPhys(ely,elx))*ZE(:,:,hv1);
        Kt(edof,edof) = Kt(edof,edof) + ((1-vxPhys(ely,elx))*Data.Emin^penal+vxPhys(ely,elx)) *KE;
        hv1 = hv1 +1;
    end
end

%Derivative of G in relaton to u - dZdu
[dgdu] = stressStiffU_derivative(B,Bg,D,t);

%% Compute pairs (lambda, phi) - (BLFs, Buckling modes)
dKt = decomposition(Kt(freedofs,freedofs),'chol','lower');
matFun = @(v) dKt\(G_global(freedofs,freedofs)*v);
[eigVec,Diag] = eigs(matFun,length(freedofs),nEig+4,'sa');
[mu,ii] = sort(diag(-Diag),'descend');
eivSort = eigVec(:,ii(1:nEig));
phi = zeros(2*(nelx+1)*(nely+1),nEig);
phi(freedofs,:) = eivSort./sqrt(diag(eivSort'*Kt(freedofs,freedofs)*eivSort)');

%% Sensitivity analysis - dmu

adjL = zeros(2*(nelx+1)*(nely+1),nEig);
phiKphi = zeros(nelx*nely,nEig);
phiGphi = zeros(nelx*nely,nEig);
% dgdu = zeros(8,8,8);
dgdx = zeros(8,8,nelx*nely);
dkdx = zeros(8,8,nelx*nely);

hv1 = 1;
for elx = 1:nelx
    for ely = 1:nely
        dgdx(:,:,hv1) = penal * ((1-vxPhys(ely,elx))*Data.Emin^(penal-1)+vxPhys(ely,elx))* ZE(:,:,hv1);
        dkdx(:,:,hv1) = penal * ((1-vxPhys(ely,elx))*Data.Emin^(penal-1)+vxPhys(ely,elx))* KE;
        hv1 = hv1 +1;
    end
end


for j = 1:nEig
    hv1 = 1;
    for elx = 1:nelx
        for ely = 1:nely
            %First term - phi * dkdx * phi
            t = phi(:,j);
            phiKphi(hv1,j) = (t(edofMat(hv1,:))'*dkdx(:,:,hv1))*t(edofMat(hv1,:));
            %Second term - phi * dgdx * phi
            phiGphi(hv1,j) = (t(edofMat(hv1,:))'*dgdx(:,:,hv1))*t(edofMat(hv1,:));
            %Setup of adjoint load
            tmp = zeros(2*(nelx+1)*(nely+1),1);
            for k = 1:8
                tmp(edofMat(hv1,k)) = tmp(edofMat(hv1,k)) + (t(edofMat(hv1,:))'*((1-vxPhys(ely,elx))*Data.Emin^(penal)+vxPhys(ely,elx))*dgdu(:,:,k))*t(edofMat(hv1,:));
            end
            adjL(edofMat(hv1,k),j) = tmp(edofMat(hv1,k));
            hv1 = hv1 +1 ;
        end
    end
end

%Solve the adjoint problem and compute wi'*dkdx*u
adj = zeros(nelx*nely,nEig);
wi = zeros(2*(nelx+1)*(nely+1),nEig);
wi(freedofs,:) = dKt\adjL(freedofs,:);
for j=1:nEig
    hv1 = 1;
    for elx = 1:nelx
        for ely = 1:nely
            vv = wi(edofMat(hv1,:),j);
            Ue = U(edofMat(hv1,:));
            adj(hv1,j) = Ue' * dkdx(:,:,hv1) * vv;
            hv1 = hv1 + 1;
        end
    end
end

dmu = -(phiGphi + mu(1:nEig)'.*phiKphi -adj);

%% KS aggregation function and sensitivity
fKs = @(p,v) max(v) + (1/p)*log(sum(exp(p*(v-max(v)))));
dfKs = @(p,v,dv) sum(exp(p*(v-max(v)))'.*dv,2)./sum(exp(p*(v-max(v))));
muKs = fKs(pAgg,mu(1:nEig));
dmuKs = dfKs(pAgg,mu(1:nEig),dmu);
b = 1 - lambda_min * (1/muKs);
if b>=0
    disp('Altert buckling');
end
dbdx = dmuKs.*lambda_min;
dbdx2 = zeros(nelx*nely,1);
dbdx(:) = Filters.H*(dbdx(:)./Filters.Hs);
%======================================================================
%>@file BucklingConstraint.m
%>@brief Runs the buckling constraint
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================