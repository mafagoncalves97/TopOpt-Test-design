% ======================================================================
%>@file Uderivative_dX.m
%>@brief Computes the derivative of U related to X, in a nonlinear FEA
%>
%>@param FeaData (@b integer) Number of elements in the horizontal direction
%>@param Data (@b integer)  Number of elements in the vertical direction
%>@param vxPhys (@b integer) Deformation type - Minimizaton (-1) or maximization
%(1)
%>@param penal (@b integer) Location of the output displacement (1,2,3)
%>@param Emin (@b integer) Location of the output displacement (1,2,3)
%>
%>@retval dU_dX
%>@retval KE_Elastoplastic
%>
%>@details
% ======================================================================
function [dU_dX, KE_Elastoplastic] = Uderivative_dX(FeaData, Data, vxPhys, penal, Emin)

Kt = FeaData.Kglobal;
U = FeaData.U;
CT = FeaData.D;
B_matrix = FeaData.B;
edofMat = Data.edofMat;
nelx = Data.nelx; nely = Data.nely;
t = Data.t; A0=Data.A0;

hv1 = 1;
for elx = 1:nelx
    for ely = 1:nely
        Ct = CT(:,:,hv1);
        Ct = Ct./((1-vxPhys(ely,elx))*Emin^penal + vxPhys(ely,elx));
        CT_solid(:,:,hv1) = Ct;
        hv1 = hv1 + 1;
    end
end

hv1 = 1;
dK_dX = zeros(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
for elx = 1:nelx
    for ely = 1:nely
        edof = edofMat(hv1,:);
        Ct = CT(:,:,hv1);
        B = B_matrix(hv1);
        [KE_ep] = computeKE_ep(B, t, Ct, A0); 
        dKt_dX(edof,edof) = dK_dX(edof,edof) + ((1-vxPhys(ely,elx))*Emin^(penal-1) + vxPhys(ely,elx))*KE_ep;
        hv1 = hv1 + 1;
        KE_Elastoplastic(:,:,hv1) = KE_ep;
    end
end


dU_dX = Kt\(dKt_dX*U);

%======================================================================
%>@file Uderivative_dX.m
%>@brief Computes the derivative of U related to X, in a nonlinear FEA
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================