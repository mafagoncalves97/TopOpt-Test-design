function [Kt] = stiffAssemb(Data,OptData,KE)

vxPhys = OptData.vxPhys;
Kt=sparse(2*(Data.nelx+1)*(Data.nely+1),2*(Data.nelx+1)*(Data.nely+1));
i=1;
for elx=1:Data.nelx
    for ely=1:Data.nely
        edof = Data.edofMat(i,:);
        Kt(edof,edof) = Kt(edof,edof) + ((1-vxPhys(ely,elx))*Data.Emin^Data.penal+vxPhys(ely,elx)) *KE; 
        i=i+1;
    end
end