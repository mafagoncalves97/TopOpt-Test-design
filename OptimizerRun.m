% ======================================================================
%>@file OptimizerRun.m
%>@brief Runs the volume constraint
%>
%>@param Data (@b integer) Number of elements in the horizontal direction
%>@param FeaData (@b integer)  Number of elements in the vertical direction
%>@param vxPhys (@b integer) Deformation type - Minimizaton (-1) or maximization
%>
%>@retval b, db Damage constraint and its sensitivity
%>
%>@details
% ======================================================================
function [OptData] = OptimizerRun(Data, OptData, Filters,loop)

    OptData.df0dx = Filters.H*(OptData.df0dx(:)./Filters.Hs);
    

    
    %%  UPDATE DESIGN VARIABLES AND FILTERED ELEMENTAL VOLUME FRACTIONS
    OptData.vxval = reshape(OptData.vx,Data.nele,1);
    [vxmma,~,~,~,~,~,~,~,~,Data] = mmasub(Data, OptData, loop);
    OptData.vxnew = reshape(vxmma,Data.nely,Data.nelx);
    OptData.vxPhys(:) = (Filters.H*OptData.vxnew(:))./Filters.Hs;
    OptData.vxold2 = OptData.vxold1(:);
    OptData.vxold1 = OptData.vx(:);
    

%======================================================================
%>@file OptimizerRun.m
%>@brief Runs the volume constraint
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================