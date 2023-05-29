% ======================================================================
%>@file VolumeConstraint.m
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
function [v,dv] = VolumeConstraint(Data, Filters,OptData)

dv = ones(Data.nely,Data.nelx);
dv(:) = Filters.H*(dv(:)./Filters.Hs);
v = sum(OptData.vxPhys(:))/(Data.vol_frac*Data.nely*Data.nelx)-1;
dv(:) = dv(:)'/(Data.vol_frac*Data.nely*Data.nelx);

%======================================================================
%>@file VolumeConstraint.m
%>@brief Runs the volume constraint
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================