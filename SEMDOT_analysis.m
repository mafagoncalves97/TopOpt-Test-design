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
function [OptData,Terr,xg,ls] = SEMDOT_analysis(Data, OptData, Filters)

    
    %% ASSIGN FILTERED ELEMENTAL VOLUME FRACTIONS TO NODAL DENSITIES
    xn = reshape((Filters.Hn*OptData.vxPhys(:)./Filters.Hns),Data.nely+1,Data.nelx+1);
    
    %% UPDATE POINT DESIGN BY A HEAVISIDE SMOOTH/STEP FUNCTION
    xg = interp2(Data.nodex,Data.nodey,xn,Data.fnx,Data.fny,'linear');
    l1 =0; l2 = 1;
    while (l2-l1) > 1.0e-5
        ls = (l1+l2)/2.0;
        xgnew = max(0.001,(tanh(OptData.beta*ls)+tanh(OptData.beta*(xg-ls)))/(tanh(OptData.beta*ls)+tanh(OptData.beta*(1-ls))));
        if sum(sum(xgnew))/((Data.ngrid*Data.nelx+1)*(Data.ngrid*Data.nely+1)) - sum(OptData.vxPhys(:))/(Data.nelx*Data.nely) > 0
            l1 = ls;
        else
            l2 = ls;
        end
    end
    
    %% ASSEMBLE GRID POINTS TO ELEMENTS
    OptData.vxPhys(:) = 0;
    Terr = 0;
    Tm=[];
    for i = 1:Data.nelx
        for j = 1:Data.nely
            e = (i-1)*Data.nely + j;
            for i1 = Data.ngrid*(i-1)+1:Data.ngrid*i+1
                for j1 = Data.ngrid*(j-1)+1:Data.ngrid*j+1
                    Tm = [Tm;xgnew(j1,i1)];
                    OptData.vxPhys(e) = OptData.vxPhys(e)+xgnew(j1,i1);
                end
            end
            if min(Tm)>0.001 && max(Tm)<1
                Terr = Terr+1;
            end
            Tm = [];
        end
    end
    OptData.vxPhys = OptData.vxPhys/(Data.ngrid+1)^2;
    

%======================================================================
%>@file OptimizerRun.m
%>@brief Runs the volume constraint
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================