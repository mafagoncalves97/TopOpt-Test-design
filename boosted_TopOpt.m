% ======================================================================
%>@file boosted_TopOpt.m
%>@brief Optimizes the topology of a heterogeneous mechanical test
%>
%>@param nelx (@b integer) Number of elements in the horizontal direction
%>@param nely (@b integer)  Number of elements in the vertical direction
%>@param def (@b integer) Deformation type - Minimizaton (-1) or maximization
%>(1)
%>@param loc (@b integer) Location of the output displacement (1,2,3)
%>@param multi (@b integer) Number of output displacements
%>@param vol (@b integer) Volume fraction for the specimen
%>@param rmin (@b integer) Minimum radius for the filter
%>@param nG (@b integer) Number of grid points
%>@param constraints (@b cell) constraints to be used in the optimization
%>problem (volume, damage and buckling)
%>
%>@retval Test topology
%>
%>@details
%>There is the option for saving gif and variables of the process. It is
%>required comment and uncomment section
% ======================================================================
function boosted_TopOpt(Input_data)
fig = figure;
%%Initialization
[Data, FeaData, OptData] = init_data(Input_data);

%%Prepare filter for elements and nodals
[Filters] = filters(Data);

%%Start loop
loop = 0; change = 1; tol = 1;
while (change > Data.tolx || tol>0.005) && loop < Data.maxloop
    loop = loop+1;
    
    %%FEA analysis | Objective-function and sensitivity computation
    [OptData, FeaData]=FEA_analysis(Filters, Data, FeaData, OptData, loop);
    
    %%Constraints
    if any(Data.constraints == 'V')
        [OptData.fval,OptData.dfdx] = VolumeConstraint(Data, Filters,OptData);
    end
    if any(Data.constraints == 'D')
        [dg,~] = DamageConstraint(Data, FeaData, OptData, Filters);
    end
    if any(Data.constraints == 'B')
        [b,db, ~] = BucklingConstraint(OptData,Data,FeaData);
    end
  
    %%Optimizer
    [OptData] = OptimizerRun(Data, OptData, Filters,loop);
    
    %%SEMDOT analysis
    [OptData,Terr,xg,ls] = SEMDOT_analysis(Data, OptData, Filters);
    
    
    %% CHECK CONVERGENCE
    OptData.xPhys2 = max(0.001, OptData.vxnew);
    change = sum(abs(OptData.xPhys2(:)-OptData.vx(:)))/(Data.vol_frac*Data.nely*Data.nelx);
    tol = Terr/Data.nele;
    OptData.vx = OptData.vxnew;
    OptData.beta = OptData.beta+OptData.ER; %Update heaviside smooth function
    
    %% PLOT RESULTS
    %Iterative results
    fprintf('It.:%d Obj.:%9.4f Vol.:%.3f ch.:%2.5f Topo.:%7.5f\n\n',loop,full(OptData.fval),mean(OptData.vxPhys(:)),change,tol);
    

    %Topologies update at each iteration
    subplot(1,2,1);contourf(Data.fnx, flipud(Data.fny), real(xg)-ls, [0 0]); set(gcf,'color','w');axis equal;  axis tight; axis off; 
    OptData.xPhys2=max(0.0010,OptData.vxPhys);
    subplot(1,2,2);    
    colormap(gray); imagesc(1-real(OptData.xPhys2)); caxis([0,1 ] ) ; axis equal; axis off; drawnow;

    %%Save file
    SaveFile(fig,loop,FeaData,OptData,Input_data.filename)
end
%======================================================================
%>@file boosted_TopOpt.m
%>@brief Script to design a heterogeneous mechanical test
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================