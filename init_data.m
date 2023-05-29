% ======================================================================
%>@file init_data.m
%>@brief Initialize problem parameters for FEA and TopOpt
%>
%>@param nelx (@b integer) Number of elements in the horizontal direction
%>@param nely (@b integer)  Number of elements in the vertical direction
%>@param rmin (@b integer) Minimum radius for the element filter
%>@param rnmin (@b integer) Minimum radius for the node filter
%>
%>@retval Data (@b struct) Struct with all the problem parameters
%>
%>@details
% ======================================================================
function [Data,FeaData,OptData] = init_data(Input_data)

%%Domain parameters
Data.Fin = Input_data.fin;
Data.Kin=Input_data.spring_x; Data.Kout=Input_data.spring_y; Data.Vh=0.238;
Data.t = 1; 
Data.ndis = Input_data.ndis;

%%FEA/SEMDOT parameters
Data.nelx = Input_data.domainx; Data.nely = Input_data.domainy; 
Data.nele = Data.nely*Data.nelx;
Data.L=1; Data.a=(Data.L)/2; Data.A0 = Data.L*Data.L;
Data.nG = 10; Data.ngrid = Data.nG -1;
nodenrs = reshape(1:(1+Data.nelx)*(1+Data.nely),1+Data.nely,1+Data.nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,Data.nelx*Data.nely,1);
Data.edofMat = repmat(edofVec,1,8)+repmat([0 1 2*Data.nely+[2 3 0 1] -2 -1],Data.nelx*Data.nely,1);
[Data.nodex,Data.nodey] = meshgrid(0:Data.nelx,0:Data.nely);
[Data.fnx,Data.fny] = meshgrid(0:1/Data.ngrid:Data.nelx,0:1/Data.ngrid:Data.nely);
Data.analysis = Input_data.analysis;
Data.matbehav = Input_data.matbehav;
Data.integration = Input_data.integration;
Data.constraints = Input_data.constraints;

%%Material parameters
Data.Emin = 0.001; Data.E0 = Input_data.E0; Data.poisson = Input_data.poisson; 

%%Optimization parameters
Data.vol_frac = Input_data.vol_frac;
Data.penal = Input_data.penal; Data.maxloop = 5000; Data.tolx = 0.005;
Data.rnmin = 1; Data.rmin = Input_data.radius;
Data.m = length(Input_data.constraints); 
Data.vxmin = 1e-3*ones(Data.nele,1); Data.vxmax = ones(Data.nele,1);


%%Optimization variables
OptData.vx = repmat(Data.vol_frac,Data.nely,Data.nelx);
OptData.vxPhys = OptData.vx;
OptData.xPhys2 = max(0.001,OptData.vxPhys);
OptData.f0val = 0;
OptData.df0dx = zeros(Data.nele,1);
OptData.df0dx2 = zeros(Data.nele,1); 
OptData.fval = 0;
OptData.dfdx = zeros(Data.nele,1);
OptData.dfdx2 = zeros(Data.m,Data.nele);
OptData.a0 = 1; OptData.ai = 0.01+zeros(Data.m,1);
OptData.ci = 1000+zeros(Data.m,1); OptData.di = zeros(Data.m,1);
OptData.low = 1e-3*ones(Data.nele,1); OptData.upp = ones(Data.nele,1);
OptData.vxval = reshape(OptData.vx,Data.nele,1);

OptData.vxold1 = reshape(OptData.vxPhys,Data.nele,1); OptData.vxold2 = OptData.vxold1;
OptData.beta = 0.5; OptData.ER = 0.5;

%%Optimization constraints
Data.sigma_vm_max_0=230;

%%FEA variables
FeaData.F = sparse(2*(Data.nely+1)*(Data.nelx+1),3*Data.ndis);
FeaData.U = zeros(2*(Data.nely+1)*(Data.nelx+1),3*Data.ndis);
FeaData.SigmaVM = zeros(Data.nely*Data.nelx,Data.ndis);
FeaData.Stresses = zeros(Data.nely*Data.nelx,Data.ndis);
FeaData.Strains = zeros(Data.nely*Data.nelx,Data.ndis);

%% Boundary conditions

if Input_data.dtype == 'sym'
    FeaData.fixeddofs = union(1:2:2*(Data.nely+1),2*(Data.nely+1):2*(Data.nely+1):2*(Data.nelx+1)*(Data.nely+1));
    FeaData.din = 2*(Data.nelx)*(Data.nely+1)+int32((1-Data.Vh)*Data.nely)*2+1:2:2*(Data.nelx+1)*(Data.nely+1);
elseif Input_data.dtype == 'tot'    
    fix1 = ceil((1-Data.Vh)*Data.nely)-1:2*(Data.nely+1)-ceil((1-Data.Vh)*Data.nely)+2;
    fix2 = 2*(Data.nelx)*(Data.nely+1)+ceil((1-Data.Vh)*Data.nely):2:2*(Data.nelx+1)*(Data.nely+1)-ceil((1-Data.Vh)*Data.nely)+2;
    FeaData.fixeddofs = union(fix1,fix2);
    FeaData.din = 2*(Data.nelx)*(Data.nely+1)+ceil((1-Data.Vh)*Data.nely)-1:2:2*(Data.nelx+1)*(Data.nely+1)-ceil((1-Data.Vh)*Data.nely)+2;
end
FeaData.alldofs = 1:2*(Data.nely+1)*(Data.nelx+1);
FeaData.freedofs = setdiff(FeaData.alldofs,FeaData.fixeddofs);

if Input_data.dtype == 'sym'
    Data.ndis=1;
    FeaData.dout = zeros(1,Data.ndis);
    FeaData.dout(1) = 2*(int32((1-Data.Vh)*Data.nely)+2); %loc2
%FeaData.dout(i) = 2*(int32((1-Data.Vh)*Data.nely)+2)+2*(Data.nely+1)*(int32(Data.nelx/2)); %loc1
%FeaData.dout(i) = 2*(int32((1-Data.Vh)*Data.nely)+2)+2*(Data.nely+1)*(int32(Data.nelx/4)); %loc3
    FeaData.F(FeaData.din,1) = Data.Fin;
    FeaData.F(FeaData.dout(1),2) = -1;
    FeaData.F(FeaData.din(1),3) = 1;
else
    Data.ndis=2;
    FeaData.dout = zeros(1,Data.ndis);
    FeaData.dout(1) = 2*(Data.nely+1)*(Data.nelx/2+1)-2*(Data.nely+1)+int32((1-Data.Vh)*Data.nely);
    FeaData.dout(2) = 2*(Data.nely+1)*(Data.nelx/2+1)-int32((1-Data.Vh)*Data.nely)+2;
    FeaData.F(FeaData.din,1) = Data.Fin;
    FeaData.F(FeaData.dout(1),2) = -1;
    FeaData.F(FeaData.din(1),3) = 1;
    FeaData.F(FeaData.din,4) = Data.Fin;
    FeaData.F(FeaData.dout(2),5) = 1;
    FeaData.F(FeaData.din(end),6) = 1;
end





end
%======================================================================
%>@file init_data.m
%>@brief Initialize problem parameters for FEA and TopOpt
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================