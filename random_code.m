%% Run multi-start topopt problem

clear; close all; clc

nelx=20; nely=20;                             %Design domain
volfrac = 0.35;                               %Volume fraction
rmin = 1.5;                                   %Filter radius
nG = 10;                                      %Number of grid points
constraints = ['V'];                          %Constraints

BCs = {{1, '-y', 2, '-y', 3, '-y'};
       {2, '-y', 3, '-x', 3, '-y'};
       {1, '-y', 2, '-y', 3 '-y'};
       {2, '-y', 3, '-x', 1, 'x'};
       {2, '-y', 3, '-x', 1, '-y'}};

for i = 1:size(BCs)
    tic 
    bc = BCs(i);    
    loc = cell2mat(bc{1}(1:2:end));
    def = bc{1}(2:2:end);
    ndis = size(bc{1},2)/2;
    %Main function
    boosted_TopOpt(nelx, nely, def, loc, ndis, volfrac, rmin, nG, constraints); 
    
    toc
end