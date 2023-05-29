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
function SaveFile(fig,loop,FeaData,OptData,filename)

    FO(loop)=OptData.fval;
    CONSTRAINTS(loop,1)=OptData.fval(1,:)';
    SENSITIVITIES(:,:,loop)=OptData.df0dx;
    DISPLACEMENTS(:,:,loop)=FeaData.U;
    STRESSES(:,:,:,loop) = FeaData.Stresses;
    SIGMAVM(:,:,:,loop)=FeaData.SigmaVM;
    STRAINS(:,:,:,loop) = FeaData.Strains;
    XPHYS(:,:,loop)=OptData.vxPhys;

    save(strcat(filename,'.mat'),'FO','CONSTRAINTS','SENSITIVITIES','DISPLACEMENTS','STRESSES','STRAINS','SIGMAVM','XPHYS');

    
        %% SAVE DATA
    %Save gif of optimization process
    name=strcat(filename,'.gif');
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    %Write to the GIF File 
    if loop == 1 
      imwrite(imind,cm,name,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,name,'gif','WriteMode','append'); 
    end 


%======================================================================
%>@file VolumeConstraint.m
%>@brief Runs the volume constraint
%>@details
%>
%>@author Mafalda Gonçalves
%>@date 20-07-2022
%======================================================================