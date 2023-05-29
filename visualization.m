% ======================================================================
%>@file visualization.m
%>@brief Script to generate graphs and figures by loading data from .mat
%>files
%>
% ======================================================================


%%  RESULTS VISUALIZATION SCRIPT %%
clear all
close all
clc

%% IMPORT DATA %%

vol = load("CM_vol_0_10_variables.mat");
dam_vol = load("CM_dam_vol_0_10_variables.mat");

%% OBJECTIVE-FUNCTION %%

figure(1)
loop_vol=1:size(vol.FO,2);
loop_dam_vol=1:size(dam_vol.FO,2);
plot(loop_vol,vol.FO,loop_dam_vol,dam_vol.FO)
legend('Only volume constraint','Volume and damage constraints','NumColumns', 1,'Location','Southeast')
xlabel('Iterations'); ylabel('FO');
title('Comparison FO');
set(gcf,'color','w');

%% DAMAGE CONSTRAINT %%

figure(2)
loop_vol=1:size(vol.FO,2);
loop_dam_vol=1:size(dam_vol.FO,2);
plot(loop_vol,vol.CONSTRAINTS(:,2),loop_dam_vol,dam_vol.CONSTRAINTS(:,2))
legend('Only volume constraint','Volume and damage constraints','NumColumns', 1,'Location','Northeast')
xlabel('Iterations'); ylabel('Damage constraint');
title('Comparison damage constraint');
set(gcf,'color','w');

%% VON MISES STRESS DISTRIBUTION %%

figure(3);
subplot(1,2,1); volM=surf(flipud(vol.SIGMAVM(:,:,loop_vol(end)))); axis off; axis equal;  colormap jet; colorbar;
caxis([0, max(max(vol.SIGMAVM(:,:,loop_vol(end))))]); view(2); title('Only volume constraint');
subplot(1,2,2); dam_volM=surf(flipud(dam_vol.SIGMAVM(:,:,loop_dam_vol(end)))); axis off; axis equal; colormap jet; colorbar;view(2)
caxis([0, max(max(vol.SIGMAVM(:,:,loop_vol(end))))]); view(2);title('Volume and damage constraints');
set(volM,'edgecolor','none')
set(dam_volM,'edgecolor','none')
sgtitle('von Mises stress distribution')
set(gcf,'color','w');

%% DISPLACEMENTS
U=DISPLACEMENTS;

% name='gif_damage_6000_250_u.gif';

%Plot displacements
fig=figure('Units','normalized','Position',[0.1 0.1 0.7 0.7]);

% for i=4
    surf(flipud(reshape(U1(2:2:end,1),51,101))); axis image; colormap jet; colorbar;     view(2)
    set(gcf,'color','w');grid off; axis([0 nelx 0 nely]); axis equal; axis off;title('Case study 1 - Vertical displacements - 1st it'); caxis([min(U1(2:2:end,1)), max(U1(2:2:end,1))]);drawnow;
    view(2)
    subplot(3,2,2);surf(flipud(reshape(U(2:2:end,1,1),nely+1,nelx+1))); colormap jet; colorbar; set(gcf,'color','w');grid off; axis([0 nelx 0 nely]); axis equal; axis off;title('Gamma-y - 1ST ITERATION'); drawnow;
    view(2)
    subplot(3,2,3);surf(flipud(reshape(U(1:2:end,3,10),nely+1,nelx+1))); colormap jet; colorbar; set(gcf,'color','w');grid off; axis([0 nelx 0 nely]); axis equal; axis off;title('Gamma-x - 10th ITERATION'); drawnow;
    view(2)
    subplot(3,2,4);surf(flipud(reshape(U(2:2:end,3,10),nely+1,nelx+1))); colormap jet; colorbar; set(gcf,'color','w');grid off; axis([0 nelx 0 nely]); axis equal; axis off;title('Gamma-y - 10th ITERATION'); drawnow;
    view(2)
    subplot(3,2,5);surf(flipud(reshape(U(1:2:end,3,20),nely+1,nelx+1))); colormap jet; colorbar; set(gcf,'color','w');grid off; axis([0 nelx 0 nely]); axis equal; axis off;title('Gamma-x - 20th ITERATION'); drawnow;
    view(2)
    subplot(3,2,6);surf(flipud(reshape(U(2:2:end,3,20),nely+1,nelx+1))); colormap jet; colorbar; set(gcf,'color','w');grid off; axis([0 nelx 0 nely]); axis equal; axis off;title('Gamma-y - 20th ITERATION'); drawnow;
    view(2)



    sgtitle('Displacements')
%     frame = getframe(fig); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if i == 1 
%       imwrite(imind,cm,name,'gif', 'Loopcount',inf); 
%     else 
%       imwrite(imind,cm,name,'gif','WriteMode','append'); 
%     end 
% end

%% STRESSES
stress = STRESSES;
name='gif_damage_reference_stresses.gif';

%Plot displacements
fig=figure('Units','normalized','Position',[0.1 0.1 0.7 0.7]);

minx = min(min(STRESSES(:,1,:)));
maxx = max(max(STRESSES(:,1,:)));
miny = min(min(STRESSES(:,2,:)));
maxy = max(max(STRESSES(:,2,:)));
minxy = min(min(STRESSES(:,3,:)));
maxxy = max(max(STRESSES(:,3,:)));

for i=1:1:it

    subplot(1,3,1);imshow(reshape(stress(:,1,i),nely,nelx)); axis off; colorbar; title('Strainxx');colormap jet;set(gcf,'color','w');grid off;caxis([minx maxx]);drawnow;
    subplot(1,3,2);imshow(reshape(stress(:,2,i),nely,nelx)); axis off; colorbar; title('Strainyy');colormap jet;set(gcf,'color','w');grid off;caxis([miny maxy]);drawnow;
    subplot(1,3,3);imshow(reshape(stress(:,3,i),nely,nelx)); axis off; colorbar; title('Strainxy');colormap jet;set(gcf,'color','w');grid off;caxis([minxy maxxy]);drawnow;
    
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,name,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,name,'gif','WriteMode','append'); 
    end 
end

%% VON MISES STRESS

it =92;
stress = STRESSES;
name='gif_ref_6000_stresses_vm.gif';

%Plot displacements
fig=figure('Units','normalized','Position',[0.1 0.1 0.7 0.7]);

for i=1:1:it
    vm_stress(:,i) = sqrt(stress(:,1,i).^2 + stress(:,2,i).^2 - stress(:,1,i).*stress(:,2,i) + 6*stress(:,3,i).^2 );
end

for i=1:92
    
    minv1 = min(min(vm_stress(:,1)));
    maxv1 = max(max(vm_stress(:,1)));

    minv10 = min(min(vm_stress(:,10)));
    maxv10 = max(max(vm_stress(:,10)));

    minv20 = min(min(vm_stress(:,20)));
    maxv20 = max(max(vm_stress(:,20)));

%     subplot(1,3,1);imshow(reshape(vm_stress(:,1),nely,nelx)); axis off; colorbar; title('von Mises stress - 1st iteration');colormap jet;set(gcf,'color','w');grid off;caxis([minv1 maxv1]);drawnow;
%     subplot(1,3,2);imshow(reshape(vm_stress(:,10),nely,nelx)); axis off; colorbar; title('von Mises stress - 10th iteration');colormap jet;set(gcf,'color','w');grid off;caxis([minv10 maxv10]);drawnow;
%     subplot(1,3,3);imshow(reshape(vm_stress(:,20),nely,nelx)); axis off; colorbar; title('von Mises stress - 20th iteration');colormap jet;set(gcf,'color','w');grid off;caxis([minv20 maxv20]);drawnow;
end

%%
name='gif_dam_6000_stresses_vm_dam_300.gif';
stress=STRESSES;nely=50;nelx=50;
%Plot displacements
fig=figure('Units','normalized','Position',[0.1 0.1 0.7 0.7]);
it=227;
for i=1:1:it
    vm_stress(:,i) = sqrt(stress(:,1,i).^2 + stress(:,2,i).^2 - stress(:,1,i).*stress(:,2,i) + 3*stress(:,3,i).^2 );
end


for i=1:it
    minv1 = min(min(vm_stress(:,i)));
    maxv1 = max(max(vm_stress(:,i)));
    subplot(1,2,1);imshow(reshape(vm_stress(:,i),nely,nelx)); axis off; colorbar; title('von Mises stress ');colormap jet;set(gcf,'color','w');grid off;caxis([minv1 maxv1]);drawnow;
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,name,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,name,'gif','WriteMode','append'); 
    end 
end

minv1 = min(min(sigma_vm));
    maxv1 = max(max(sigma_vm));
 subplot(1,3,1);imshow(sigma_vm); axis off; colorbar; title('von Mises stress - 1st iteration');colormap jet;set(gcf,'color','w');grid off;caxis([minv1 maxv1]);drawnow;

%% PRINCIPAL STRAINS

estrain = STRAINS;

strain_1=0.5*(estrain(:,1,end)+estrain(:,2,end))+( (0.5*(estrain(:,1,end)-estrain(:,2,end))).^2 + (1/2*estrain(:,3,end)).^2 ).^0.5;
strain_2=0.5*(estrain(:,1,end)+estrain(:,2,end))-( (0.5*(estrain(:,1,end)-estrain(:,2,end))).^2 + (1/2*estrain(:,3,end)).^2 ).^0.5;

plot(strain_2, strain_1,'.')

%% PRINCIPAL STRESSES

load('ref_6000_250_variables.mat');
estress = STRESSES;

estress_1=0.5*(estress(:,1,end)+estress(:,2,end))+( (0.5*(estress(:,1,end)-estress(:,2,end))).^2 + (1/2*estress(:,3,end)).^2 ).^0.5;
estress_2=0.5*(estress(:,1,end)+estress(:,2,end))-( (0.5*(estress(:,1,end)-estress(:,2,end))).^2 + (1/2*estress(:,3,end)).^2 ).^0.5;

figure;
X1 = linspace(-500,500,1001);
X2 = linspace(-500,500,1001);
[XX1,XX2] = meshgrid(X1,X2);
Z = XX1.^2+XX2.^2-XX1.*XX2-250^2;
contour(XX1,XX2,Z,[250,250])
hold on;
plot(estress_2, estress_1,'.')


load('damage_6000_250_variables.mat')
estress = STRESSES;

estress_1=0.5*(estress(:,1,end)+estress(:,2,end))+( (0.5*(estress(:,1,end)-estress(:,2,end))).^2 + (1/2*estress(:,3,end)).^2 ).^0.5;
estress_2=0.5*(estress(:,1,end)+estress(:,2,end))-( (0.5*(estress(:,1,end)-estress(:,2,end))).^2 + (1/2*estress(:,3,end)).^2 ).^0.5;

figure;
X1 = linspace(-500,500,1001);
X2 = linspace(-500,500,1001);
[XX1,XX2] = meshgrid(X1,X2);
Z = XX1.^2+XX2.^2-XX1.*XX2-300^2;
contour(XX1,XX2,Z,[300,300],'k')
set(gcf,'color','w')
hold on;
xlabel('Minor stress')
ylabel('Major stress')
title('Damage solution')
plot(estress_2, estress_1,'.')


% h=0.9:0.1:1.1;
clear all
delta=50;tau=0.01;phi=0.1;

ratio=0.1:0.01:1.2;
h=ratio-1;
for i=1:size(ratio,2)
    if h(i) <= -phi
        dbeta_dhf(i) = 0;
    elseif h(i) > -phi & h(i) < 0
        dbeta_dhf(i) = exp(((delta*tau^2)/(phi^(2*phi/tau)))*(h(i)+phi)^(2*phi/tau))*((2*phi*delta*tau^2)/(tau*phi^(2*phi/tau)))*(h(i)+phi)^((2*phi/tau)-1);
    elseif h(i) >= 0
        dbeta_dhf(i) = exp(delta*(h(i)+tau)^2)*2*delta*(h(i)+tau);
    end
end
for i=1:size(ratio,2)
    if ratio(i) <= 0.9
        dbeta_dh(i) = 0;
    elseif ratio(i) > 0.9 & ratio(i) <= 1
        dbeta_dh(i) = 20*exp((7.672705*(ratio(i)-0.9))^20)*(7.672705^20)*(ratio(i)-0.9)^19;
    elseif ratio(i) >1
        dbeta_dh(i) = 2*exp(50*(ratio(i)-0.99)^2)*50*(ratio(i)-0.99);
    end
end
plot(ratio,dbeta_dh,ratio, dbeta_dhf,'.')

%% 
clear all
delta=50;tau=0.01;phi=0.1;

ratio=0.1:0.01:1.2;
h=ratio-1;
for i=1:size(ratio,2)
        if ratio(i) <= (1-phi)
            BetaS_f(i) = 1;
        elseif ratio(i) > (1-phi) & ratio(i) <1
            BetaS_f(i) = exp((delta*tau^2/(phi^(2*phi/tau)))*((h(i)+phi)).^(2*phi/tau));
        elseif ratio(i) >=1
            BetaS_f(i) = exp(50*((ratio(i)-0.99)).^2);
        end
end

for i=1:size(ratio,2)
        if ratio(i) <= (1-phi)
            BetaS(i) = 1;
        elseif ratio(i) > (1-phi) && ratio(i) <1
            BetaS(i) = exp((7.672705*(ratio(i)-0.9)).^20);
        elseif ratio(i) >=1
            BetaS(i) = exp(delta*(h(i)+tau).^2);
        end
end

plot(ratio,BetaS,ratio, BetaS_f,'.')
%% FO

plot(FO)


%% SENSITIVITIES

dc = DSENSITIVITIES;

minv1 = min(min(dc(:,:,1)));
maxv1 = max(max(dc(:,:,1)));

minv20 = min(min(dc(:,:,20)));
maxv20 = max(max(dc(:,:,20)));
minv40 = min(min(dc(:,:,40)));
maxv40 = max(max(dc(:,:,40)));
minv60 = min(min(dc(:,:,60)));
maxv60 = max(max(dc(:,:,60)));
minv80 = min(min(dc(:,:,80)));
maxv80 = max(max(dc(:,:,80)));

subplot(1,5,1);imshow(dc(:,:,1)); axis off; colorbar; title('sensitivities - 1st iteration');colormap jet;set(gcf,'color','w');grid off;drawnow;
subplot(1,5,2);imshow(dc(:,:,20)); axis off; colorbar; title('sensitivities- 20th iteration');colormap jet;set(gcf,'color','w');grid off;caxis([minv20 maxv20]);drawnow;
subplot(1,5,3);imshow(dc(:,:,40)); axis off; colorbar; title('sensitivities - 40th iteration');colormap jet;set(gcf,'color','w');grid off;caxis([minv40 maxv40]);drawnow;
subplot(1,5,4);imshow(dc(:,:,60)); axis off; colorbar; title('sensitivities - 60th iteration');colormap jet;set(gcf,'color','w');grid off;caxis([minv60 maxv60]);drawnow;
subplot(1,5,5);imshow(dc(:,:,80)); axis off; colorbar; title('sensitivities - 80th iteration');colormap jet;set(gcf,'color','w');grid off;caxis([minv80 maxv80]);drawnow;





subplot(1,3,1);imshow(mafa); axis off; colorbar; title('sensitivities damage constraint - 1st iteration');colormap jet;set(gcf,'color','w');grid off;caxis([min(min(mafa)) max(max(mafa))]);drawnow;


%% FO and constraints evolution
loop=1:105;
plot(loop,FO,loop,CONSTRAINTS(:,1),loop,CONSTRAINTS(:,2))
legend('FO','Volume constraint','Damage constraint','NumColumns', 1,'Location','Southwest')
xlabel('Iterations');
set(gcf,'color','w');

%% Stress evolution at specific points
loop=1:size(FO,2);
for i=1:size(loop,2)
    str(i,1)=SIGMAVM(40,1,i);
    str(i,2)=SIGMAVM(39,1,i);
    str(i,3)=SIGMAVM(41,2,i);
    str(i,4)=SIGMAVM(41,3,i);
    tt(i,1)=XPHYS(40,1,i);
    tt(i,2)=XPHYS(39,1,i);
    tt(i,3)=XPHYS(41,2,i);
    tt(i,4)=XPHYS(42,3,i);
end



figure;plot(loop,str(:,1),loop,str(:,2),loop,str(:,3),loop,str(:,4))
figure;plot(loop,tt(:,1),loop,tt(:,2),loop,tt(:,3),loop,tt(:,4))
legend('Point 1','Point 2','Point 3','Point 4','NumColumns', 4,'Location','southeast')
xlabel('Iterations'); ylabel('von Mises stress')
set(gcf,'color','w');

%% SENSITIVITIES
dg=DSENSITIVITIES;

for i=1:105
surf(flipud(SIGMAVM(:,:,end))); colorbar;axis equal; axis off; grid off;view(2); drawnow; colormap jet;
end
