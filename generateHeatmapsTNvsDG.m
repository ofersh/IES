%A script to generate heatmaps of 2D populations of integer vectors sampled
%via the TN and DG distributions.
close all; clear all;
Nsample=10000; LB=-10; UB=10;
sigma1=1;%
sigma2=3; %
c12 = 1.2; %-1.2 %-0.8; %0
disp([sigma1^2 c12;c12 sigma2^2]);
baseC = [sigma1^2 c12;c12 sigma2^2];
baseC = enforceSymmetryAndPSD(baseC)
if baseC(1,1)~=baseC(2,2)
    alpha12=0.5*(atan(2*baseC(1,2)/(baseC(1,1)-baseC(2,2))));
else
    alpha12=0.0;
end

C = anglesToCovariance([sqrt(baseC(1,1)) alpha12; alpha12 sqrt(baseC(2,2))])
ZR=zeros(Nsample,2);
for i=1:Nsample,
   ZR(i,:)=round(corr_sigma([sqrt(baseC(1,1)) alpha12; alpha12 sqrt(baseC(2,2))])); 
end

% DOUBLE-GEOMETRIC
ZZ=zeros(Nsample,2);
for i=1:Nsample
   gc1=geom_correlated([sqrt(baseC(1,1)) alpha12; alpha12 sqrt(baseC(2,2))]); 
   gc2=geom_correlated([sqrt(baseC(1,1)) alpha12; alpha12 sqrt(baseC(2,2))]); 
   ZZ(i,:) = round(gc1) - round(gc2);
end

[G4,x1,y1,CL] = plot_multigauss([0;0],C, 0.2, 0.2,LB,UB);
Mlb=-6; Mub=6;
[binCounts, xEdges, yEdges] = histcounts2(ZZ(:,1), ZZ(:,2));
figure; imagesc(xEdges, yEdges, binCounts');
colormap(flipud(hot)); %colormap('coolwarm');
colorbar;
axis ([Mlb Mub Mlb Mub]); % Ensure proper Cartesian layout
title('Heatmap of DG Sample Histogram');
hold on; 
[cg,hg] = contour(x1,y1,G4,CL,'ShowText','off'); clim([0 800]);
set(gcf,'color','white');
set(gca,'xtick',[Mlb:1:Mub]); 
set(gca,'ytick',[Mlb:1:Mub]); %set(gca,'ytick',[Mub:-1:Mlb]);
set(gca, 'YDir', 'normal');  % Fix y-axis direction
grid on;
hg.LineWidth = 1; hg.LineColor = 'c'; %hg.FaceColor = 'r';
xlabel('$x_1$','interpreter','latex','FontSize',16); ylabel('$x_2$','interpreter','latex','FontSize',16);
saveas(gcf, 'DG2.svg');
%%%%%%%%%%%%%%%%
%figure(20); %TN
[binCounts, xEdges, yEdges] = histcounts2(ZR(:,1), ZR(:,2));
figure; imagesc(xEdges, yEdges, binCounts');
colormap(flipud(hot)); %colormap(flipud(gray)); 
colorbar;
axis ([Mlb Mub Mlb Mub]); % Ensure proper Cartesian layout
% Labels and title
title('Heatmap of TN Sample Histogram');
hold on; 
[cg,hg] = contour(x1,y1,G4,CL,'ShowText','off'); clim([0 800]);
set(gcf,'color','white');
set(gca,'xtick',[Mlb:1:Mub]); 
set(gca,'ytick',[Mlb:1:Mub]);
set(gca, 'YDir', 'normal');  % Fix y-axis direction
grid on;
hg.LineWidth = 1; hg.LineColor = 'c'; %hg.LineColor = 'y'; hg.FaceColor = 'r';
xlabel('$x_1$','interpreter','latex','FontSize',16); ylabel('$x_2$','interpreter','latex','FontSize',16);
saveas(gcf, 'TN2.svg');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%