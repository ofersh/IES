function [G,x1,y1,CL] = plot_multigauss(mue,Sigma,dx,dy,LB,UB)
% This function plots a 2D multivariate gaussian distribution given the
% mean vector and the covariance matrix, as well as the x-y intervals 
% Usage example: G = plot_multigauss([0.5;-0.7],[2.0 0.5;0.5 0.5], 0.1, 0.1);
CL=20; %contour lines
x1 = LB:dx:UB;
y1 = LB:dy:UB;
[x y] = meshgrid(x1,y1);
XY = [x(:) y(:)]';
mn = repmat(mue,1,size(XY,2));
mulGauss = dx*dy/((2*pi)*sqrt(det(Sigma)))*exp(-0.5.*(XY-mn)'*inv(Sigma)*(XY-mn));
G = reshape(diag(mulGauss),length(x),length(y));
% figure;
% surf(x,y,G);
% xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
% set(gcf,'color','white'); colorbar;
figure;
contour(x1,y1,G,CL,'ShowText','off');
set(gcf,'color','white'); grid on;
end