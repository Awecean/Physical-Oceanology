% Physical Oceanology
% Final Project
% B11501037 Po-Tao, Lin
% Problem 1.(3) Problem 1.(4)
clear; close all;

%% Problem 1.(3)
% parameter setting
Lx = 6000e3;    %[m] length of x range
Ly = 3000e3;    %[m] length of y range
tau0 = 0.1;     %[N/m^2] wind stress
rho = 1000;     %[kg/m^3] density of seawater
beta =2e-11;    %[1/m-s] beta-parameter

% Domain setting
Nx = 100;
Ny = 100;
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% V(y) Sverdrup transport
coeff = (tau0 * pi) / (beta * rho * Ly); 
Vy = coeff * sin(pi * Y / Ly);  % V

% psi(x,y) V = \partial \pst/\patial x
% and \psi x = Lx = 0;
psi = (Lx-X) .* Vy;

% U and V
[psi_x, psi_y] = gradient(psi,x(2)-x(1),y(2)-y(1));
U = -psi_y;
V = psi_x;

% Figure of problem 1.(3)
% combine contourf and quiver
figure('Position',[100,100,1200,500]);
set(gcf,'Color','White')
contourf(X/1e3, Y/1e3, psi/1e6,10,'linecolor','none');  % 
colormap('jet')
c = colorbar;
c.Label.String = '10^6 m^3/s';
c.Label.FontSize = 14;
hold on;
Nscale = 10;
quiver(X(1:5:end,1:5:end)/1e3, Y(1:5:end,1:5:end)/1e3, U(1:5:end,1:5:end)*Nscale, V(1:5:end,1:5:end)*Nscale,0, 'w','filled');  
quiver(X(10,90)/1e3+200, Y(10,90)/1e3, 10*Nscale, 0,0, 'r','filled');  
text(X(10,90)/1e3+200, Y(10,90)/1e3+100,'10 m^2/s','Color','r','FontSize',12),
xlabel('x (km)','FontSize',14);
ylabel('y (km)','FontSize',14);
title('Transport Streamfuction','FontSize',14);
axis([0,Lx/1e3,0,Ly/1e3])
exportgraphics(gcf,'Fig1streamfunction.pdf','ContentType', 'vector');
%% Problem 1.(4)
% plot the V section of x = 0.5Lx and x = 0.75Lx
x1 = 0.5*Lx;
x2 = 0.75*Lx;
[~,idx1] = min(abs(x-x1));
[~,idx2] = min(abs(x-x2));

% draw three figure in one picture.
figure('Position',[100,100,1300,500])
set(gcf,'Color','White')
t = tiledlayout(1,5);
% plot V profile at x = 0.75 L_x
ax1 = nexttile(1,[1 1]);
plot(V(:,idx1),Y(:,idx1)/1e3);
xlabel('V (m^2/s)','FontSize',12);
ylabel('y (km)','FontSize',12);
title('x = 0.5L_x','FontSize',14);
% plot V profile at x = 0.75 L_x
ax2 = nexttile(2, [1 1]);
plot(V(:,idx2),Y(:,idx2)/1e3);
xlabel('V (m^2/s)','FontSize',12);
yticklabels('')
title('x = 0.75L_x','FontSize',14);
% plot the streamfunction
ax3 = nexttile(3,[1 3]);
contourf(X/1e3, Y/1e3, psi/1e6,10,'linecolor','none');  % 
colormap('jet')
c = colorbar;
c.Label.String = '10^6 m^3/s';
c.Label.FontSize = 14;
hold on;
Nscale = 10;
quiver(X(1:5:end,1:5:end)/1e3, Y(1:5:end,1:5:end)/1e3, U(1:5:end,1:5:end)*Nscale, V(1:5:end,1:5:end)*Nscale,0, 'w','filled');  
quiver(X(10,90)/1e3-200, Y(10,90)/1e3, 10*Nscale, 0,0, 'r','filled');  
text(X(10,90)/1e3-200, Y(10,90)/1e3+100,'10 m^2/s','Color','r','FontSize',12),
xlabel('x (km)','FontSize',14);
yticklabels('')
title('Transport Streamfuction','FontSize',14);
axis tight manual
axis([0 Lx/1e3 0 Ly/1e3])
exportgraphics(gcf,'Fig2Vprofile.pdf','ContentType', 'vector');
