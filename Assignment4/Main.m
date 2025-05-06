% Physical Oceanography HW4
% B11501037 CE3
% Po-Tao, Lin
clear; close all;
%% Part 1 parameter set of problem 1.(3),(4)
rho1 = 1000; % seawater density (kg/m^3)
f1 = 5e-5;   % Coriolis parameter (1/s)               
Az1 = 1e-2;  % vertical dymaic viscosity (m^2/s)
tau_wind_x_1 = 0.1; % wind stress (Pa)
             
Az2 = 1e-3;  % vertical dymaic viscosity (m^2/s)
Az3 = 1e-1;

zn = 300;
z = linspace(0,-300,zn); % z-set (m)
%% Part 2 function set 
function [u,v, delta, flux_x, flux_y, flux_x_theory, flux_y_theory] = Ekman(rho, f, Az, tau, z)
    % u: the x-direction velocity
    % v: the y-direction velocity
    % delta: the Ekman layer thickness
    % flux_x: Ekman transport in x-direction by program
    % flux_y: Ekman transport in y-direction by program
    % flux_x_theory: Ekmantransport in x-direction in theory
    % flux_y_theory: Ekmantransport in y-direction in theory


    delta = sqrt(2*Az/f); % depth of Ekman layer (m)
    C = exp(z/delta);
    factor = tau / (rho*sqrt(Az*f));
    % Ekman sprial
    u = C.*factor.*sin(z/delta + pi/4);
    v = -C.*factor.*cos(z/delta + pi/4);
    % Ekman transport
    temp = -cumtrapz(z,u); flux_x = temp(end); % consider the integral order
    temp = -cumtrapz(z,v); flux_y = temp(end);
    flux_x_theory = 0;
    flux_y_theory = -tau/rho/f;
end
%% Problem 1.(3), 1.(4)
figure('Position',[100,100,500,500]);
set(gcf, 'Color','White');
[u1, v1, delta1, flux_x1, flux_y1, flux_x_theory1, flux_y_theory1] = Ekman(rho1, f1, Az1, tau_wind_x_1, z);
hold on
for i = 1:length(z)
    plot3([0 u1(i)], [0 v1(i)], [z(i) z(i)], 'b');
end
grid on;
xlabel('u (m/s)');ylabel('v (m/s)');zlabel('z (m)');
zlim([-100,0])
title('Ekman Spiral of Problem 1.(3)');
view(3);
exportgraphics(gcf, 'Fig1.pdf', 'ContentType', 'vector');
disp('Problem 1.(3), 1.(4)')
fprintf('Ekman transport in program: flux in x: %.2f, in y: %.2f\n', flux_x1, flux_y1)
fprintf('Ekman transport in theory: flux in x: %.2f, in y: %.2f\n', flux_x_theory1, flux_y_theory1);
%% Problem 1.(5)
%%% Az = 1e-3
figure('Position',[100,100,500,500]);
set(gcf, 'Color','White');
[u2, v2, delta2, flux_x2, flux_y2, flux_x_theory2, flux_y_theory2] = Ekman(rho1, f1, Az2, tau_wind_x_1, z);
hold on
for i = 1:length(z)
    plot3([0 u2(i)], [0 v2(i)], [z(i) z(i)], 'b');
end
grid on;
xlabel('u (m/s)');ylabel('v (m/s)');zlabel('z (m)');
zlim([-35,0])
title(sprintf('Ekman Spiral with A_z = %.1e', Az2));
view(3);
exportgraphics(gcf, 'Fig2a.pdf', 'ContentType', 'vector');
fprintf('Problem 1.(5) with with A_z = %.1e\n', Az2)
fprintf('Ekman layer thickness %.2f\n', delta2);
fprintf('Ekman transport in program: flux in x: %.2f, in y: %.2f\n', flux_x2, flux_y2);
fprintf('Ekman transport in theory: flux in x: %.2f, in y: %.2f\n', flux_x_theory2, flux_y_theory2);
%%% Az = 1e-2
figure('Position',[100,100,500,500]);
set(gcf, 'Color','White');
[u2, v2, delta2, flux_x2, flux_y2, flux_x_theory2, flux_y_theory2] = Ekman(rho1, f1, Az1, tau_wind_x_1, z);
hold on
for i = 1:length(z)
    plot3([0 u2(i)], [0 v2(i)], [z(i) z(i)], 'b');
end
grid on;
xlabel('u (m/s)');ylabel('v (m/s)');zlabel('z (m)');
zlim([-100,0])
title(sprintf('Ekman Spiral with A_z = %.1e', Az1));
view(3);
exportgraphics(gcf, 'Fig2b.pdf', 'ContentType', 'vector');
fprintf('Problem 1.(5) with with A_z = %.1e\n', Az1)
fprintf('Ekman layer thickness %.2f\n', delta2);
fprintf('Ekman transport in program: flux in x: %.2f, in y: %.2f\n', flux_x2, flux_y2);
fprintf('Ekman transport in theory: flux in x: %.2f, in y: %.2f\n', flux_x_theory2, flux_y_theory2);
%%% Az = 1e-1
figure('Position',[100,100,500,500]);
set(gcf, 'Color','White');
[u2, v2, delta2, flux_x2, flux_y2, flux_x_theory2, flux_y_theory2] = Ekman(rho1, f1, Az3, tau_wind_x_1, z);
hold on
for i = 1:length(z)
    plot3([0 u2(i)], [0 v2(i)], [z(i) z(i)], 'b');
end
grid on;
xlabel('u (m/s)');ylabel('v (m/s)');zlabel('z (m)');
zlim([-300,0])
fprintf('Problem 1.(5) with with A_z = %.1e\n', Az3)
view(3);
title(sprintf('Ekman Spiral with A_z = %.1e', Az3));
exportgraphics(gcf, 'Fig2c.pdf', 'ContentType', 'vector');
fprintf('Problem 1.(5) with with A_z = %.1e\n', Az3)
fprintf('Ekman layer thickness %.2f\n', delta2);
fprintf('Ekman transport in program: flux in x: %.2f, in y: %.2f\n', flux_x2, flux_y2);
fprintf('Ekman transport in theory: flux in x: %.2f, in y: %.2f\n', flux_x_theory2, flux_y_theory2);
