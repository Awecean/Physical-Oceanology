% Main 
% This is the main program to call other function to solve the Problems
% B11501037 Po-Tao, Lin
clear; close all;
%% Part 0 load data
filename = 'mercatorglorys12v1_gl12_mean_1993_2016_01.nc';
lon = ncread(filename, 'longitude');    % longitude (degree)
lat = ncread(filename, 'latitude');     % latitude (degree)
depth  = ncread(filename, 'depth');     % depth (m)
zos = ncread(filename, 'zos');          % SSH 
uo = ncread(filename, 'uo');            % eastward velocity
vo = ncread(filename, 'vo');            % north velocity

% sort the longitude by 180
lon360 = mod(lon, 360);
[lon360_sorted, idx] = sort(lon360);

zos_sorted = zos(idx, :);
uo_sorted = uo(idx, :,:);
vo_sorted = vo(idx, :,:);

%% Part 1
close all
figure('Position', [100, 100, 700, 400])

set(gcf, 'Color', 'white');
contourf(lon360_sorted, lat, zos_sorted.', 'Levelstep',0.1);
xlabel('lon(degree)','FontSize',14);
ylabel('lat(degree)','FontSize',14);
title('SSH(m) in January','FontSize',14);
axis([0 360 -90 90])
%ylim([-90, 90]);
colorbar;
clim([-1,1]);
exportgraphics(gcf, 'Fig1a1.pdf', 'ContentType', 'image');

%%% make similar picture with a black region
%lon = 120 to 260, lat = 5 to 50
hold on 
rectangle('Position', [120, 5, 140, 45], 'LineWidth',2);
hold off
exportgraphics(gcf, 'Fig1a2.pdf', 'ContentType', 'image');

%% Part 2-(i)
%onlyplot the region
[~, index_lon_w] = min(abs(lon360_sorted-120));
[~, index_lon_e] = min(abs(lon360_sorted-260));
[~, index_lat_s] = min(abs(lat-5));
[~, index_lat_n] = min(abs(lat-50));
figure('Position', [100, 100, 1200, 550])
set(gcf, 'Color', 'White');
lon_region = lon360_sorted(index_lon_w:index_lon_e);
lat_region = lat(index_lat_s:index_lat_n);
eta_region = squeeze(zos_sorted(index_lon_w:index_lon_e,index_lat_s:index_lat_n)).';
contourf(lon_region, lat_region, eta_region, 'Levelstep',0.1);
hold on
dx = 20; % take 1 data every 20 point
[lon_grid, lat_grid] = meshgrid(lon360_sorted(index_lon_w:dx:index_lon_e), lat(index_lat_s:dx:index_lat_n));
u_data = squeeze(uo_sorted(index_lon_w:dx:index_lon_e,index_lat_s:dx:index_lat_n,1)).';
v_data = squeeze(vo_sorted(index_lon_w:dx:index_lon_e,index_lat_s:dx:index_lat_n,1)).';
quiver(lon_grid, lat_grid, u_data, v_data,6,'k-')
colorbar
clim([-1,1]);
axis([120 260 5 50])
title(sprintf('SSH with velocity at depth = %.2f [m] (in January)', depth(1)),'FontSize',14);
xlabel('lon (degree)','FontSize',14);
ylabel('lat (degree)','FontSize',14);
exportgraphics(gcf, 'Fig2a1.pdf', 'ContentType', 'image');
%%% plot the estimated gyre
%%% the answer if problem 2.(i),estimated the range of gyre
rectangle('Position', [120, 15, 120, 25],'LineWidth',2,'LineStyle','--', 'EdgeColor','r');
exportgraphics(gcf, 'Fig2b1.pdf', 'ContentType', 'image');
hold off


%% Part 2-(ii) estimate the velocity magnitude of gyre
V_mag = sqrt(uo_sorted(:,:, 1).^2+vo_sorted(:,:, 1).^2);
[~, index_lon_gyre_w] = min(abs(lon360_sorted-120));
[~, index_lon_gyre_e] = min(abs(lon360_sorted-240));
[~, index_lat_gyre_s] = min(abs(lat-5));
[~, index_lat_gyre_n] = min(abs(lat-40));
figure('Position', [100, 100, 1200, 550])
set(gcf, 'Color', 'White');
contourf(lon360_sorted(index_lon_gyre_w:index_lon_gyre_e), ...
    lat(index_lat_gyre_s:index_lat_gyre_n),...
    squeeze(V_mag(index_lon_gyre_w:index_lon_gyre_e,index_lat_gyre_s:index_lat_gyre_n)).',...
    'Levelstep',0.2);
colorbar
title(sprintf('velocity magnitude [m/s]'),'FontSize',14);
xlabel('lon (degree)','FontSize',14);
ylabel('lat (degree)','FontSize',14);
exportgraphics(gcf, 'Fig2b2.pdf', 'ContentType', 'image');

%% Part 2-(iii) estimate the kuroshio
V_mag = sqrt(uo_sorted(:,:, 1).^2+vo_sorted(:,:, 1).^2);
[~, index_lon_kuroshio_w] = min(abs(lon360_sorted-124));
[~, index_lon_kuroshio_e] = min(abs(lon360_sorted-131));
[~, index_lat_kuroshio_s] = min(abs(lat-25));
[~, index_lat_kuroshio_n] = min(abs(lat-31));
%create the region of Kuroshio
lon_kuroshio = lon360_sorted(index_lon_kuroshio_w:index_lon_kuroshio_e);
lat_kuroshio = lat(index_lat_kuroshio_s:index_lat_kuroshio_n);
V_mag_kuroshio = V_mag(index_lon_kuroshio_w:index_lon_kuroshio_e,...
    index_lat_kuroshio_s:index_lat_kuroshio_n);

figure('Position', [100, 100, 1200, 550])
set(gcf, 'Color', 'White');
contourf(lon_kuroshio, lat_kuroshio, squeeze(V_mag_kuroshio).',...
    'Levelstep',0.2);
colorbar
title(sprintf('velocity magnitude at Kuroshio [m/s]'),'FontSize',14);
xlabel('lon (degree)','FontSize',14);
ylabel('lat (degree)','FontSize',14);
hold on
exportgraphics(gcf, 'Fig2b3.pdf', 'ContentType', 'image');

V_mag_threshold = 0.18;% determine the kuroshio velocity magnitude

% H_lat
target_lon = 127.5;
[~, ilon] = min(abs(lon_kuroshio-target_lon));
vsection = V_mag_kuroshio(ilon,:);
k_great_lat = find(vsection>V_mag_threshold); % velocity greater than 0.2 in lat
k_s = k_great_lat(1); k_n = k_great_lat(end); % latitude range
% W_lon
target_lat = 28;
[~, ilat] = min(abs(lat_kuroshio-target_lat));
usection = V_mag_kuroshio(:,ilat);
k_great_lon = find(usection>V_mag_threshold); % velocity greater than 0.2 in lon
k_w = k_great_lon(1); k_e = k_great_lon(end); % latitude range


line([lon_kuroshio(k_w), lon_kuroshio(k_e)], [lat_kuroshio(ilat), lat_kuroshio(ilat)],'Color','red','Linestyle','-','Linewidth', 2);
text(126.5, 27.6,'W_{lon}','Color','r','FontSize',10,'BackgroundColor','white')
line([lon_kuroshio(ilon), lon_kuroshio(ilon)], [lat_kuroshio(k_s), lat_kuroshio(k_n)],'Color','red','Linestyle','-','Linewidth',2);
text(127.6, 28.9,'H_{lat}','Color','r','FontSize',10,'BackgroundColor','white')
hold off
exportgraphics(gcf, 'Fig2b4.pdf', 'ContentType', 'image');

% estimate the chracteristic width
RE = 6371e3; %earth radius in [m]
W_lon = RE*cos(deg2rad(lat_kuroshio(ilat)))*(deg2rad(lon_kuroshio(k_e)-lon_kuroshio(k_w)));
H_lat = RE*(deg2rad(lat_kuroshio(k_n)-lat_kuroshio(k_s)));
L_kuroshio = W_lon*H_lat/sqrt(W_lon^2+H_lat^2); % width of Kuroshio
fprintf('W_lon = %.2f[km]\n', W_lon/1e3);
fprintf('H_lat = %.2f[km]\n', H_lat/1e3);
fprintf('The width of Kuroshio W_Kuroshio = %.2f [km]\n', L_kuroshio/1e3);

% estimate the characteristic velocity
U_kuroshio = mean(V_mag_kuroshio(ilon,k_great_lat));
fprintf('The characteristic velocity of Kuroshio is %.2f [m/s]\n', U_kuroshio);

% estimate the Rossby number
f_Kuroshio = 2*(2*pi/86400)*sin(deg2rad(lat_kuroshio(ilat))); %Coriolis frequency
RO_kuroshio = U_kuroshio/(f_Kuroshio*L_kuroshio);
fprintf('The Rossby number of Kuroshio is %.3f\n', RO_kuroshio);
%% Part 3 
g = 9.81; % gravitational accleration [m/s^2]
% using same range in Part 2.(i)
lon_sub = deg2rad(lon360_sorted(index_lon_w:index_lon_e));  
lat_sub = deg2rad(lat(index_lat_s:index_lat_n));  

[lon_grid, lat_grid] = meshgrid(lon_sub, lat_sub);  % M x N

eta = zos_sorted(index_lon_w:index_lon_e, index_lat_s:index_lat_n).';

omega = 2*pi/86400;
[dx,dy, deta_inx, deta_iny, f, ug, vg] = deal(zeros(size(eta)));
for i =1:size(eta,1)-1
    for j = 1:size(eta,2)-1
        dx(i,j) = RE*cos(lat_sub(i))*(lon_sub(j+1)-lon_sub(j));
        dy(i,j) = RE*(lat_sub(i+1)-lat_sub(i));
        deta_inx(i,j) = eta(i,j+1)-eta(i,j);
        deta_iny(i,j) = eta(i+1,j)-eta(i,j);
        f(i,j) = 2*omega*sin(lat_sub(i));
    end
end
for i =1:size(eta,1)-1
    for j = 1:size(eta,2)-1
        ug(i,j) = -g/f(i,j)*deta_iny(i,j)/dy(i,j);
        vg(i,j) = g/f(i,j)*deta_inx(i,j)/dx(i,j);
    end
end
ug_data = squeeze(ug(1:20:end,1:20:end));
vg_data = squeeze(vg(1:20:end,1:20:end));

figure('Position', [100, 100, 1200, 550])
set(gcf, 'Color', 'White');
contourf(lon360_sorted(index_lon_w:index_lon_e), lat(index_lat_s:index_lat_n),...
    eta,'LevelStep',0.2);
colorbar
clim([-1,1]);
hold on
quiver(lon360_sorted(index_lon_w:20:index_lon_e), lat(index_lat_s:20:index_lat_n), ug_data, vg_data,6,'k-')
axis([120,240,5,50])
title(sprintf('SSH with geostophic velocity(in January)'),'FontSize',14);
xlabel('lon (degree)','FontSize',14);
ylabel('lat (degree)','FontSize',14);
exportgraphics(gcf, 'Fig3a1.pdf', 'ContentType', 'image');
hold off
%% Part 4


[~, index_lon_w] = min(abs(lon360_sorted-120));
[~, index_lon_e] = min(abs(lon360_sorted-260));
[~, index_lat_s] = min(abs(lat-5));
[~, index_lat_n] = min(abs(lat-50));
figure('Position', [100, 100, 1200, 550])
set(gcf, 'Color', 'White');
lon_region = lon360_sorted(index_lon_w:index_lon_e);
lat_region = lat(index_lat_s:index_lat_n);
eta_region = squeeze(zos_sorted(index_lon_w:index_lon_e,index_lat_s:index_lat_n)).';
contourf(lon_region, lat_region, eta_region, 'Levelstep',0.1);
hold on
dx = 20; % take 1 data every 20 point
[lon_grid, lat_grid] = meshgrid(lon360_sorted(index_lon_w:dx:index_lon_e), lat(index_lat_s:dx:index_lat_n));

u_data_depth = squeeze(uo_sorted(index_lon_w:dx:index_lon_e,index_lat_s:dx:index_lat_n,23)).';
v_data_depth = squeeze(vo_sorted(index_lon_w:dx:index_lon_e,index_lat_s:dx:index_lat_n,23)).';
quiver(lon_grid, lat_grid, u_data_depth, v_data_depth,6,'k-')
colorbar
clim([-1,1]);
axis([120 260 5 50])
title(sprintf('SSH with velocity at depth = %.2f [m] (in January)', depth(23)),'FontSize',14);
xlabel('lon (degree)','FontSize',14);
ylabel('lat (degree)','FontSize',14);
exportgraphics(gcf, 'Fig4a1.pdf', 'ContentType', 'image');

figure
uo_surf_origin = uo_sorted(index_lon_w:index_lon_e,index_lat_s:index_lat_n,1).';
vo_surf_origin = vo_sorted(index_lon_w:index_lon_e,index_lat_s:index_lat_n,1).';
contour(lon360_sorted(index_lon_w:index_lon_e), lat(index_lat_s:index_lat_n), (ug-uo_surf_origin).^2+(vg-vo_surf_origin).^2);
colorbar
surfaceerror = nansum(nansum((ug-uo_surf_origin).^2)+sum((vg-vo_surf_origin).^2));
fprintf('the error of surface layer is %.2e\n', surfaceerror);
figure
uo_depth_origin = uo_sorted(index_lon_w:index_lon_e,index_lat_s:index_lat_n,23).';
vo_depth_origin = vo_sorted(index_lon_w:index_lon_e,index_lat_s:index_lat_n,23).';
contourf(lon360_sorted(index_lon_w:index_lon_e), lat(index_lat_s:index_lat_n),  (ug-uo_depth_origin).^2+(vg-vo_depth_origin).^2);
colorbar
deptherror = nansum(nansum((ug-uo_depth_origin).^2)+sum((vg-vo_depth_origin).^2));
fprintf('the error of depth layer is %.2e\n', deptherror);