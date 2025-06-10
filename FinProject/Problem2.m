% Physical Oceanology
% Final Project
% B11501037 Po-Tao, Lin
% Problem 2
% This Code is base on the "read_nc_data_SNC.m" teacher provided.


% ECCO v4 ocean state estimate, release 3 -- 1992-2015 climatology  [ftp://ecco.jpl.nasa.gov/Version4/Release3/]
%    nx = 720 : per 0.5 deg. lon
%    ny = 360 : per 0.5 deg. lat
%    nz = 50  : depth(m), uneven-grid
%    nt = 12  : monthly
%  SSH       [nx,ny,nt]     Surface Height Anomaly adjusted with global steric height change and sea-ice load (m)
%  EVEL      [nx,ny,nz,nt]  East-West Comp. of Velocity (m/s)
%  NVEL      [nx,ny,nz,nt]  North-South Comp. of Velocity (m/s)
%  WVELMASS  [nx,ny,nz,nt]  Vertical Comp of Velocity (m/s)
%  oceTAUE   [nx,ny,nt]     Eastward Comp. surface wind stress  (N/m^2)
%  oceTAUN   [nx,ny,nt]     Northward Comp. surface wind stress (N/m^2)
clear; close all;
%%
path_temp = fileparts(mfilename('fullpath'));  % read the filepath of this code
f_nc  = [path_temp,'/DATA_V4R3/NVEL.0001.nc'];
%%
time  = ncread(f_nc,'tim');   % [12x1] :15,46,75,...
depth = ncread(f_nc,'dep');   % [50x1] : 5 15,..
lat2d = ncread(f_nc,'lat');   % [720x360]
lon2d = ncread(f_nc,'lon');   % [720x360]
NVEL  = ncread(f_nc,'NVEL');
f_nc = [path_temp,'/DATA_V4R3/EVEL.0001.nc'];        EVEL  = ncread(f_nc,'EVEL');
f_nc = [path_temp,'/DATA_V4R3/WVELMASS.0001.nc'];    WVEL  = ncread(f_nc,'WVELMASS');
f_nc = [path_temp,'/DATA_V4R3/SSH.0001.nc'];          SSH  = ncread(f_nc,'SSH');
f_nc = [path_temp,'/DATA_V4R3/oceTAUN.0001.nc'];  oceTAUN  = ncread(f_nc,'oceTAUN');
f_nc = [path_temp,'/DATA_V4R3/oceTAUE.0001.nc'];  oceTAUE  = ncread(f_nc,'oceTAUE');


% calc time.mean
vo   = nanmean(    NVEL,4 );     % annual mean meridional velocity
uo   = nanmean(    EVEL,4 );     %             zonal velocity
wo   = nanmean(    WVEL,4 );     %             vertical velocity
zos  = nanmean(     SSH,3 );     %             SSH
Estr = nanmean( oceTAUE,3 );     %             zonal wind stress (\tau_x)
Nstr = nanmean( oceTAUN,3 );     %             meridional wind stress (\tau_y)

  clearvars  EVEL NVEL WVEL SSH oceTAUN oceTAUE;

[nx,ny,nz] = size(vo);

% Pacific Ocean, lon : 120E~90W : 120~270
% re-arrange longitudinal grid from -180~180 to 0~360

lon2d_360 = lon2d; lon2d_360( lon2d<0 ) = lon2d( lon2d<0 )+360;
[tmp_val, tmp_idx] = sort( lon2d_360(:,1) );
lon2d_360 = lon2d_360( tmp_idx,:,:);
lon2d = lon2d( tmp_idx,:,:);
lat2d = lat2d( tmp_idx,:,:);
zos   = zos(   tmp_idx,:,:);
Estr  = Estr(  tmp_idx,:,:);
Nstr  = Nstr(  tmp_idx,:,:);
vo    = vo(    tmp_idx,:,:);
uo    = uo(    tmp_idx,:,:);
wo    = wo(    tmp_idx,:,:);

%---------------------------------
% Earth parameters
g      = 9.8;        % gravity             [m/s^2]
omega  = 7.29e-5;    % Earth rotation rate [1/s]
Rearth = 6378137;    % Earth radius        [m]
f      = 2*omega*sind( lat2d );          % Coriolis parameter          [1/s]
beta   = 2*omega*cosd( lat2d )/Rearth;   % beta = df/dy                [1/(ms)]
rho0   = 1025;                           % reference density of ocean  [kg/m^3]

%---------------------------------
% parameters for vertical integration
nthz1 = 12;  % nth depth-grid ~ depth beneath the Ekman layer
nthz2 = 31;  % nth depth-grid ~ depth of min. w (i.e. level of no motion)

%---------------------------------
% region to check Sverdrup Balance
reg_lx = [150 230];  % Pacific Ocean & away from Eastern & Western boundaries

%---------------------------------
% note for calculating dy & dx
% Example:
% for 1-degree N-S displacement -- dy = 1 * pi/180 * R_earth ~ 111 km
% for 1-degree E-W displacement (along 30N) -- dx = 1 * pi/180 * (R_earth * cos(30*pi/180) ) ~ 96 km 


%---------------------------------
%  plotting example


%% Solve the Problem
figure
contourf(lon2d_360,lat2d,Nstr,-1:0.1:1); hold on;
%contourf(lon2d_360,lat2d,zos,-1:0.1:1); hold on;
clim([-1 1]); colorbar;
line([1 1]*reg_lx(1),[0 90],'color',[1 0 1],'linestyle','--','linewidth',2);
line([1 1]*reg_lx(2),[0 90],'color',[1 0 1],'linestyle','--','linewidth',2);
title('annual mean ssh');
axis equal;  xlabel('lat'); ylabel('lon');

%% Problem 2.(1)
dx = 1 * pi/180 * (Rearth * cos(30*pi/180) );  % E-W spacing (m)
dy = 1 * pi/180 * Rearth;                % N-S spacing (m)

[dtauy_dx, ~] = gradient(Nstr, dx, dy); %
[~, dtaux_dy] = gradient(Estr, dx, dy);

curl_tau = dtauy_dx - dtaux_dy;
% Ekman pumping velocity
Wek = curl_tau ./ (rho0 * f);
%Wek(abs(f) < 1e-5) = NaN;
figure('Position',[50,100,1200,650]);
set(gcf,'Color','White')
t1 = tiledlayout(2,1);

ax1 = nexttile(1,[1 1]);

pcolor(lon2d_360, lat2d, wo(:,:,nthz1));  
c = colorbar;
clim([-1e-6,1e-6])
c.Label.String = 'm^2/s';
c.Label.FontSize = 14;
shading interp;
xticklabels('');
ylabel('latitude (^\circ N)','FontSize',14);
title(sprintf('w at depth z = %.2f m',depth(nthz1)),'FontSize',14);
axis([120, 240, 0, 60])

ax2 = nexttile(2,[1 1]);

pcolor(lon2d_360, lat2d, Wek);
c = colorbar;
clim([-1e-6,1e-6])
c.Label.String = 'm^2/s';
c.Label.FontSize = 14;
shading interp;
xlabel('longitude(^\circ)','FontSize',14);
ylabel('latitude (^\circ N)','FontSize',14);
title('W_{EK}', 'FontSize', 14);
axis([120, 240, 0, 60])

exportgraphics(gcf,'Fig3DistributionOfVelocityAndWindCurl.pdf','ContentType', 'vector');
%% Problem 2.(2)
dz = diff(depth);
dz = [dz; dz(end)];  


dz_slice = reshape(dz(nthz1:nthz2), 1, 1, []);
vg_obs = squeeze(nansum(vo(:,:,nthz1:nthz2) .* dz_slice, 3)); 

%
WEK  = squeeze(wo(:,:,nthz1));  
Wmin = squeeze(wo(:,:,nthz2));  
Wdif = WEK - Wmin;
Wtemp = wo(:,:,nthz1)-wo(:,:,nthz2);
Vtheory = (f ./ beta) .*Wtemp;
vg_theory = (f ./ beta) .* Wdif;

% in here, only take the region of our target(150~230, 20~40N) by mask
region_mask = (lon2d_360 >= 150 & lon2d_360 <= 230) & ...
              (lat2d >= 20 & lat2d <= 40);
vg_obs(~region_mask) = NaN;
vg_theory(~region_mask) = NaN;

v1 = vg_theory(region_mask);    % V_g^{theory}
v2 = vg_obs(region_mask);       % V_g^{obs}

R = corr(v1, v2, 'rows','complete');  
RMSE = sqrt(mean((v1 - v2).^2));     
fprintf('Correlation: %.2f, RMSE: %.2e m^2/s\n', R, RMSE);

%% Plot figure
% Figure 4 the distribution of Vg 
figure('Position',[50,20,1200,650])
set(gcf,'Color','White') 
t2 = tiledlayout(2,1);
ax3 = nexttile(1, [1 1]);
pcolor(lon2d_360, lat2d, Vtheory);
xticklabels('');
ylabel('Latitude(^\circ N)','FontSize',14);
shading interp;
c = colorbar;
c.Label.String = 'm^2/s';
c.Label.FontSize = 14;
title('V_g^{theory}')
clim([-1e1,1e1])
axis([150, 230, 20, 40])


ax4 = nexttile(2, [1 1]);
pcolor(lon2d_360, lat2d, vg_obs);
shading interp;
c = colorbar;
c.Label.String = 'm^2/s';
c.Label.FontSize = 14;
clim([-1e1,1e1])
title('V_g^{obs}')
xlabel('Longitude(^\circ)','FontSize',14);
ylabel('Latitude(^\circ N)','FontSize',14);
axis([150, 230, 20, 40])
exportgraphics(gcf,'Fig4Vdistribution.pdf','ContentType', 'vector');





%% Figure 5 Compare result
figure('Position',[100,100,800,500]);
set(gcf,'Color','White')
scatter(v1, v2, 3, 'filled');
xlabel('V_g^{theory} (m^2/s)');
ylabel('V_g^{obs} (m^2/s)');
title('Sverdrup Balance: Theory vs Observed');
text(-10,10,sprintf('R = %.2f',R),'FontSize',12);
lsline; %least-square line
axis equal; grid on;
exportgraphics(gcf,'Fig5TheResultComparison.pdf','ContentType', 'vector');

