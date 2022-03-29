% clear
% clc
% close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))

% ========================================================================
%%% Run switches
% ========================================================================
plot_secondary         = 0;
plot_secondaryGradient = 0;
plot_system            = 0;
plot_zvSurface         = 1;
plot_L1                = 0;
plot_L2                = 0;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Setting up system
% ========================================================================

primary = bodies.jupiter;
secondary = bodies.europa;

prms.u = secondary.MR;
prms.n = 1;
prms.R2_n = secondary.R_n;

L123 = EquilibriumPoints(prms.u, prms.n,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

% ========================================================================
%%% Plots
% ========================================================================
% -------------------------------------------------
% Plot Near Secondary
% -------------------------------------------------
if plot_secondary == 1
dvLp_mps = 300; % 450?, 525, 650

% [JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(2,:),[0,0,0]);
JC_Lp = getJacobiConstant_ZH([L123(2,:),0,0,0], prms);
dJC_vel_kps = dvLp_mps/1000;
dJC_Lp = (dJC_vel_kps/vNorm)^2;
JC_scInitial = JC_Lp-dJC_Lp;

%%% Near Secondary
figure; hold all
if isfield(secondary,'img') == 1
    plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
else
    plotBody3(secondary.R_n,[1-secondary.MR,0,0],secondary.color,1)
end
% PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
PlotBoi3_CR3Bn(26)
axis equal
plot3(L123(2,1),0,0,'^','markersize',5,'markeredgecolor',colors.black,'markerfacecolor',colors.grn);
plot3(L123(1,1),0,0,'^','markersize',5,'markeredgecolor',colors.black,'markerfacecolor',colors.grn);
prms.R2_n = secondary.R_n;
if plot_L1 == 1
    plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 1, 0, prms, colors.black, 1.5)
end
if plot_L2 == 1
    plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 2, 0, prms, colors.black, 1.5)
end
% plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 1, 0, prms, colors.black, 1.5)
plotCR3BP_Neck(prms,L123,JC_scInitial,600,200,colors.black,1.5)
view(-38,30)
end % plot_secondary

% -------------------------------------------------
% Plot gradient of energies near secondary
% -------------------------------------------------
if plot_secondaryGradient == 1
dvLp_mps_vec = [50, 150, 350];

gradColors = colorScale([colors.ltred; colors.ltblue], 3);
    
%%% Near Secondary
figure; hold all

for kk = 1:length(dvLp_mps_vec)
    dvLp_mps = dvLp_mps_vec(kk);
    [JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(2,:),[0,0,0]);
    dJC_vel_kps = dvLp_mps/1000;
    dJC_Lp = (dJC_vel_kps/vNorm)^2;
    JC_scInitial = JC_Lp-dJC_Lp;
    
    prms.R2_n = secondary.R_n;
    plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 2, 0, prms, gradColors(kk,:), 1.5)
    plotCR3BP_Neck(secondary,L123,JC_scInitial,600,200,gradColors(kk,:),1.5)
    
end

PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
axis equal
plot3(L123(2,1),0,0,'^','markersize',5,'markeredgecolor',colors.black,'markerfacecolor',colors.grn);
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
view(-38,30)

end % plot_secondary
% -------------------------------------------------
% Plot Bigger System
% -------------------------------------------------
if plot_system == 1
dvLp_mps = 350; % 450?, 525, 650

[JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(2,:),[0,0,0]);
dJC_vel_kps = dvLp_mps/1000;
dJC_Lp = (dJC_vel_kps/vNorm)^2;
JC_scInitial = JC_Lp-dJC_Lp;

%%% System
figure('position',[67 323 560 420]); hold all
%%% Creating ZV contours
% Contour bounds
xCont_min = -1.1;
xCont_max = L123(2,1)+5*secondary.R_n;
yCont_min = -1.1;
yCont_max = 1.1;

% Creating x-y grid
xs = linspace(xCont_min,xCont_max,1500); % 600
ys = linspace(yCont_min,yCont_max,1500); % 200
[X_xy, Y_xy] = meshgrid(xs,ys);
clear xs ys zs

% Calculating JCs across x-y grid
JCs_xy = zeros(size(X_xy));
for xk = 1:size(X_xy,1)
    for yk = 1:size(X_xy,2)
        %%% Zero-Velocity Curve
        zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0] ,[0, 0, 0]);
        JCs_xy(xk,yk) = zv;
    end
end

[xyContourPoints,href] = contour(X_xy,Y_xy,JCs_xy,[JC_scInitial, JC_scInitial],...
    'color',colors.black,'linewidth',1.5);

plot3(L123(2,1),0,0,'^','markersize',5,'markeredgecolor',colors.black,'markerfacecolor',colors.grn);
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
plotBodyTexture3(primary.R/rNorm, [-secondary.MR, 0, 0], primary.img)
plot([L123(2,1) L123(2,1)],[-0.15 0.15],'r--','linewidth',1)
% PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
PlotBoi3_CR3Bn(26)
axis equal
view(0,90)


end % plot_system

% -------------------------------------------------
% Plot full surface of neck around secondary
% -------------------------------------------------
if plot_zvSurface == 1
mu = secondary.MR;

xCont_min = L123(1,1)-2*secondary.R_n;
xCont_max = L123(2,1)+2*secondary.R_n;
yCont_min = -9*secondary.R_n;
yCont_max = 9*secondary.R_n;
zCont_min = -9*secondary.R_n;
zCont_max = 9*secondary.R_n;

% Calculate Jacobi constant
dvLp_mps = 150; % 450?, 525, 650

JC_Lp = getJacobiConstant_ZH([secondary.MR,L123(2,:),0,0,0],prms);
dJC_vel_kps = dvLp_mps/1000;
dJC_Lp = (dJC_vel_kps/vNorm)^2;
JC_scInitial = JC_Lp-dJC_Lp;

% Meshgrid across contour space
n = 1e2;
xbnd = linspace(xCont_min,xCont_max,n);
ybnd = linspace(yCont_min,yCont_max,n);
zbnd = linspace(zCont_min,zCont_max,n);
[x,y,z] = meshgrid(xbnd,ybnd,zbnd);

% Grid of potential 
r1 = sqrt((x+mu).^2 + y.^2 + z.^2);
r2 = sqrt((x-1+mu).^2 + y.^2 + z.^2);
U = 1/2*(x.^2 + y.^2) + (1-mu)./r1 + mu./r2;

% Contour where Cj = 2*U + 0
v = JC_scInitial / 2;

figure
hold on

% plotMoon(1 / c.length, [1-mu, 0 0])


% colormap(jet)
p = patch(isosurface(x, y, z, U, v));
% p.EdgeColor = 0.8 * ones(3,1);
% p.EdgeColor = 'None';
p.EdgeColor = colors.black;
% p.FaceLighting = 'gouraud';
 
ver = p.Vertices;
nVer = length(p.Vertices);
c = zeros(nVer, 1);
for ii = 1:nVer
    c(ii) =  -ver(ii, 3);  
end
p.FaceVertexCData = c;
% p.FaceColor = 'interp';
p.FaceColor = colors.blue;

% plot3([-2 2],[0 0],[0 0], 'k-')
% plot3([0 0],[-2 1],[0 0], 'k-')
% plot3([0 0],[0 0],[-1.5 2], 'k-')

% PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
% axis('equal')
% % camva(7)
% view([30,30])
% % light('Position',[1 -1 0])
% light('Position',[-1 -1 0.5])

end % plot_neckSurface




