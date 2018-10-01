clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Testing
% ========================================================================

% a = icoSphereMesh(1);
% length(a.x)
% 
% a = icoSphereMesh(2);
% length(a.x)
% 
% a = icoSphereMesh(3);
% length(a.x)
% 
% 
% a = icoSphereMesh(10);
% length(a.x)
% figure
% plot3(a.x,a.y,a.z,'.')
% axis equal


% % % % % % % % n_tgt = 1297; % ~ 5 deg spacing
% % % % % % % n_tgt = 352;  % ~ 10 deg spacing
% % % % % % % % n_tgt = 145;  % ~ 15 deg spacing - should be used
% % % % % % % % n_tgt = 37;   % ~ 30 deg spacing - decent
% % % % % % % % n_tgt = 17;   % ~ 45 deg spacing
% % % % % % % % n_tgt = 5;    % ~ 90 deg spacing
% % % % % % % 
% % % % % % % [x,y,z] = mySphere(n_tgt*2);
% % % % % % % figure
% % % % % % % plot3(x,y,z,'.','markersize',16)
% % % % % % % axis equal
% % % % % % % PlotBoi3('x','y','z',14)
% % % % % % % view(0,90)
% % % % % % % 
% % % % % % % points = find(y<=1e-15);
% % % % % % % x_hem = x(points);
% % % % % % % y_hem = y(points);
% % % % % % % z_hem = z(points);
% % % % % % % pts = [x_hem',y_hem',z_hem'];
% % % % % % % n_actual = length(points);
% % % % % % % n_tgt
% % % % % % % n_actual
% % % % % % % 
% % % % % % % figure
% % % % % % % plot3(x_hem,y_hem,z_hem,'.','markersize',16)
% % % % % % % axis equal
% % % % % % % PlotBoi3('x','y','z',14)
% % % % % % % view(-90,90)

[vHats] = vHatHemisphere(352,'-z');
figure
plot3(vHats(:,1),vHats(:,2),vHats(:,3),'.','markersize',16)
axis equal
PlotBoi3('x','y','z',14)

% 
% primary = bodies.jupiter;
% secondary = bodies.europa;
% 
% %%% Normalizing constants
% rNorm = secondary.a;         % n <-> km
% tNorm = 1/secondary.meanMot; % n <-> sec
% vNorm = rNorm / tNorm;       % n <-> km/sec
% 
% %%% Acquire Collinear Lagrange points
% L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
% 
% %%% Jacobi constant of Lagrange point
% [JC_L1] = JacobiConstantCalculator(secondary.MR,L123(1,:),[0,0,0]);
% [JC_L2] = JacobiConstantCalculator(secondary.MR,L123(2,:),[0,0,0]);
% [JC_L3] = JacobiConstantCalculator(secondary.MR,L123(3,:),[0,0,0]);
% 
% %%% Creating ZV contours
% % Contour bounds
% xCont_min = L123(1,1)-5*secondary.R_n;
% xCont_max = L123(2,1)+5*secondary.R_n;
% yCont_min = -secondary.R_n*8;
% yCont_max = secondary.R_n*8;
% 
% % Creating x-y grid
% xs = linspace(xCont_min,xCont_max,750);
% ys = linspace(yCont_min,yCont_max,750);
% [X_xy, Y_xy] = meshgrid(xs,ys);
% clear xs ys
% 
% % Calculating JCs across x-y grid
% JCs_xy = zeros(size(X_xy));
% for xk = 1:size(X_xy,1)
%     for yk = 1:size(X_xy,2)
%         %%% Zero-Velocity Curve
%         zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0] ,[0, 0, 0]);
%         JCs_xy(xk,yk) = zv;
%     end
% end
% 
% figure; hold all
% plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
% % plotBodyTexture3(primary.R./rNorm, [-secondary.MR, 0, 0], primary.img)
% [xyContourPoints,href] = contourf(X_xy,Y_xy,JCs_xy,[JC_L2-.00001, JC_L2-.00001],...
%     'color',colors.std.black,'linewidth',1.5);
% colormap(colors.sch.r6(6,:))
% PlotBoi3('X','Y','Z',14)
% view(0,90)
% camva(9)
% axis equal
% set(gcf,'Position',[440 146 839 652]); % fullscreen mac







