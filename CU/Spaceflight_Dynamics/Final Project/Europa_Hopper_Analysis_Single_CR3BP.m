clear
clc
close all
addpath('ProjectBin')

% ------------------------------------------------------------------------
%%% Plot Switches
% ------------------------------------------------------------------------
% Plot Europa ECEF?
ECEFplot = 1; % no = 0, yes = 1
scale1 = 8; % Plot Scale (limits = scale1 x E_radius)

% Plot Europa ECI?
ECIplot = 1; % no = 0, yes = 1

% Plot Jupiter System Inertial?
JCIplot = 0; % no = 0, yes = 1

% Plot ECI distance vs time?
rECIplot = 0; % no = 0, yes = 1

% Run movie?
scale2 = 1.4;
framespeed = 100; % Higher is faster
runECEFMovie = 0; % no = 0, yes = 1
runECIMovie = 0; % no = 0, yes = 1
% ------------------------------------------------------------------------
%%% Setting Initial Values
% ------------------------------------------------------------------------
%%% Time Constraints
ti = 0;
tf = 520000; % sec
time = ti:tf;

%%% Jupiter Parameters
uJ = 126672520; % km^3 / s^2
J_pos = [0, 0, 0]; % km
J_radius = 69911; % km

%%% Europa Parameters
E_radius = 1560.8; % km
E_a = 671100; % km
uE = 3203.413216; % km^3 / s^2
nE = sqrt(uJ / (E_a^3)); % rad/s
% Circular Velocity
% rbc = States(1,1:3)*uE/(uE + uJ);
vE = sqrt(uJ/E_a); % km/s
wE = [0; 0; nE]; % Rot. velocity (tidally locked)

%%% Intial Europa State (Jupiter-Centric)
rE0 = [E_a, 0, 0]; % km
vE0 = [0, vE, 0]; % km

%%% Initial Hopper State (Jupiter-Centric) 
% Surface Position (latitude / longitude)
lat1 = 0; % deg
lon1 = -45; % deg
[rH01] = latlon2surfECEF(lat1, lon1, E_radius); % km

%%% Radial Velocity (Europa relative)
% v_mag = 1.9; % km/s (-45, 0)
% v_mag = 1.95; % km/s (-45,0)
v_mag = 1.9; % km/s
vH01 = (rH01/norm(rH01))*v_mag;
% .013, -85 lon 0 lat
% 
% %%% Vector Velocity (Europa relative)
% vH01 = [.2, -vE/7.9, .8]; % km/s (-45,0)
% % vH01 = [-1.7, -.49, .6]; % km/s (-180,0)
% % vH01 = [-1.7, -.2, .6]; % km/s (-180,0)
% % vH01 = [.2, -1.6, 1]; % km/s (-180,0)

%%% Correcting Hopper Velocity 
vH02 = vH01 + cross(wE,rH01); % Adding rotational velocity component
vH03 = vH02 + vE0; % Making velocity relative to Jupiter
%%% Correcting Hopper Position
rH02 = rH01 + rE0; % Making position relative to Jupiter

%%% Setting Initial State Vector
X0 = [rH02 vH03]; % km km/s km km/s

% ------------------------------------------------------------------------
%%% Propagating the State with Numerically Integration
% ------------------------------------------------------------------------
%%% Setting integrator accuracy
tol = 1E-13;
options = odeset('Events',@impactEvent_CR3BP,'RelTol',tol,'AbsTol',tol);

%%% Propagating the State
[Times,States] = ode45(@EJ_EOMIntegrator_CR3BP,time,X0,options,E_radius,uE,uJ,nE,E_a);

% ------------------------------------------------------------------------
%%% Calculating Europa states
% ------------------------------------------------------------------------
r_Europa = zeros(size(Times,1),3);
for k = 1:size(Times,1)
    ta = nE * Times(k); % rads
    r_Europa(k,:) = R3(rE0,ta)'; % km, JCI Europa position vectors
end

%%% Assigning final states of Europa (for simplicity)
rEf = r_Europa(end,:); % km
taf = nE * Times(end); % rad
vEf = R3(vE0,taf)'; % km/s

% ------------------------------------------------------------------------
%%% Calculating rECI positions
% ------------------------------------------------------------------------
rECI_Hopper = States(:,1:3) - r_Europa;

% ------------------------------------------------------------------------
%%% Calculating rECEF positions
% ------------------------------------------------------------------------
rECEF_Hopper = zeros(size(States,1),3);
for k = 1:size(Times,1)
    ta = nE * Times(k); % rad
    rECEF_Hopper(k,:) = R3(rECI_Hopper(k,:),-ta);
end

% ------------------------------------------------------------------------
%%% Calculating Hopper Initial and Final Energies in ECI
% ------------------------------------------------------------------------
v0_m = norm(States(1,4:6) - vE0)*1000; % m/s, ECI
vf_m = norm(States(end,4:6) - vEf)*1000; % m/s, ECI
dVelocity = vf_m - v0_m; % m/s
KE0 = .5*v0_m^2;
KEf = .5*vf_m^2;
dKE = KEf - KE0; % J/kg .... m^2/s^2

% ------------------------------------------------------------------------
%%% Calculating Change in Hopper Surface Angle & Distance Traveled
% ------------------------------------------------------------------------
%%% Defining beginning and ending states
Pos0 = States(1,1:3) - rE0; % ECI / ECEF
Posf = States(end,1:3) - rEf; % ECI
% Rotating Final Position into ECEF
theta = -nE * Times(end);
Posf = R3(Posf,theta)'; % ECEF

%%% Angle change and surface distance traveled
AngChange = atan2(norm(cross(Pos0,Posf)),dot(Pos0,Posf))*180/pi; % deg
Traveled = (AngChange*pi/180)*E_radius*1000; % m

%%% Finding direction traveled
[lat2, lon2] = ECEF2latlon(Posf); % rad
lat2 = lat2 * 180/pi; % deg
lon2 = lon2 * 180/pi; % deg
az = azimuth(lat1,lon1,lat2,lon2); % deg
[azColor] = azimuthColorWheel(az); % color [c1 c2 c3]
if lat1 == 90
    azColor = [0 0 1]; % Can only go south...
elseif lat1 == -90
    azColor = [1 0 0]; % Can only go north...
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
%%% Plotting the Europa-Centered-Europa-Fixed Frame (ECEF)
% ------------------------------------------------------------------------
if ECEFplot == 1
    figure
    hold all
    axis equal
    trackWidth = 2;
    lineWidth = 3;
    PlotBoi3('X, km','Y, km','Z, km',16)

    %%% Plotting Hopper motion
    plot3(rECEF_Hopper(:,1),rECEF_Hopper(:,2),rECEF_Hopper(:,3),'m','linewidth',trackWidth);
    
%     %%% Plotting Europa equator
%     th = 0:.01:2*pi;
%     x = E_radius * cos(th);
%     y = E_radius * sin(th);
%     plot(x, y,'b','linewidth',lineWidth);

    %%% Coloring in Europa (2D)
    n = 5000;
    THETA=linspace(0,2*pi,n);
    RHO=ones(1,n)*(E_radius);
    [X,Y] = pol2cart(THETA,RHO);
    fill(X, Y, 'c');

    %%% Plotting Europa Surface (3D)
    bodySurface3(E_radius, [0 0 0], [0 1 1]);

    %%% Plotting Europa-Jupiter Vector
    quiver3(0,0,0,-2*E_radius,0,0,...
        'linewidth',2,'color',[1 .5 0]);

    %%% Plotting Intial Velocity
    sc1 = 1000; % scalar
    quiver3(rH01(1),rH01(2),rH01(3),...
        sc1*(vH01(1)),sc1*(vH01(2)),sc1*(vH01(3)),...
        'linewidth',2,'color',[0 0 0])

    xlim([-E_radius*scale1 E_radius*scale1])
    ylim([-E_radius*scale1 E_radius*scale1])
    zlim([-E_radius*scale1 E_radius*scale1])
    
    %%% Frame
    title('ECEF')
end

% ------------------------------------------------------------------------
%%% Plotting the Europa Inertial Frame (ECI)
% ------------------------------------------------------------------------
if ECIplot == 1
    figure
    hold all
    axis equal
    trackWidth = 2;
    lineWidth = 3;
    PlotBoi3('X, km','Y, km','Z, km',16)

    %%% Plotting Hopper motion
    plot3(rECI_Hopper(:,1),rECI_Hopper(:,2),rECI_Hopper(:,3),'m','linewidth',trackWidth);

    %%% Plotting Europa equator
    th = 0:.01:2*pi;
    x = E_radius * cos(th) + rE0(1);
    y = E_radius * sin(th) + rE0(2);
    plot(x, y,'b','linewidth',lineWidth);

    %%% Coloring in Europa (2D)
    n = 5000;
    THETA=linspace(0,2*pi,n);
    RHO=ones(1,n)*(E_radius);
    [X,Y] = pol2cart(THETA,RHO);
    fill(X + rE0(1),Y+ rE0(2),'c');

    %%% Plotting Europa Surface (3D)
    bodySurface3(E_radius, [0 0 0], [0 1 1]);

    %%% Plotting Initial Europa-Jupiter Vector
    quiver3(0,0,0,-2*E_radius,0,0,...
        'linewidth',2,'color',[.7 .7 .7])
    
    %%% Plotting Final Europa-Jupiter Vector
    E2Jhat_f = -r_Europa(end,:)./norm(r_Europa(end,:));
    jpECI = quiver3(0, 0, 0, 2*E_radius*E2Jhat_f(end,1), 2*E_radius*E2Jhat_f(end,2), 2*E_radius*E2Jhat_f(end,3),...
        'linewidth',2,'color',[1 .5 0]);
    legend([jpECI], 'Jupiter-Pointing (Final)')
    
    %%% Plotting Intial Velocity
    sc1 = 1000; % scalar
    quiver3(rH01(1),rH01(2),rH01(3),...
        sc1*(vH02(1)),sc1*(vH02(2)),sc1*(vH02(3)),...
        'linewidth',2,'color',[0 0 0])

    %%% To Focus on Europa
    xlim([-E_radius*scale1 E_radius*scale1])
    ylim([-E_radius*scale1 E_radius*scale1])
    zlim([-E_radius*scale1 E_radius*scale1])
    
    %%% Frame
    title('ECI')
    
end

% ------------------------------------------------------------------------
%%% Plotting the Jupiter System (JCI)
% ------------------------------------------------------------------------
if JCIplot == 1
    figure
    hold all
    grid on
    axis equal
    trackWidth = 2;
    lineWidth = 3;
    PlotBoi3('X, km','Y, km','Z, km',16)

    %%% Plotting Hopper motion
    % p3 = plot3(r_Hopper(:,1),r_Hopper(:,2),r_Hopper(:,3),'m','linewidth',trackWidth,'markersize',5);
    plot3(States(:,1),States(:,2),States(:,3),'m','linewidth',trackWidth,'markersize',5);

    %%% Plotting Europa motion
    plot3(r_Europa(:,1),r_Europa(:,2),r_Europa(:,3),'--r','linewidth',trackWidth);

    %%% Plotting Europa equator
    th = 0:.01:2*pi;
    x = E_radius * cos(th) + rEf(1);
    y = E_radius * sin(th) + rEf(2);
    plot(x, y,'b','linewidth',lineWidth);

    %%% Coloring in Europa (2D)
    n = 5000;
    THETA=linspace(0,2*pi,n);
    RHO=ones(1,n)*(E_radius);
    [X,Y] = pol2cart(THETA,RHO);
    fill(X + rEf(1),Y+ rEf(2),'c');

    %%% Plotting Europa-Jupiter Vector
    EJ_hat = -[rEf(1),rEf(2),rEf(3)]/...
        norm(rEf);
    scale = 2*E_radius;
    quiver3(rEf(1),rEf(2),rEf(3),...
        scale*EJ_hat(1),scale*EJ_hat(2),scale*EJ_hat(3),...
        'linewidth',2,'color',[1 .5 0])

    %%% Plotting Jupiter Equator
    th = 0:.001:2*pi;
    x = J_radius * cos(th);
    y = J_radius * sin(th);
    plot(x, y,'k','linewidth',lineWidth);

    %%% Coloring in Jupiter (2D)
    n = 5000;
    THETA=linspace(0,2*pi,n);
    RHO=ones(1,n)*(J_radius);
    [X,Y] = pol2cart(THETA,RHO);
    fill(X,Y,[1 .5 0]);
    
    %%% Frame
    title('JCI')
end

% ------------------------------------------------------------------------
%%% ECI Distance vs Time plot
% ------------------------------------------------------------------------
if rECIplot == 1
figure
hold all
plot(Times,ones(size(Times)).*E_radius,'--b','linewidth',1.5)
plot(Times,arrayfun(@(x) norm(rECI_Hopper(x,:)), 1:size(rECI_Hopper,1))','m','linewidth',trackWidth)
PlotBoi2('Time, sec','Distance to Europa Center, km',16)
legend('Europa Mean Radius')
end

% ------------------------------------------------------------------------
%%% Movies
% ------------------------------------------------------------------------
%%% ECEF Movie
if runECEFMovie == 1
    figure
    for i = 1:size(States,1) %size(States,1)
        if rem(i,framespeed) == 0 || i == size(States,1)
            clf
            hold all
            %%% Plotting Europa motion
            trackWidth = 2;
            plot3(rECEF_Hopper(i,1),rECEF_Hopper(i,2),rECEF_Hopper(i,3),'rX','linewidth',trackWidth,'markersize',10);
            PlotBoi3('X, km','Y, km','Z, km',16)
            grid on
            axis square
            xlim([-E_radius*scale2 E_radius*scale2])
            ylim([-E_radius*scale2 E_radius*scale2])
            zlim([-E_radius*scale2 E_radius*scale2])

            %%% Plotting Past
            plot3(rECEF_Hopper(1:i-1,1),rECEF_Hopper(1:i-1,2),rECEF_Hopper(1:i-1,3),'-m')

            %%% Europa Equator 
            th = 0:.01:2*pi;
            x = E_radius * cos(th);
            y = E_radius * sin(th);
            p1 = plot(x, y,'b','linewidth',lineWidth);

            %%% Plotting Europa Surface (3D)
            bodySurface3(E_radius, [0 0 0], [0 1 1]);


            %%% Plotting Initial Europa-Jupiter Vector
            quiver3(0,0,0,-2*E_radius,0,0,...
                'linewidth',2,'color',[.7 .7 .7])
            
            %%% View(azimuth,elevation)
%             view(-80,60)
            
            %%% Draw frame
            drawnow limitrate
        end
    end
end

%%% ECI Movie
if runECIMovie == 1
    figure
    for i = 1:size(States,1) %size(States,1)
        if rem(i,framespeed) == 0 || i == size(States,1)
            clf
            hold all
            %%% Plotting Europa motion
            trackWidth = 2;
            plot3(rECI_Hopper(i,1),rECI_Hopper(i,2),rECI_Hopper(i,3),'rX','linewidth',trackWidth,'markersize',10);
            PlotBoi3('X, km','Y, km','Z, km',16)
            grid on
            axis square
            xlim([-E_radius*scale2 E_radius*scale2])
            ylim([-E_radius*scale2 E_radius*scale2])
            zlim([-E_radius*scale2 E_radius*scale2])

            %%% Plotting Past
            plot3(rECI_Hopper(1:i-1,1),rECI_Hopper(1:i-1,2),rECI_Hopper(1:i-1,3),'-m')

            %%% Europa Equator 
            th = 0:.01:2*pi;
            x = E_radius * cos(th);
            y = E_radius * sin(th);
            p1 = plot(x, y,'b','linewidth',lineWidth);

            %%% Plotting Europa Surface (3D)
            bodySurface3(E_radius, [0 0 0], [0 1 1]);


            %%% Plotting Initial Europa-Jupiter Vector
            quiver3(0,0,0,-2*E_radius,0,0,...
                'linewidth',2,'color',[.7 .7 .7])

            %%% Plotting Final Europa-Jupiter Vector
            E2Jhat_f = -r_Europa(i,:)./norm(r_Europa(i,:));
            quiver3(0, 0, 0, 2*E_radius*E2Jhat_f(1), 2*E_radius*E2Jhat_f(2), 2*E_radius*E2Jhat_f(3),...
                'linewidth',2,'color',[1 .5 0]);
            
            %%% View(azimuth,elevation)
            view(5,30)
            
            %%% Draw frame
            drawnow limitrate
        end
    end
end

% ------------------------------------------------------------------------
%%% Printing Impact Time
% ------------------------------------------------------------------------
if Times(end) ~= tf
    fprintf('Time of Flight:    %1.4f sec\n', Times(end))
elseif Times(end) == tf
    warning off backtrace
    warning('No Impact in Simulation')
else
    fprintf('***Impact Time Error***\n\n')
end

% ------------------------------------------------------------------------
%%% Printing Energy Comparison
% ------------------------------------------------------------------------
fprintf('Initial Velocity:  %4.3f m/s\n', v0_m)
fprintf('Final Velocity:    %4.3f m/s\n', vf_m)
fprintf('Delta-V:           %4.3f m/s\n', dVelocity)
fprintf('Delta-E:           %4.3f J/kg\n', dKE)
fprintf('Angle Change:      %3.4f�\n', AngChange)
fprintf('Distance Traveled: %6.3f m\n', Traveled)
fprintf('Azimuth:           %3.2f�\n',az)

% ------------------------------------------------------------------------
%%% Stating Assumptions
% ------------------------------------------------------------------------
fprintf('\n\n\n***************************\n')
fprintf('ASSUMPTIONS:\n')
fprintf('-Europa starts on +X JCI axis (CURRENTLY NECESSARY - I think from lat/lon fnc)\n')
fprintf('-Europa and Jupiter are point masses at centers\n')
fprintf('-Europa is spherical\n')
fprintf('-Europa in circular orbit\n')
fprintf('-NO SMALLER THAN 1 m/s UNLESS SMALLER TIME STEP!\n')
fprintf('-Plenty of others...\n')

