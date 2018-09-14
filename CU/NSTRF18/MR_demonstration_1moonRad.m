clear
clc
close all
addpath('/Users/lukebury/Documents/MATLAB/mbin')
colors = get_color_palettes();

% ------------------------------------------------------------------------
%%% Constants / ICs
% ------------------------------------------------------------------------
% Getting mass ratios
[earth, mars, jupiter, saturn, uranus, neptune, pluto] = getMassRatios(); % e.m, j.e, etc

%%% Defining cases
sys{1}.u = earth.moon;           % MR, 1.2e-02
sys{1}.a = 384748;
sys{1}.w = 2*pi/(27.321661*86400); % rad/s
sys{1}.Srad = 1737;
sys{1}.Srad_n = sys{1}.Srad/sys{1}.a;
sys{1}.tag = 'u = 1.2e-02';

sys{2} = sys{1}; sys{3} = sys{1}; sys{4} = sys{1}; sys{5} = sys{1}; sys{6} = sys{1};

sys{2}.u = 1.2e-03;
sys{2}.tag = '1.2e-03';
sys{3}.u = 1.2e-05;
sys{3}.tag = '1.2e-05';
sys{4}.u = 1.2e-07;
sys{4}.tag = '1.2e-07';




% ------------------------------------------------------------------------
%%% Setting up simulation
% ------------------------------------------------------------------------
% ---------------------------
%%% Plotting normalized secondary body
% ---------------------------
figure; hold all
axis equal
p = {};
for kk = 1:length(sys)
    % Normalizing constants ... (real/norm = normalized)
    rNorm = sys{kk}.a;
    tNorm = 1/sys{kk}.w;
    vNorm = rNorm/tNorm;
    
    %%% Rotating frame cooridinates
    rP_BCR_n = [-sys{kk}.u, 0, 0];
    rS_BCR_n = [1-sys{kk}.u, 0, 0];

    % ---------------------------
    %%% Defining Particle State
    % ---------------------------
    %%% Initial Particle Position
    startPos_n = [sys{kk}.Srad_n, 0, 0];
    rsc_BCR_n = rS_BCR_n + startPos_n;

    %%% Initial Partical Velocity
    vsc_BCR_n = [1,.05,0]./vNorm;
    
    % ---------------------------
    %%% Propagating State
    % ---------------------------
    %%% Setting normalized time vector (seconds/tNorm)
    ti = 0./tNorm;
    dt = 1./tNorm;
    tf = 1.15*3600./tNorm;
    time = ti:dt:tf;

    %%% Choosing ode45 tolerance
    tol = 1e-13;

    %%% Setting integrator options
    options = odeset('Events',@normalCR3BP_impactEvent,'RelTol',tol,'AbsTol',tol);

    %%% Setting Initial State Vector (ECEF)
    X0_n = [rsc_BCR_n, vsc_BCR_n]'; % km, km/s

    %%% Propagating the State
    [Times_n,States_BCR_n] = ode45(@normalCR3BP_Int,time,X0_n,options,sys{kk}.u,rP_BCR_n,rS_BCR_n,sys{kk}.Srad_n);
    
    % Turning positions real
    position_BCR = States_BCR_n(:,1:3).*rNorm;
    position_BCR = position_BCR - position_BCR(1,1:3);
    position_BCR = position_BCR + startPos_n*rNorm;
    
    % ---------------------------
    %%% Plotting to normalized secondary 
    % ---------------------------
    p{kk} = plot(position_BCR(:,1),position_BCR(:,2),'linewidth',2);
    kk
    % Plotting secondary surface
    x = sys{kk}.Srad * cos(0:.001:2*pi);
    y = sys{kk}.Srad * sin(0:.001:2*pi);
    plot(x, y, 'linewidth', 1.5);

    PlotBoi2('X','Y',16)


end

% % Plotting secondary surface
% x = cos(0:.001:2*pi);
% y = sin(0:.001:2*pi);
% plot(x, y, 'color', colors.new.blue, 'linewidth', 1.5);

% Legend
legend([p{1},p{2},p{3}, p{4}],sys{1}.tag,sys{2}.tag,sys{3}.tag,sys{4}.tag)