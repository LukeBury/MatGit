clear
clc
close all
addpath('/Users/lukebury/Documents/MATLAB/mbin')
addpath('/Users/lukebury/Documents/MATLAB/CU/bin/LowEnergy')
addpath('/Users/lukebury/Documents/MATLAB/CU/bin')
colors = get_color_palettes();

% ------------------------------------------------------------------------
%%% Got a lot of options
% ------------------------------------------------------------------------
run_MR_Demonstration = 1;
run_J2_Demonstration = 0;


% ------------------------------------------------------------------------
%%% Constants / ICs
% ------------------------------------------------------------------------
% Getting mass ratios
[earth, mars, jupiter, saturn, uranus, neptune, pluto] = getMassRatios(); % e.m, j.e, etc

% Setting colors
% cols = [colors.new.red;colors.new.grn;colors.new.ltblue;colors.new.maglt];
% cols = [38 84 124; 239 71 111; 255 209 102; 6 214 160]./255;
% cols = [6 214 160; 27 154 170; 239 71 111; 255 196 61]./255;
cols = [6 214 160; 38 84 124; 239 71 111; 255 196 61]./255;


%%% Defining cases
sys{1}.u = earth.moon;           % MR, 1.2e-02
sys{1}.a = 384748; % km
sys{1}.w = 2*pi/(27.321661*86400); % rad/s
sys{1}.j2 = 0.001082626;
sys{1}.Prad = 6378; % km
sys{1}.Srad = 1737; % km
sys{1}.Srad_n = sys{1}.Srad/sys{1}.a;
sys{1}.tag = 'Moon, u ~ 1e-02';
sys{1}.moon = 'Moon';
sys{1}.ub1ub2 = [5.972e24, 7.347e22].*6.67408e-20;

sys{2}.u = neptune.triton;       % MR, 2.1e-04
sys{2}.a = 354759;
sys{2}.w = 2*pi/(5.876854*86400);
sys{2}.j2 = 0.004;
sys{2}.Prad = 24622; % km
sys{2}.Srad = 1353;
sys{2}.Srad_n = sys{2}.Srad/sys{2}.a;
sys{2}.tag = 'Triton, u ~ 2e-04';
sys{2}.moon = 'Triton';
sys{2}.ub1ub2 = [1.024e26, 2.14e22].*6.67408e-20;

sys{3}.u = jupiter.europa;       % MR, 2.5e-05
sys{3}.a = 671100;
sys{3}.w = 2*pi/(3.551181*86400);
sys{3}.j2 = 0.01475;
sys{3}.Prad = 69911; % km
sys{3}.Srad = 1560;
sys{3}.Srad_n = sys{3}.Srad/sys{3}.a;
sys{3}.tag = 'Europa, u ~ 3e-05';
sys{3}.moon = 'Europa';
sys{3}.ub1ub2 = [1.898e27,4.799e22].*6.67408e-20;

sys{4}.u = saturn.enceladus;     % MR, 1.9e-07
sys{4}.a = 237948;
sys{4}.w = 2*pi/(1.370218 * 86400);
sys{4}.j2 = 0.01645;
sys{4}.Prad = 58232; % km
sys{4}.Srad = 252;
sys{4}.Srad_n = sys{4}.Srad/sys{4}.a;
sys{4}.tag = 'Enceladus, u ~ 2e-07';
sys{4}.moon = 'Enceladus';
sys{4}.ub1ub2 = [5.683e26,1.08e20].*6.67408e-20;

% ------------------------------------------------------------------------
%%% Jump from different mass ratios
% ------------------------------------------------------------------------
% % ---------------------------
% %%% Plotting normalized secondary body
% % ---------------------------
% figure; hold all
% axis equal
% p = {};
% 
% for kk = 1:length(sys)
%     % Normalizing constants ... (real/norm = normalized)
%     rNorm = sys{kk}.a;
%     tNorm = 1/sys{kk}.w;
%     vNorm = rNorm/tNorm;
%     
%     %%% Rotating frame cooridinates
%     rP_BCR_n = [-sys{kk}.u, 0, 0];
%     rS_BCR_n = [1-sys{kk}.u, 0, 0];
% 
%     % ---------------------------
%     %%% Defining Particle State
%     % ---------------------------
%     %%% Initial Particle Position
%     startPos_n = [sys{kk}.Srad_n, 0, 0];
%     rsc_BCR_n = rS_BCR_n + startPos_n;
% 
%     %%% Initial Partical Velocity
%     vsc_BCR_n = [1.5,.4,0]./vNorm;
%     
%     % ---------------------------
%     %%% Propagating State
%     % ---------------------------
%     %%% Setting normalized time vector (seconds/tNorm)
%     ti = 0./tNorm;
%     dt = 1./tNorm;
% %     tf = 1.15*3600./tNorm;
%     tf = 2.5*3600./tNorm;
%     time = ti:dt:tf;
% 
%     %%% Choosing ode45 tolerance
%     tol = 1e-13;
% 
%     %%% Setting integrator options
%     options = odeset('Events',@normalCR3BP_impactEvent,'RelTol',tol,'AbsTol',tol);
% 
%     %%% Setting Initial State Vector (ECEF)
%     X0_n = [rsc_BCR_n, vsc_BCR_n]'; % km, km/s
% 
%     %%% Propagating the State
%     [Times_n,States_BCR_n] = ode45(@normalCR3BP_Int,time,X0_n,options,sys{kk}.u,rP_BCR_n,rS_BCR_n,sys{kk}.Srad_n);
%     
%     % Turning positions real
%     position_BCR = States_BCR_n(:,1:3).*rNorm;
%     position_BCR = position_BCR - position_BCR(1,1:3);
%     position_BCR = position_BCR + startPos_n*rNorm;
%     
%     % ---------------------------
%     %%% Plotting to normalized secondary 
%     % ---------------------------
%     p{kk} = plot(position_BCR(:,1),position_BCR(:,2),'linewidth',2,'color',cols(kk,:));
%     
%     % Plotting secondary surface
%     x = sys{kk}.Srad * cos(0:.001:2*pi);
%     y = sys{kk}.Srad * sin(0:.001:2*pi);
%     plot(x, y, 'color', cols(kk,:), 'linewidth', 1.5);
% 
%     
% 
% 
% end
% % Plot format
% PlotBoi2('X, km','Y, km',16)
% xlim([-0.1487,4.5030]*1.0e+03)
% ylim([-1.8261,1.8427]*1.0e+03)
% % Legend
% legend([p{1},p{2},p{3}, p{4}],sys{1}.tag,sys{2}.tag,sys{3}.tag,sys{4}.tag)

% ------------------------------------------------------------------------
%%% L1/L2 of different mass ratios
% ------------------------------------------------------------------------
% figure; hold all
% axis equal
for kk = 1:length(sys)
    % Acquire L1/L2 points and center at 0
    L = EquilibriumPoints(sys{kk}.u);
    sys{kk}.L1L2 = [L(1,1), L(2,1)] - (1-sys{kk}.u);
    
%     % Plot L1/L2 points
%     p{kk} = plot(sys{kk}.L1L2, [0,0], 'x','markersize',7,'linewidth',1.5,'color',cols(kk,:));
%     
%     % Plotting secondary surface
%     x = sys{kk}.Srad_n * cos(0:.001:2*pi);
%     y = sys{kk}.Srad_n * sin(0:.001:2*pi);
%     plot(x, y, 'color', cols(kk,:), 'linewidth', 1.5);
end
% % Plot format
% PlotBoi2('X','Y',16)
% title('Moons and L1/L2 points normalized to CR3BP')
% % Legend
% legend([p{1},p{2},p{3}, p{4}],sys{1}.tag,sys{2}.tag,sys{3}.tag,sys{4}.tag)

% ------------------------------------------------------------------------
%%% L1/L2 radii away from surface
% ------------------------------------------------------------------------
for kk = 1:length(sys)
    sys{kk}.L1radiiFromSurface = -sys{kk}.L1L2(1) / sys{kk}.Srad_n - 1;
end

fprintf('-------------------- L1 radii from surface \n')
fprintf('Earth-Moon\t\t%1.2e\t%2.2f\n',sys{1}.u,sys{1}.L1radiiFromSurface)
fprintf('Neptune-Triton\t\t%1.2e\t%2.2f\n',sys{2}.u,sys{2}.L1radiiFromSurface)
fprintf('Jupiter-Europa\t\t%1.2e\t%2.2f\n',sys{3}.u,sys{3}.L1radiiFromSurface)
fprintf('Saturn-Enceladus\t%1.2e\t%2.2f\n',sys{4}.u,sys{4}.L1radiiFromSurface)


% ------------------------------------------------------------------------
%%% Playing w/ enceladus
% ------------------------------------------------------------------------
if run_MR_Demonstration == 1
% % % % velocities = [1.8, .5, .4, .05];
% % % % % velocities = [1.1, .5, .4, .05];

% For const normalized velocity (also need to change velocity)
simTimes = [660,137,79,18];
ylimMult = [0.63,0.63,0.7,0.7];

% % For constant real velocity
% simTimes = [550,130,85,18];
figure
    
for ss = 1:length(sys)
    subplot(2,2,ss); 
    hold all; axis equal
    angles = [0:5:170].*(pi/180);
    sys{ss}.crashCount = 0;

    for kk = 1:length(angles)
        % Normalizing constants ... (real/norm = normalized)
        rNorm = sys{ss}.a;
        tNorm = 1/sys{ss}.w;
        vNorm = rNorm/tNorm;

        %%% Rotating frame cooridinates
        rP_BCR_n = [-sys{ss}.u, 0, 0];
        rS_BCR_n = [1-sys{ss}.u, 0, 0];

        % ---------------------------
        %%% Defining Particle State
        % ---------------------------
        %%% Initial Particle Position
    %     startPos_n = [sys{kk}.Srad_n, 0, 0];
        startPos_n = [sys{ss}.L1L2(2), 0, 0];
        rsc_BCR_n = rS_BCR_n + startPos_n;

        %%% Initial Partical Velocity
% % %         vsc_BCR_n = R3([0,velocities(ss),0],angles(kk))./vNorm;
        vsc_BCR_n = R3([0,0.0040,0],angles(kk)); % const normalized velocity (also need to change time)
%         vsc_BCR_n = R3([0,.05,0],angles(kk))./vNorm; % const real velocity

        % ---------------------------
        %%% Propagating State
        % ---------------------------
        %%% Setting normalized time vector (seconds/tNorm)
        ti = 0./tNorm;
        dt = 10*60./tNorm;
        tf = simTimes(ss)*3600./tNorm;
        time = ti:dt:tf;

        %%% Choosing ode45 tolerance
        tol = 1e-13;

        %%% Setting integrator options
        options = odeset('Events',@normalCR3BP_impactEvent,'RelTol',tol,'AbsTol',tol);

        %%% Setting Initial State Vector (ECEF)
        X0_n = [rsc_BCR_n, vsc_BCR_n]'; % km, km/s

        %%% Propagating the State
        [Times_n,States_BCR_n] = ode45(@normalCR3BP_Int,time,X0_n,options,sys{ss}.u,rP_BCR_n,rS_BCR_n,sys{ss}.Srad_n);
        
        % Checking if trajectory landed
        if Times_n(end) ~= tf
            sys{ss}.crashCount = sys{ss}.crashCount + 1;
        end
%         % Turning positions real
%         position_BCR = States_BCR_n(:,1:3).*rNorm;
%         position_BCR = position_BCR - position_BCR(1,1:3);
%         position_BCR = position_BCR + startPos_n*rNorm;
        
        % Keeping positions normalized
        position_BCR = States_BCR_n(:,1:3);
%         position_BCR = position_BCR - position_BCR(1,1:3);
%         position_BCR = position_BCR + startPos_n;

        % ---------------------------
        %%% Plotting to normalized secondary 
        % ---------------------------
        p{4} = plot(position_BCR(:,1),position_BCR(:,2),'linewidth',1,'color',cols(ss,:));

        % Plotting secondary surface
        x = sys{ss}.Srad_n * cos(0:.001:2*pi) + rS_BCR_n(1);
        y = sys{ss}.Srad_n * sin(0:.001:2*pi);
        plot(x, y, 'color', 'k', 'linewidth', 1.5);

        % Plotting L1/L2
        plot((sys{ss}.L1L2)+rS_BCR_n(1), [0,0], 'k+','markersize',7,'linewidth',2);

        
    end
    % Plot Formatting
    if ss == 1
        PlotBoi2('','y',12)
    elseif ss == 3
        PlotBoi2('x','y',12)
    elseif ss == 4
        PlotBoi2('x','',12)
    end
    
    xlim(([sys{ss}.L1L2(1),sys{ss}.L1L2(2)].*1.3)+rS_BCR_n(1))
    ylim([sys{ss}.L1L2(1),sys{ss}.L1L2(2)].*ylimMult(ss))
    title(sys{ss}.tag)
end

for ss = 1:length(sys)
    sys{ss}.crashPerc = sys{ss}.crashCount / length(angles);
end


fprintf('\n\n-------------------- Crash percentage \n')
fprintf('Earth-Moon\t\t%2.0f\n',sys{1}.crashPerc * 100)
fprintf('Neptune-Triton\t\t%2.0f\n',sys{2}.crashPerc * 100)
fprintf('Jupiter-Europa\t\t%2.0f\n',sys{3}.crashPerc * 100)
fprintf('Saturn-Enceladus\t%2.0f\n',sys{4}.crashPerc * 100)

end
% ========================================================================
%%% J2 Demonstration
% ------------------------------------------------------------------------
if run_J2_Demonstration == 1
% % % syms rsb2 rb2b1 Rsb1 ub1 ub2 Rb1 J2
% % % 
% % % U_2B = ub2/rsb2;
% % % 
% % % U_3B = ub1/rb2b1 - ub1/rsb1;
% % % 
% % % % U_J2 = -(3*u*RE*RE*J2*z*z)/(2*(r^5)) + (u*RE*RE*J2)/(2*(r^3))
% % % U_J2_noZ = (ub1*Rb1*Rb1*J2)/(2*(r^3));
fprintf('\n\n')
for ss = 4
    ub1 = sys{ss}.ub1ub2(1);
    ub2 = sys{ss}.ub1ub2(2);
    J2 = sys{ss}.j2;
    Rb1 = sys{ss}.Prad;
    
    % Looking at L1 position
    rb2b1 = sys{ss}.a; % km
    rsb2 = -sys{ss}.L1L2(1) * sys{ss}.a; % km, L2 wrt B2 un-normalized
    rsb1 = rb2b1 + rsb2; % km

    % Calculating potentials
    U_2B = ub2/rsb2;
    U_3B = ub1*(1/rb2b1 - 1/rsb1);
    U_J2_noZ = (ub1*Rb1*Rb1*J2)/(2*(rsb1^3));
    
    % Calculating difference made by J2
    diff = U_J2_noZ / (U_2B + U_3B) *100;
    fprintf('B1 J2 effect on L1 = %2.6f%%\n',diff)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% negative
    % Calculating potentials
    U_2B = ub2/rsb2;
    U_3B = -ub1*(1/rb2b1 - 1/rsb1);
    U_J2_noZ = (ub1*Rb1*Rb1*J2)/(2*(rsb1^3));
    
    % Calculating difference made by J2
    diff2 = U_J2_noZ / (U_2B + U_3B) *100;
    fprintf('222: B1 J2 effect on L1 = %2.6f%%\n',diff2)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% accelerations
    % Accelertion at L1
    a_2B = -ub2/(rsb2^2);
    a_3B = -ub1 * (1/(rsb1^2) - 1/(rb2b1^2));
    a_J2 = -(3/2)*ub1*Rb1*Rb1*J2/(rsb1^4);
    
    diff3 = a_J2/(a_2B + a_3B);
    fprintf('ACC: B1 J2 effect on L1 = %2.6f%%\n\n',diff3)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normalized potential (L1) [koon p8]
    u = sys{ss}.u;
    r1 = 1 - sys{ss}.L1L2(1);
    r2 = sys{ss}.L1L2(1);
    x = (1-u) - r2;
    y = 0;
    rNorm = sys{ss}.a;
    tNorm = 1/sys{ss}.w;
    % Calculating potentials (normalized -> unnormalized)
    U_xy_n = (1-u)/r1 + u/r2 + (x^2 + y^2)/2;
%     U_xy = U_xy_n * (rNorm^2)/(tNorm^2);
    U_J2_noZ = (ub1*Rb1*Rb1*J2)/(2*(rsb1^3));
    U_J2_noZ_n = U_J2_noZ * (tNorm^2)/(rNorm^2);
    
    % Calculating difference made by J2
    diffn = U_J2_noZ_n / (U_xy_n) *100;
%     U_J2_noZ*100/U_xy_n;
    fprintf('nnn: B1 J2 effect on L1 = %2.6f%%\n',diffn)
    
    
    
%     % Looking at L2 position
%     rb2b1 = sys{ss}.a; % km
%     rsb2 = sys{ss}.L1L2(2) * sys{ss}.a; % km, L2 wrt B2 un-normalized
%     rsb1 = rb2b1 + rsb2; % km
% 
%     % Calculating potentials
%     U_2B = ub2/rsb2;
%     U_3B = ub1*(1/rb2b1 - 1/rsb1);
%     U_J2_noZ = (ub1*Rb1*Rb1*J2)/(2*(rsb1^3));
%     
%     % Calculating difference made by J2
%     diff = U_J2_noZ / (U_2B + U_3B) *100;
%     fprintf('B1 J2 effect on L2 = %2.6f%%\n\n',diff)

    
end


% U(x,y) = -(1-u)/r1 - u/r2
% where r1 = dist to primary and r2 = dist to secondary






% % Utot_uJ2J3 = -u/r + U_J2J3; % Total potential
% % 
% % EQM = [dx; dy; dz; diff(Utot_uJ2J3, x); diff(Utot_uJ2J3, y); diff(Utot_uJ2J3,z)];

end
