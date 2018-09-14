clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

%% =======================================================================
%%% Setting Up
% ========================================================================
% -----------------------------------------------------------------
%%% Run Switches
% -----------------------------------------------------------------
run_symAnalysis        = 0; % 3.1
run_initialSim         = 0; % 3.2 -should be run together
run_geoModeling        = 0; % 3.3 -should be run together
run_controlSim         = 0; % 3.4
run_inclinationSim     = 0; % 3.5
run_monteCarloAnalysis = 1; % 3.6
% -----------------------------------------------------------------
%%% General Data
% -----------------------------------------------------------------
%%% Color options/schemes
colors = get_colors();

%%% Figures folder path
figPath = '/Users/lukebury/Documents/School/CU/Courses/4-5010-Attitude_Dynamics/Project/Project_Latex/Figures/';
% -----------------------------------------------------------------
%%% Relating to Earth
% -----------------------------------------------------------------
RE = 6378; % km, radius of earth
uE = 398600; % km^2 s^-2, grav parameter of earth
wE = 7.2921159e-5; % rad/s, rotation rate of earth about n3

% -----------------------------------------------------------------
%%% Relating to satellite motion
% -----------------------------------------------------------------
%%% Position
altitude = 450; % km, altitude of s/c
r = RE + altitude; % km, position magnitude of s/c

%%% Mean Motion
th_dot = sqrt(uE/(r^3)); % rad/s, angular rate of s/c, constant
Tp = 2*pi/th_dot;
% -----------------------------------------------------------------
%%% Relating to Dipole
% -----------------------------------------------------------------
M = 7.838e6; % T km^3, dipole moment
gm = 17*pi/180; % rad, tilt angle of dipole 
Bm0 = 0; % rad, magnetic field angle

% -----------------------------------------------------------------
%%% Initial states
% -----------------------------------------------------------------
sig0 = [0.3; 0.2; 0.4]; % MRP
w0_BN = [15; 8; 12].*pi/180; % rad/s, initial angular velocity


I1 = 3.5; % kg m^2
I2 = 5;   % kg m^2
I3 = 8;   % kg m^2
I = blkdiag(I1, I2, I3); % Inertia tensor of s/c in body frame,  kg m^2

%%% Euler Angles
raan0 = 0; % rad
in0 = 45 * pi/180; % rad
th0 = 0; % rad, true anomaly

%% =======================================================================
%%% Problems - Midway Report
% ========================================================================
% -----------------------------------------------------------------
%%% 3.1 - Rotation Matrices
% -----------------------------------------------------------------
if run_symAnalysis == 1
syms v i raan gamma K real

HN = R3(eye(3),-v)*R1(eye(3),-i)*R3(eye(3),-raan0)

MN = R3(eye(3),-0)*R1(eye(3),-gamma)*R3(eye(3),-K)

HM = HN * MN'

return
end
% -----------------------------------------------------------------
%%% 3.2 - Integration via RK4
% -----------------------------------------------------------------
if run_initialSim == 1
%%% Initial state
X0 = [sig0;w0_BN];

%%% Time vector
t0 = 0; % sec
dt = .1; % sec
tf = Tp; % sec
times = t0:dt:tf;

%%% Preallocating
states = zeros(6,length(times));
controls = zeros(3,length(times));
Hs = zeros(1,length(times));
Ts = zeros(1,length(times));

%%% Storing ICs
states(:,1) = X0;
Ts(1) = (1/2) * states(4:6,1)' * I * states(4:6,1);
H = I * states(4:6,1);
Hs(1) = sqrt(H'*H); % taking norm

%%% Integrating and calculating H and T
for kk = 1:length(times)-1    
    % ---------------------
    % Integrating states
    % ---------------------
    %%% Checking MRPs for shadow set switch
    s = norm(states(1:3,kk));
    if s > 1
        states(1:3,kk) = -states(1:3,kk)./(s^2);
    end
    
    [X_p] = rk4_MRPw(states(:,kk),dt,I,controls(:,kk));
    states(:,kk+1) = X_p;
    
    % ---------------------
    % Computing T and H
    % ---------------------
    Ts(kk+1) = (1/2) * states(4:6,kk+1)' * I * states(4:6,kk+1);
    H = I * states(4:6,kk+1);
    Hs(kk+1) = sqrt(H'*H); % taking norm
end

% ---------------------
% Plotting MRPs
% ---------------------
figure; hold all
plot(times(1:1000)./60,states(1,1:1000),'.','linewidth',2)
plot(times(1:1000)./60,states(2,1:1000),'.','linewidth',2)
plot(times(1:1000)./60,states(3,1:1000),'.','linewidth',2)
PlotBoi2('Time, min','MRP Value',14)
[legh,objh] = legend('\sigma1', '\sigma2', '\sigma3');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-');
ylim([-1 1])
print(gcf,'-dpng',[figPath,'MRPs.png']);

% ---------------------
% Plotting Angular Velocity
% ---------------------
figure; hold all
plot(times(1:1000)./60,states(4,1:1000),'.','linewidth',3)
plot(times(1:1000)./60,states(5,1:1000),'.','linewidth',3)
plot(times(1:1000)./60,states(6,1:1000),'.','linewidth',3)
PlotBoi2('Time, min','Angular Velocity, rad/s',14)
[legh,objh] = legend('\omega1', '\omega2', '\omega3');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-');
print(gcf,'-dpng',[figPath,'AngVelocities.png']);

% ---------------------
% Plotting Angular Momentun
% ---------------------
figure
plot(times(1:1000)./60,percentchange(Hs(1:1000)),'linewidth',1,'color',colors.std.purp)
PlotBoi2('Time, min','H %-Change',14)
print(gcf,'-dpng',[figPath,'AngMom.png']);

% ---------------------
% Plotting Rotational Energy
% ---------------------
figure
plot(times(1:1000)./60,percentchange(Ts(1:1000)),'linewidth',1,'color',colors.std.purp)
PlotBoi2('Time, min','T %-Change',14)
print(gcf,'-dpng',[figPath,'RotEnergy.png']);

end % run_initialSim

% -----------------------------------------------------------------
%%% 3.3 - Geomagnetic Modeling
% -----------------------------------------------------------------
if run_geoModeling == 1
%%% Preallocating
bs_H = zeros(3,length(times));
bs_B = zeros(3,length(times));

%%% Looping through times
for kk = 1:length(times)
    %%% Finding current ellapsed time
    ti = times(kk); % sec
    
    %%% Calculate Bm
    Bm = Bm0 + wE*ti; % rad
    
    %%% Calculate b_H
    zeta = acos(cos(in0)*cos(gm) + sin(in0)*sin(gm)*cos(raan0-Bm));
    eta = asin(sin(gm)*sin(raan0-Bm)/sin(zeta)); 
    b_H = (M/r^3)*[cos(th_dot*ti-eta)*sin(zeta); cos(zeta); -2*sin(th_dot*ti-eta)*sin(zeta)];
    

    bs_H(:,kk) = b_H;
    
    %%% Calculate [BH]
    [ BN ] = mrp2DCM( states(1:3,kk) );
    HN = [cos(raan0)*cos(th_dot*ti) - cos(in0)*sin(raan0)*sin(th_dot*ti), cos(th_dot*ti)*sin(raan0) + cos(in0)*cos(raan0)*sin(th_dot*ti), sin(in0)*sin(th_dot*ti);...
          -cos(raan0)*sin(th_dot*ti) - cos(in0)*cos(th_dot*ti)*sin(raan0), cos(in0)*cos(raan0)*cos(th_dot*ti) - sin(raan0)*sin(th_dot*ti), cos(th_dot*ti)*sin(in0);...
          sin(in0)*sin(raan0), -cos(raan0)*sin(in0), cos(in0)];
    BH = BN*HN';
    
    %%% Calculate b_B
    b_B = BH*b_H;
    bs_B(:,kk) = b_B;
    
end

figure
subplot(3,1,1)
plot(times./60,bs_H(1,:),'linewidth',2,'color',colors.std.purp)
PlotBoi2('','$$^H\hat{b}_1$$',14,'Latex')
subplot(3,1,2)
plot(times./60,bs_H(2,:),'linewidth',2,'color',colors.std.purp)
PlotBoi2('','$$^H\hat{b}_2$$',14,'Latex')
subplot(3,1,3)
plot(times./60,bs_H(3,:),'linewidth',2,'color',colors.std.purp)
PlotBoi2('Time, min','$$^H\hat{b}_3$$',14,'Latex')
print(gcf,'-dpng',[figPath,'b_H.png']);

figure
subplot(3,1,1)
plot(times(1:1000)./60,bs_B(1,1:1000),'linewidth',2,'color',colors.std.purp)
PlotBoi2('','$$^B\hat{b}_1$$',14,'Latex')
subplot(3,1,2)
plot(times(1:1000)./60,bs_B(2,1:1000),'linewidth',2,'color',colors.std.purp)
PlotBoi2('','$$^B\hat{b}_2$$',14,'Latex')
subplot(3,1,3)
plot(times(1:1000)./60,bs_B(3,1:1000),'linewidth',2,'color',colors.std.purp)
PlotBoi2('Time, min','$$^B\hat{b}_3$$',14,'Latex')
print(gcf,'-dpng',[figPath,'b_B.png']);

end % run_geoModeling

% -----------------------------------------------------------------
%%% 3.4 - Control Law Implementation
% -----------------------------------------------------------------
if run_controlSim == 1
%%% Time vector
t0 = 0; % sec
dt = .1; % sec

%%% Titles for cases
controlTitles = {'Modulating B-dot Control Law','Bang-Bang B-dot Control Law'};
%%% Max control
mMax = 3; % Am^2

%%% Setting initial conditions
X0 = {[0.3; 0.2; 0.4; [15;8;12].*pi/180]...
    [0.1; 0.1; 0.4; [1;12;1].*pi/180]...
    [0.35; 0.2; 0.15; [6;4;13].*pi/180]};

%%% Setting integration times so that 3 deg/s is met
tfs = [3.5, 2.5, 4.5];
989
989
for controlLaw = 1:1 % 1:2
    for ICi = 1:1 % 1:3
        %%% Setting time vector
        tf = Tp*tfs(ICi); % sec
        times = t0:dt:tf;
        
        %%% Preallocating
        bs_H = zeros(3,length(times));
        bs_B = zeros(3,length(times));
        states = zeros(6,length(times));
        controls = zeros(3,length(times));
        commandDipoles = zeros(3,length(times));
        
        %%% Storing ICs
        states(:,1) = X0{ICi};

        for kk = 1:length(times)-1
            ti = times(kk); % sec

            % ---------------------
            % Magnetic Field Info
            % ---------------------
            %%% Calculate Bm
            Bm = Bm0 + wE*ti; % rad

            %%% Calculate b_H
            zeta = acos(cos(in0)*cos(gm) + sin(in0)*sin(gm)*cos(raan0-Bm));
            eta = asin(sin(gm)*sin(raan0-Bm)/sin(zeta)); 
            b_H = (M/r^3)*[cos(th_dot*ti-eta)*sin(zeta); cos(zeta); -2*sin(th_dot*ti-eta)*sin(zeta)];

            bs_H(:,kk) = b_H;

            %%% Calculate [BH]
            [ BN ] = mrp2DCM( states(1:3,kk) );
            HN = [cos(raan0)*cos(th_dot*ti) - cos(in0)*sin(raan0)*sin(th_dot*ti), cos(th_dot*ti)*sin(raan0) + cos(in0)*cos(raan0)*sin(th_dot*ti), sin(in0)*sin(th_dot*ti);...
                  -cos(raan0)*sin(th_dot*ti) - cos(in0)*cos(th_dot*ti)*sin(raan0), cos(in0)*cos(raan0)*cos(th_dot*ti) - sin(raan0)*sin(th_dot*ti), cos(th_dot*ti)*sin(in0);...
                  sin(in0)*sin(raan0), -cos(raan0)*sin(in0), cos(in0)];
            BH = BN*HN';

            %%% Calculate b_B
            b_B = BH*b_H;
            bs_B(:,kk) = b_B;
            
            % ---------------------
            % Control
            % ---------------------
            %%% Set gain
            gain = 1;
            % c1_IC1_g50.0 - Threshold crossing: 296.81 min
            % c1_IC1_g1.0 - No threshold crossing
            
            %%% Calculate control
            if controlLaw == 1
                [u,m] = control_modulatingBdot(th_dot,zeta,I1,bs_B(:,kk),states(4:6,kk),mMax,gain);
            elseif controlLaw == 2
                if kk > 1
                    [u,m] = control_bangBang(mMax,bs_B(:,kk-1),bs_B(:,kk),dt);
                else
                    u = zeros(3,1);
                    m = zeros(3,1);
                end
            end
            
            controls(:,kk) = u;
            commandDipoles(:,kk) = m;
            % ---------------------
            % Integrating states
            % ---------------------
            %%% Checking MRPs for shadow set switch
            s = norm(states(1:3,kk));
            if s > 1
                states(1:3,kk) = -states(1:3,kk)./(s^2);
            end

            [X_p] = rk4_MRPw(states(:,kk),dt,I,controls(:,kk));
            states(:,kk+1) = X_p;

        end % times loop
        
        %%% Finding index where error threshold is crossed
        controlError = colnorm(states(4:6,:));
        errorIndex = find(controlError<=(3*pi/180));
        if isempty(errorIndex) == 0
            errorIndex = errorIndex(1);
            fprintf('c%1.0f_IC%1.0f_g%1.1f - Threshold crossing: %2.2f min\n',controlLaw,ICi,gain*10,times(errorIndex)/60)
        else
            fprintf('c%1.0f_IC%1.0f_g%1.1f - No threshold crossing\n',controlLaw,ICi,gain*10)
        end
        % ---------------------
        % Plotting case
        % ---------------------
        %%% Plotting MRPs
        figure
        subplot(3,1,1)
        plot(times./Tp,states(1,:),'.','linewidth',0.75,'color',colors.std.purp)
        title(sprintf('%s  ...  IC #%1.0f  ...  Gain: %1.1f',controlTitles{controlLaw},ICi,gain));
        ylim([-1 1])
        PlotBoi2('','\sigma_1',18)
        subplot(3,1,2)
        plot(times./Tp,states(2,:),'.','linewidth',0.75,'color',colors.std.purp)
        ylim([-1 1])
        PlotBoi2('','\sigma_2',18)
        subplot(3,1,3)
        plot(times./Tp,states(3,:),'.','linewidth',0.75,'color',colors.std.purp)
        ylim([-1 1])
        PlotBoi2('Time, period','\sigma_3',18)
        print(gcf,'-dpng',[figPath,'MRPs',sprintf('_c%1.0f_x%1.0f_g%1.0f',controlLaw,ICi,gain*10),'.png']);
%         figure; hold all
%         plot(times./Tp,states(1,:),'.','linewidth',1)
%         plot(times./Tp,states(2,:),'.','linewidth',1)
%         plot(times./Tp,states(3,:),'.','linewidth',1)
%         PlotBoi2('Time, period','MRP Value',18)
%         title(sprintf('%s  ...  IC #%1.0f  ...  Gain: %1.1f',controlTitles{controlLaw},ICi,gain));
%         [legh,objh] = legend('\sigma1', '\sigma2', '\sigma3');
%         lineh = findobj(objh,'type','line');
%         set(lineh,'linestyle','-');
%         ylim([-1 1])
%         print(gcf,'-dpng',[figPath,'MRPs',sprintf('_c%1.0f_x%1.0f_g%2.0f',controlLaw,ICi,gain*10),'.png']);

        %%% Plotting Angular Velocity
        figure; hold all
        plot(times./Tp,states(4,:),'linewidth',1.5)
        plot(times./Tp,states(5,:),'linewidth',1.5)
        plot(times./Tp,states(6,:),'linewidth',1.5)
        plot(times./Tp,controlError,'linewidth',1.5,'color','k')
        plot([times(1),times(end)]./Tp,[3*pi/180, 3*pi/180],'--r','linewidth',2)
        if isempty(errorIndex) == 0
            plot([times(errorIndex),times(errorIndex)]./Tp,[-max(controlError), max(controlError)],'--m','linewidth',2)
        end
        plot([times(1),times(end)]./Tp,[-3*pi/180, -3*pi/180],'--r','linewidth',2)
        PlotBoi2('Time, period','Angular Velocity, rad/s',18)
        title(sprintf('%s  ...  IC #%1.0f  ...  Gain: %1.1f',controlTitles{controlLaw},ICi,gain));
        [legh,objh] = legend('\omega1', '\omega2', '\omega3','|\omega| - control error','3 °/s threshold','threshold crossing');
%         lineh = findobj(objh,'type','line');
%         set(lineh,'linestyle','-');
        print(gcf,'-dpng',[figPath,'AngVelocities',sprintf('_c%1.0f_x%1.0f_g%1.0f',controlLaw,ICi,gain*10),'.png']);
        
        %%% Plotting control torque
        figure
        subplot(3,1,1)
        plot(times./Tp,controls(1,:),'linewidth',0.7,'color',colors.std.purp)
        PlotBoi2('','u$$_x$$, Nm',18,'Latex')
        title(sprintf('%s  ...  IC #%1.0f  ...  Gain: %1.1f',controlTitles{controlLaw},ICi,gain));
        subplot(3,1,2)
        plot(times./Tp,controls(2,:),'linewidth',0.7,'color',colors.std.purp)
        PlotBoi2('','u$$_y$$, Nm',18,'Latex')
        subplot(3,1,3)
        plot(times./Tp,controls(3,:),'linewidth',0.7,'color',colors.std.purp)
        PlotBoi2('Time, period','u$$_z$$, Nm',18,'Latex')
        print(gcf,'-dpng',[figPath,'controlTorque',sprintf('_c%1.0f_x%1.0f_g%1.0f',controlLaw,ICi,gain*10),'.png']);
        
        %%% Plotting command dipole
        figure
        subplot(3,1,1)
        plot(times./Tp,commandDipoles(1,:),'linewidth',0.7,'color',colors.std.purp)
        PlotBoi2('','m$$_x$$, Am$$^2$$',18,'Latex')
        title(sprintf('%s  ...  IC #%1.0f  ...  Gain: %1.1f',controlTitles{controlLaw},ICi,gain));
        subplot(3,1,2)
        plot(times./Tp,commandDipoles(2,:),'linewidth',0.7,'color',colors.std.purp)
        PlotBoi2('','m$$_y$$, Am$$^2$$',18,'Latex')
        subplot(3,1,3)
        plot(times./Tp,commandDipoles(3,:),'linewidth',0.7,'color',colors.std.purp)
        PlotBoi2('Time, period','m$$_z$$, Am$$^2$$',18,'Latex')
        print(gcf,'-dpng',[figPath,'commandDipole',sprintf('_c%1.0f_x%1.0f_g%1.0f',controlLaw,ICi,gain*10),'.png']);

    end % initial conditions loop
end % control law loop


end % run_controlSim




warning('Need to try various gains')

% -----------------------------------------------------------------
%%% 3.5 - Orbit Inclinations vs Control Performance
% -----------------------------------------------------------------
if run_inclinationSim == 1
%%% Min Torque switch
run_minTorque = 0;

%%% Time vector
t0 = 0; % sec
dt = .1; % sec

%%% Titles for cases
controlTitles = {'Modulating B-dot Control Law','Bang-Bang B-dot Control Law'};

%%% Inclinations
incs = [15, 105].*pi/180; % rad
incLabels = {'15','105'};

%%% Setting integration times so that 3 deg/s is met
tfs = [4.5, 4];

if run_minTorque == 0
    incRange = 2;
elseif run_minTorque == 1
    incRange = 1;
end

for incs_i = 1:incRange % 1:2
    in0 = incs(incs_i);
    
    %%% Setting times vector
    if run_minTorque == 0
        tf = Tp*tfs(incs_i); % sec
    elseif run_minTorque == 1
        tf = Tp*3; % sec
    end
    times = t0:dt:tf;
    
    %%% Preallocating
    bs_H = zeros(3,length(times));
    bs_B = zeros(3,length(times));
    states = zeros(6,length(times));
    controls = zeros(3,length(times));
    commandDipoles = zeros(3,length(times));

    for controlLaw = 1:2 % 1:2
        %%% Storing ICs
        states(:,1) = [0.3; 0.2; 0.4; [15;8;12].*pi/180];

        for kk = 1:length(times)-1
            ti = times(kk); % sec

            % ---------------------
            % Magnetic Field Info
            % ---------------------
            %%% Calculate Bm
            Bm = Bm0 + wE*ti; % rad

            %%% Calculate b_H
            zeta = acos(cos(in0)*cos(gm) + sin(in0)*sin(gm)*cos(raan0-Bm));
            eta = asin(sin(gm)*sin(raan0-Bm)/sin(zeta)); 
            b_H = (M/r^3)*[cos(th_dot*ti-eta)*sin(zeta); cos(zeta); -2*sin(th_dot*ti-eta)*sin(zeta)];

            bs_H(:,kk) = b_H;

            %%% Calculate [BH]
            [ BN ] = mrp2DCM( states(1:3,kk) );
            HN = [cos(raan0)*cos(th_dot*ti) - cos(in0)*sin(raan0)*sin(th_dot*ti), cos(th_dot*ti)*sin(raan0) + cos(in0)*cos(raan0)*sin(th_dot*ti), sin(in0)*sin(th_dot*ti);...
                  -cos(raan0)*sin(th_dot*ti) - cos(in0)*cos(th_dot*ti)*sin(raan0), cos(in0)*cos(raan0)*cos(th_dot*ti) - sin(raan0)*sin(th_dot*ti), cos(th_dot*ti)*sin(in0);...
                  sin(in0)*sin(raan0), -cos(raan0)*sin(in0), cos(in0)];
            BH = BN*HN';

            %%% Calculate b_B
            b_B = BH*b_H;
            bs_B(:,kk) = b_B;

            % ---------------------
            % Control
            % ---------------------
            %%% Set gain
            gain = 1;

            %%% Calculate control
            if run_minTorque == 0
                mMax = 3; % Am^2
            elseif run_minTorque == 1
                if controlLaw == 1
                    mMax = 4.24; % Am^2 ... modulated
                elseif controlLaw == 2
                    mMax = 4.05; % Am^2 .... bang bang
                end                
            end

            if controlLaw == 1
                [u,m] = control_modulatingBdot(th_dot,zeta,I1,bs_B(:,kk),states(4:6,kk),mMax,gain);
            elseif controlLaw == 2
                if kk > 1
                    [u,m] = control_bangBang(mMax,bs_B(:,kk-1),bs_B(:,kk),dt);
                else
                    u = zeros(3,1);
                    m = zeros(3,1);
                end
            end

            controls(:,kk) = u;
            commandDipoles(:,kk) = m;
            % ---------------------
            % Integrating states
            % ---------------------
            %%% Checking MRPs for shadow set switch
            s = norm(states(1:3,kk));
            if s > 1
                states(1:3,kk) = -states(1:3,kk)./(s^2);
            end

            [X_p] = rk4_MRPw(states(:,kk),dt,I,controls(:,kk));
            states(:,kk+1) = X_p;

        end % times loop
        
        
        %%% Finding index where error threshold is crossed
        controlError = colnorm(states(4:6,:));
        if run_minTorque == 0
            errorThresh = 3; % deg/s
        elseif run_minTorque == 1
            errorThresh = 1; % deg/s
        end
        errorIndex = find(controlError<=(errorThresh*pi/180));
        if isempty(errorIndex) == 0
            errorIndex = errorIndex(1);
            fprintf('c%1.0f_inc%s_g%1.1f - Threshold crossing: %2.2f min\n',controlLaw,incLabels{incs_i},gain*10,times(errorIndex)/60)
        else
            fprintf('c%1.0f_inc%s_g%1.1f - No Threshold crossing\n',controlLaw,incLabels{incs_i},gain*10)
        end

        %%% Plotting control error
        figure; hold all
        plot(times./Tp,controlError,'linewidth',2,'color','k')
        plot([times(1),times(end)]./Tp,[errorThresh*pi/180, errorThresh*pi/180],'--r','linewidth',2)
        if isempty(errorIndex) == 0
            plot([times(errorIndex),times(errorIndex)]./Tp,[0, max(controlError)],'--m','linewidth',2)
        end
        PlotBoi2('Time, period','Control Error, rad/s',18)
        title(sprintf('%s  ...  Inclination: %s°  ...  Gain: %1.1f',controlTitles{controlLaw},incLabels{incs_i},gain));
        legend('|\omega| - control error',sprintf('%1.0f°/s threshold',errorThresh),'threshold crossing');
        if run_minTorque == 0
            print(gcf,'-dpng',[figPath,'controlError',sprintf('_c%1.0f_i%s_g%2.0f',controlLaw,incLabels{incs_i},gain*10),'.png']);
        elseif run_minTorque == 1
            print(gcf,'-dpng',[figPath,'controlErrorMin',sprintf('_c%1.0f_i%s_g%2.0f',controlLaw,incLabels{incs_i},gain*10),'.png']);
        end
        
    end % control law loop
end % inclination loop

end % run_inclinationSim

% -----------------------------------------------------------------
%%% 3.6 - Monte Carlo Analysis
% -----------------------------------------------------------------
if run_monteCarloAnalysis == 1
tic
%%% Time vector
t0 = 0; % sec
dt = .1; % sec
tf = Tp*5; % sec
times = t0:dt:tf;

%%% Preallocating
bs_H = zeros(3,length(times));
bs_B = zeros(3,length(times));
states = zeros(6,length(times));
controls = zeros(3,length(times));
commandDipoles = zeros(3,length(times));

%%% Titles for cases
controlTitles = {'Modulating B-dot Control Law','Bang-Bang B-dot Control Law'};

%%% Max control
mMax = 4; % Am^2
    
for controlLaw = 1:2 % 1:2
    figure
    
    for mc = 1:25 % 1:25
        %%% Setting ICs
        w0_rand = rand(3,1)*6 + 10;
        signs_rand = 2*round(rand(3,1))-1;
        w0_rand = w0_rand .* signs_rand;
        states(:,1) = [0.3; 0.2; 0.4; w0_rand.*pi/180];
        
        for kk = 1:length(times)-1
            ti = times(kk); % sec

            % ---------------------
            % Magnetic Field Info
            % ---------------------
            %%% Calculate Bm
            Bm = Bm0 + wE*ti; % rad

            %%% Calculate b_H
            zeta = acos(cos(in0)*cos(gm) + sin(in0)*sin(gm)*cos(raan0-Bm));
            eta = asin(sin(gm)*sin(raan0-Bm)/sin(zeta)); 
            b_H = (M/r^3)*[cos(th_dot*ti-eta)*sin(zeta); cos(zeta); -2*sin(th_dot*ti-eta)*sin(zeta)];

            bs_H(:,kk) = b_H;

            %%% Calculate [BH]
            [ BN ] = mrp2DCM( states(1:3,kk) );
            HN = [cos(raan0)*cos(th_dot*ti) - cos(in0)*sin(raan0)*sin(th_dot*ti), cos(th_dot*ti)*sin(raan0) + cos(in0)*cos(raan0)*sin(th_dot*ti), sin(in0)*sin(th_dot*ti);...
                  -cos(raan0)*sin(th_dot*ti) - cos(in0)*cos(th_dot*ti)*sin(raan0), cos(in0)*cos(raan0)*cos(th_dot*ti) - sin(raan0)*sin(th_dot*ti), cos(th_dot*ti)*sin(in0);...
                  sin(in0)*sin(raan0), -cos(raan0)*sin(in0), cos(in0)];
            BH = BN*HN';

            %%% Calculate b_B
            b_B = BH*b_H;
            bs_B(:,kk) = b_B;

            % ---------------------
            % Control
            % ---------------------
            %%% Set gain
            gain = 1;

            %%% Calculate control
            if controlLaw == 1
                [u,m] = control_modulatingBdot(th_dot,zeta,I1,bs_B(:,kk),states(4:6,kk),mMax,gain);
            elseif controlLaw == 2
                if kk > 1
                    [u,m] = control_bangBang(mMax,bs_B(:,kk-1),bs_B(:,kk),dt);
                else
                    u = zeros(3,1);
                    m = zeros(3,1);
                end
            end

            controls(:,kk) = u;
            commandDipoles(:,kk) = m;
            % ---------------------
            % Integrating states
            % ---------------------
            %%% Checking MRPs for shadow set switch
            s = norm(states(1:3,kk));
            if s > 1
                states(1:3,kk) = -states(1:3,kk)./(s^2);
            end

            [X_p] = rk4_MRPw(states(:,kk),dt,I,controls(:,kk));
            states(:,kk+1) = X_p;

        end % times loop
        
        %%% Plotting control error
        controlError = colnorm(states(4:6,:));
        plot(times./Tp,controlError,'linewidth',1,'color','k'); hold all
        
        controlLaw
        mc
    end % monte carlo loop
    
    %%% Touching up plot
    PlotBoi2('Time, period','Control Error, rad/s',18)
    title(sprintf('Monte Carlo - %s  ...  Gain: %1.1f',controlTitles{controlLaw},gain));
    legend('|\omega| - control error');
    set(gca,'yscale','log')
    print(gcf,'-dpng',[figPath,'monteCarlo',sprintf('_c%1.0f_g%2.0f',controlLaw,gain*10),'.png']);
    
        
end % control Law loop

toc
end % run_monteCarloAnalysis















fprintf('\n\nInclude: EOM stuff like u/I; Tp in minutes, effects of gain, control probably unrealistic since doing it at every .1 sec \n')
fprintf('-Question on P5; First part of Inclination problem, thresholds are for 1 deg/s!!!;  \n')
fprintf('-What is control error?\n')
%% =======================================================================
%%% Functions
% ========================================================================


function [u,m] = control_modulatingBdot(th_dot,zeta,Imin,b,w,mMax,gain)
%%% Inputs
%    1) th_dot - [1x1] mean motion
%    2) zeta - [1x1]
%    3) Imin - [1x1] minimum principal inertia
%    4) b - [3x1] local magnetic field in the body frame
%    5) w - [3x1] angular velocity
%    6) mMax - [1x1] max command dipole
%    7) gain - gain for control

%%% Output
%    1) u - control vector
% ============================================================
%%% Calculate gain
kw = gain*2*th_dot*(1+sin(zeta))*Imin;

%%% Calculate command dipole
bhat = b./norm(b);
m = cross(-kw*bhat./norm(b),(eye(3)-bhat*bhat')*w);

%%% Checking that m is feasible
for kk = 1:3
    if abs(m(kk)) > mMax
        m(kk) = sign(m(kk))*mMax;
    end
end

%%% Calculate control
u = cross(m,b);

end

function [u,m] = control_bangBang(mMax,bOld,bNew,dt)
%%% Inputs
%    1) mMax - [1x1] max command dipole
%    2) bNew - [3x1] local magnetic field in the body frame at time t
%    3) bOld - [3x1] local magnetic field in the body frame at time t-1
%    4) dt   - time between bNew and bOld
%%% Output
%    1) u - control vector
% ============================================================
%%% Time derivative of magnetic field in body frame
bPrime = (bNew-bOld)/dt;

%%% Command dipole
m = -mMax *sign(bPrime);

%%% Control
u = cross(m,bNew);

end










