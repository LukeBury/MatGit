clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
get_colors();
% ========================================================================
%%% Assumptions
% ========================================================================
% - Stationary secondary

% ========================================================================
%%% Run-Switches
% ========================================================================
run_altitudeShiftAnalysis  = 1;
run_flybyEffectOnPeriapsis = 1;
% ========================================================================
%%% Importing Data .. Choosing Settings
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Selecting Bodies
primary = bodies.mars;
secondary = bodies.phobos;

%%% Constant altitude of flybys
alt_flyby = 0.5; % km

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Analysis
% ========================================================================
% -------------------------------------------------
% Target Orbit and Intercept Parameters
% -------------------------------------------------
target.ra = 33793 + primary.R; % 37183 km
target.rp = 250 + primary.R; % 3640 km
target.a = (target.ra + target.rp)/2;
target.e = (target.ra - target.rp)/(target.ra + target.rp);
target.energy = -primary.u/(2*target.a); % km^2/s^2
target.va = visviva_v( target.ra, target.a, primary.u); % km/s
target.vp = visviva_v( target.rp, target.a, primary.u); % km/s

%%% r & v of spacecraft at target orbit flyby
r_tgt_int = secondary.a;
v_tgt_int = sqrt(primary.u*(2/r_tgt_int - 1/target.a));

%%% velocity of secondary in circular orbit
secondary.v = sqrt(primary.u/secondary.a);

%%% Target Intercept Anomalies
E_tgt_int = acos(1/target.e - r_tgt_int/(target.a*target.e)); % rad
ta_tgt_int = 2*atan2(sqrt(1+target.e)*tan(E_tgt_int/2),sqrt(1-target.e)); % rad

%%% Stationary secondary positon (for target flyby)
rSecondary_tgt = R3([secondary.a, 0,0],ta_tgt_int);

%%% Orbital period of secondary
secondary.P = 2*pi*sqrt(secondary.a^3/secondary.u); % sec

%% =======================================================================
% Analysis on how dV at periapsis can affect flyby altitude
% ========================================================================
if run_altitudeShiftAnalysis == 1
    
% -------------------------------------------------
% Setting up and Running Simulation
% -------------------------------------------------
%%% Initial States
% Nominal (target)
[r0, v0] = OE2ECI(target.a, target.e, 0, 0, 0, 0, primary.u);
X0 = [r0;v0];

if secondary.name == 'Phobos'
    % Plus 
    X0_p = [r0;v0 + [0;.000095;0]];

    % Minus
    X0_m = [r0;v0 - [0;.000095;0]];
    
    % Final time
    tf = 1*3600; % sec
elseif secondary.name == 'Deimos'
    % Plus 
    X0_p = [r0;v0 + [0;.000013;0]];
    % Minus
    X0_m = [r0;v0 - [0;.000013;0]];
    
    tf = 4*3600; % sec
end
%%% Setting time vector and normalizing 
ti = 0; % sec
% dt = .001; % sec
dt = .01; % sec

time = ti:dt:tf;

%%% Choosing ode45 tolerance
tol = 1e-9;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating the States
[timef, X_tgt] = ode45(@Int_2BI, time, X0, options, primary.u);
[timef, X_p] = ode45(@Int_2BI, time, X0_p, options, primary.u);
[timef, X_m] = ode45(@Int_2BI, time, X0_m, options, primary.u);

% -------------------------------------------------
% Plotting Mars-Centric Simulation Results
% -------------------------------------------------
hold all
plot(X_tgt(:,1),X_tgt(:,2),'k')
plot(X_p(:,1),X_p(:,2),'g')
plot(X_m(:,1),X_m(:,2),'r')
plotBody2( secondary.a, [0,0], [1,1,1,0],colors.std.red,1) % secondary orbit
plotBody2( primary.R, [0,0], primary.color,primary.color,1) % mars
plotBody2( secondary.R, rSecondary_tgt(1:2), secondary.color,secondary.color,1) % secondary
axis equal

% -------------------------------------------------
% Performing Altitude analysis
% -------------------------------------------------
alts_tgt = zeros(size(X_tgt,1),1);
alts_p = zeros(size(X_tgt,1),1);
alts_m = zeros(size(X_tgt,1),1);
for kk = 1:size(X_tgt,1)
    alts_tgt(kk) = sqrt((X_tgt(kk,1:3)-rSecondary_tgt)*(X_tgt(kk,1:3)-rSecondary_tgt)'); % km
    alts_p(kk) = sqrt((X_p(kk,1:3)-rSecondary_tgt)*(X_p(kk,1:3)-rSecondary_tgt)'); % km
    alts_m(kk) = sqrt((X_m(kk,1:3)-rSecondary_tgt)*(X_m(kk,1:3)-rSecondary_tgt)'); % km
end

minAlt_tgt = min(alts_tgt) % minimum altitude of flyby
minAlt_p   = min(alts_p)   % minimum altitude of flyby
minAlt_m   = min(alts_m)   % minimum altitude of flyby

desired = .5  % aiming for altitude of 500 m for (plus) and (minus) trajs


ti_tgt = find(alts_tgt==min(alts_tgt)); % time index of C/A
ti_p = find(alts_p==min(alts_p));       % time index of C/A
ti_m = find(alts_m==min(alts_m));       % time index of C/A

t_tgt = timef(ti_tgt)  % time of C/A
t_p = timef(ti_p)      % time of C/A
t_m = timef(ti_m)      % time of C/A

t_tgt-t_p % time difference of C/A
t_tgt-t_m % time difference of C/A

end % run_altitudeShiftAnalysis

%% =======================================================================
% Analysis on how flybys can affect periapsis velocity
% ========================================================================
if run_flybyEffectOnPeriapsis == 1
% -------------------------------------------------
% Preparing for analysis
% -------------------------------------------------
%%% Assigning independent variables
rApos = linspace(target.ra,secondary.a+.001,1000)'; % km, from target orbit to last possible flyby

%%% interception distance of spacecraft and secondary
rInt = secondary.a; % km

%%% Constant periapse radius
rp = target.rp; % km

%%% Preallocating
vInfs     = zeros(size(rApos)); % km/s
alphas    = zeros(size(rApos)); % rad; angle of traj. wrt secondary
betas     = zeros(size(rApos)); % rad
gammas    = zeros(size(rApos)); % rad
deltas    = zeros(size(rApos)); % rad; turning angle
vScs_m    = zeros(size(rApos)); % km/s; s/c velocity before flyby
vScs_p    = zeros(size(rApos)); % km/s; s/c velocity after flyby
dVs       = zeros(size(rApos)); % km/s; mars-centric dV from flyby
vPers_m   = zeros(size(rApos)); % km/s; periapsis vel. of trajecotry before flyby
vPers_p   = zeros(size(rApos)); % km/s; periapsis vel. of trajecotry after flyby
dVpers    = zeros(size(rApos)); % km/s; mars-centric periapsis dV from flyby
periods_m = zeros(size(rApos)); % s; period of full orbit before flyby
% -------------------------------------------------
% Looping through Apoapses, find dV, relating that to periapsis speed
% -------------------------------------------------
betaSwitch = 0;
betaFulcrumIndex = 0;
for ii = 1:length(rApos)
    %%% Assigning current apoapsis
    ra_m = rApos(ii); % km
    
    %%% Calculating incoming alpha angle
    a_m = (ra_m+rp)./2; % km
    e_m = (ra_m - rp)/(ra_m + rp);
    E_m = acos(1/e_m - rInt/(a_m*e_m)); % rad
    ta_m = 2*atan2(sqrt(1+e_m)*tan(E_m/2),sqrt(1-e_m)); % rad
    [ fpa_m ] = oe2fpa( e_m, ta_m); % rad
    alphas(ii) = fpa_m; % rad
    
    %%% Calculating period of full orbit before flyby
    periods_m(ii) = 2*pi*sqrt(a_m^3/secondary.u); % sec
    
    %%% Calculating incoming vInf
    vScs_m(ii) = visviva_v( rInt, a_m, primary.u); % km/s
    vInfs(ii) = sqrt(secondary.v^2 + vScs_m(ii)^2 - 2*secondary.v*vScs_m(ii)*cos(alphas(ii))); % km/s 
    
    %%% Calculating Beta and Gamma angles
    betas(ii) = acos((vInfs(ii)^2 + vScs_m(ii)^2 - secondary.v^2)/(2*vInfs(ii)*vScs_m(ii)));
    gammas(ii) = pi - betas(ii) - alphas(ii); % rad
    
    %%% Calculating turning angle
    deltas(ii) = calcTurningAngle( alt_flyby, vInfs(ii), secondary.u ); % rad
    
    %%% Calculating final mars-centric spacecraft velocity
    vScs_p(ii) = sqrt(vInfs(ii)^2 + secondary.v^2 - 2*vInfs(ii)*secondary.v*cos(gammas(ii)-deltas(ii))); % km/s
    
    %%% Calculating change in mars-centric velocity
    dVs(ii) = vScs_p(ii) - vScs_m(ii);
    
    %%% Calculating change in periapsis velocity
    a_p = visviva_a( rInt, vScs_p(ii), primary.u); % km
    vPers_m(ii) = visviva_v( rp, a_m, primary.u); % km/s
    vPers_p(ii) = visviva_v( rp, a_p, primary.u); % km/s
    dVpers(ii) = vPers_p(ii) - vPers_m(ii); % km/s
end

% -------------------------------------------------
% Finding rApos with periods near resonance of secondary
% -------------------------------------------------
if secondary.name == 'Phobos'
    resonantPeriods = [secondary.P, secondary.P*2, secondary.P*3];
    resonantPeriodsIndices = [find(abs(periods_m-secondary.P)==min(abs(periods_m-secondary.P))), find(abs(periods_m-secondary.P*2)==min(abs(periods_m-secondary.P*2))),find(abs(periods_m-secondary.P*3)==min(abs(periods_m-secondary.P*3)))];
    resonantPeriodsLegend = {'Sec:Traj, 1:1','Sec:Traj, 2:1','Sec:Traj, 3:1'};
elseif secondary.name == 'Deimos'
    resonantPeriods = [secondary.P/2];
    resonantPeriodsIndices = [find(abs(periods_m-secondary.P/2)==min(abs(periods_m-secondary.P/2)))];
    resonantPeriodsLegend = {'Sec:Traj, 1:2'};
end

% -------------------------------------------------
% Plotting
% -------------------------------------------------
figure; hold all % rApos vs alpha
[p] = plotResPeriods(secondary.name, resonantPeriodsIndices, rApos, alphas.*180/pi);
plot(rApos,alphas.*180/pi,'linewidth',2)
title(sprintf('%s, flyby altitude: %1.1f km',secondary.name,alt_flyby))
legend(p,resonantPeriodsLegend)
PlotBoi2('rApo, km','alpha, °',14)
print(gcf,'-dpdf','rAposVsAlpha.pdf');

figure; hold all % rApos vs beta
[p] = plotResPeriods(secondary.name, resonantPeriodsIndices, rApos, betas.*180/pi);
plot(rApos,betas.*180/pi,'linewidth',2)
title(sprintf('%s, flyby altitude: %1.1f km',secondary.name,alt_flyby))
legend(p,resonantPeriodsLegend)
PlotBoi2('rApo, km','beta, °',14)
% print(gcf,'-dpdf','rAposVsBeta.pdf');

figure; hold all % rApos vs vInf
[p] = plotResPeriods(secondary.name, resonantPeriodsIndices, rApos, vInfs);
plot(rApos,vInfs,'linewidth',2)
title(sprintf('%s, flyby altitude: %1.1f km',secondary.name,alt_flyby))
legend(p,resonantPeriodsLegend)
PlotBoi2('rApo, km','vInf, km/s',14)
print(gcf,'-dpdf','rAposVsVinf.pdf');

figure; hold all % rApos vs vSc_minus
[p] = plotResPeriods(secondary.name, resonantPeriodsIndices, rApos, vScs_m);
plot(rApos,vScs_m,'linewidth',2)
title(sprintf('%s, flyby altitude: %1.1f km',secondary.name,alt_flyby))
legend(p,resonantPeriodsLegend)
PlotBoi2('rApo, km','vSc_{minus}, km/s',14)
% print(gcf,'-dpdf','rAposVsVsc_minus.pdf');

% figure; hold all % rApos vs vSc_plus
% [p] = plotResPeriods(secondary.name, resonantPeriodsIndices, rApos, vScs_p);
% plot(rApos,vScs_p,'linewidth',2)
% title(sprintf('%s, flyby altitude: %1.1f km',secondary.name,alt_flyby))
% legend(p,resonantPeriodsLegend)
% PlotBoi2('rApo, km','vSc_{plus}, km/s',14)
% print(gcf,'-dpdf','rAposVsVsc_plus.pdf');

figure; hold all % rApos vs turning angle
[p] = plotResPeriods(secondary.name, resonantPeriodsIndices, rApos, deltas.*180/pi);
plot(rApos,deltas.*180/pi,'linewidth',2)
title(sprintf('%s, flyby altitude: %1.1f km',secondary.name,alt_flyby))
legend(p,resonantPeriodsLegend)
PlotBoi2('rApo, km','turning angle, °',14)
print(gcf,'-dpdf','rAposVsTurningAngle.pdf');

figure; hold all % rApos vs dV
[p] = plotResPeriods(secondary.name, resonantPeriodsIndices, rApos, dVs.*1000);
plot(rApos,dVs.*1000,'linewidth',2)
title(sprintf('%s, flyby altitude: %1.1f km',secondary.name,alt_flyby))
legend(p,resonantPeriodsLegend)
PlotBoi2('rApo, km','mars-centric dV, m/s',14)
print(gcf,'-dpdf','rAposVsDV.pdf');

figure; hold all % rApos vs dVper
[p] = plotResPeriods(secondary.name, resonantPeriodsIndices, rApos, dVpers.*1000);
plot(rApos,dVpers.*1000,'linewidth',2)
title(sprintf('%s, flyby altitude: %1.1f km',secondary.name,alt_flyby))
legend(p,resonantPeriodsLegend)
PlotBoi2('rApo, km','periapsis dV, m/s',14)
print(gcf,'-dpdf','rAposVsDVper.pdf');

figure; hold all % rApos vs periods
[p] = plotResPeriods(secondary.name, resonantPeriodsIndices, rApos, periods_m./3600);
plot(rApos,periods_m./3600,'linewidth',2)
title(sprintf('%s, flyby altitude: %1.1f km',secondary.name,alt_flyby))
legend(p,resonantPeriodsLegend)
PlotBoi2('rApo, km','periods, hr',14)
% print(gcf,'-dpdf','rAposVsPeriods.pdf');

%%% Combining PDFs and deleting individuals
outputPDF_name = sprintf('MarsAnalysis_%s_%1.1fkm.pdf',secondary.name,alt_flyby);
if exist(outputPDF_name, 'file') == 2
    delete(outputPDF_name);
end
append_pdfs(outputPDF_name,'rAposVsDVper.pdf','rAposVsDV.pdf','rAposVsTurningAngle.pdf','rAposVsVinf.pdf','rAposVsAlpha.pdf')
delete('rAposVsDVper.pdf','rAposVsDV.pdf','rAposVsTurningAngle.pdf','rAposVsVinf.pdf','rAposVsAlpha.pdf')

distFig('Screen','Main') % 'Main' or 'Secondary'

end % run_flybyEffectOnPeriapsis



%% =======================================================================
% Function for plotting dashed lines at rApo values that signify the orbit
% has a period in some resonance with the secondary
% ========================================================================
function [p] = plotResPeriods(secondaryName, resonantPeriodsIndices, rApos, yAxisData)
if secondaryName == 'Deimos'
    p1 = plot([rApos(resonantPeriodsIndices),rApos(resonantPeriodsIndices)], [min(yAxisData)-0.1*abs(max(yAxisData)-min(yAxisData)),max(yAxisData)+0.1*abs(max(yAxisData)-min(yAxisData))],'--xr','linewidth',1,'markersize',8);
    p = [p1];
elseif secondaryName == 'Phobos'
    p1 = plot([rApos(resonantPeriodsIndices(1)),rApos(resonantPeriodsIndices(1))], [min(yAxisData)-0.1*abs(max(yAxisData)-min(yAxisData)),max(yAxisData)+0.1*abs(max(yAxisData)-min(yAxisData))],'--xr','linewidth',1,'markersize',8);
    p2 = plot([rApos(resonantPeriodsIndices(2)),rApos(resonantPeriodsIndices(2))], [min(yAxisData)-0.1*abs(max(yAxisData)-min(yAxisData)),max(yAxisData)+0.1*abs(max(yAxisData)-min(yAxisData))],'--or','linewidth',1,'markersize',8);
    p3 = plot([rApos(resonantPeriodsIndices(3)),rApos(resonantPeriodsIndices(3))], [min(yAxisData)-0.1*abs(max(yAxisData)-min(yAxisData)),max(yAxisData)+0.1*abs(max(yAxisData)-min(yAxisData))],'--^r','linewidth',1,'markersize',8);
    p = [p1,p2,p3];
end
end














