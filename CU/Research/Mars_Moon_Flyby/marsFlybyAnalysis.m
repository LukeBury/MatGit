clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
get_colors();
% ========================================================================
%%% Information
% ========================================================================
% -------------------------------------------------
% Questions?
% -------------------------------------------------
% -How much can we pump down energy per flyby?
% -Flyby altitudes? (down to ~1 km)
%
% -------------------------------------------------
% Target
% -------------------------------------------------
% -Target orbit after aerocapture: 33793 km x 250 km centered on mars
%
% -------------------------------------------------
% Assumptions
% -------------------------------------------------
% -For secondary intersect, ignoring eccentricity
%
% ========================================================================
%%% Run-Switches
% ========================================================================
run_tangentAnalysis      = 0;
run_nonTangentAnalysis   = 1;
run_loopingFlybyAnalysis = 1;
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Selecting Bodies
primary = bodies.mars;
secondary = bodies.phobos;

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Analysis
% ========================================================================
% -------------------------------------------------
% Target Orbit Parameters
% -------------------------------------------------
target.ra = 33793 + primary.R; % 37183 km
target.rp = 250 + primary.R; % 3640 km
target.a = (target.ra + target.rp)/2;
target.e = (target.ra - target.rp)/(target.ra + target.rp);
target.energy = -primary.u/(2*target.a); % km^2/s^2

% -------------------------------------------------
% Finding orbit parameters where target orbit intersects moon orbit
% -------------------------------------------------
%%% r & v of spacecraft
r_int = secondary.a;
v_int = sqrt(primary.u*(2/r_int - 1/target.a));

%%% v of phobos
v_secondary = sqrt(primary.u*(2/r_int - 1/secondary.a));

% -------------------------------------------------
% Setting independent variables (v_sc and r_alt)
% -------------------------------------------------
%%% Setting 'n' for nxn square mesh
% n for tangent analysis
n1 = 30;

% n for non-tangent analysis
n2 = 50; 

%%% Mars-centric inertial velocity of spacecraft at flyby
v_infs = linspace(abs(v_secondary-v_int), abs(v_secondary-v_int)+.5, n1); % km/s

%%% Altitude of spacecraft above secondary at flyby
altitudes = linspace(.5, 10, n1); % km

%%% Angle between incoming trajectory and circular secondary orbit (180 degrees if tangent)
% alphas = fliplr(linspace(180,130,n2)).*pi/180; % rad
alphas = fliplr(linspace(0,50,n2)).*pi/180; % rad
% -------------------------------------------------
% Running Tangent Analyses: Turning angle, velocity change, energy change
% -------------------------------------------------
if run_tangentAnalysis == 1
%%% Creating mesh grid for values
[X, Y] = meshgrid(v_infs,altitudes);

%%% Preallocating results
turningAngles_tan = zeros(size(X));
velocityChanges_tan = zeros(size(X));
energyChanges_tan = zeros(size(X));

for v_i = 1:size(X,1)
    for alt_i = 1:size(X,2)
        %%% Turning angle analysis
        turningAngles_tan(v_i,alt_i) = calcTurningAngle( altitudes(alt_i), v_infs(v_i), secondary.u ); % rad
        
        %%% Inertial, Mars-centric velocity change analysis
        v1 = v_infs(v_i) + v_secondary;
        v2 = sqrt(v_infs(v_i)^2 + v_secondary^2 - 2*v_infs(v_i)*v_secondary*cos(pi-turningAngles_tan(v_i,alt_i)));
        velocityChanges_tan(v_i,alt_i) = (v2-v1); % km/s
        
        %%% Energy change analysis
        energy1 = v1^2/2 - primary.u/(altitudes(alt_i) + secondary.a);
        energy2 = v2^2/2 - primary.u/(altitudes(alt_i) + secondary.a);
        energyChanges_tan(v_i,alt_i) = v2^2/2 - v1^2/2; % km^2/s^2
        
    end
end

% Converting to desired unites
turningAngles_tan = turningAngles_tan.*180/pi;   % rad  -> degrees
velocityChanges_tan = velocityChanges_tan.*1000; % km/s -> m/s

% Plotting surface of turning angles
figure
surf(v_infs,altitudes,turningAngles_tan)
title([secondary.name,' - Turning Angles, deg'])
xlabel('Spacecraft v-infinity, km/s')
ylabel('Flyby altitude, km')
view(0,90)
c = colorbar; ylabel(c, 'degrees','fontsize',14)
caxis([0,2])
print(gcf,'-dpdf','turningAngleAnalysis.pdf');

% Plotting surface of velocity changes
figure
surf(v_infs,altitudes,velocityChanges_tan)
title([secondary.name,' - Inertial, Mars-Centric \DeltaV from flyby, m/s'])
xlabel('Spacecraft v-infinity, km/s')
ylabel('Flyby altitude, km')
view(0,90)
c = colorbar; ylabel(c, 'm/s','fontsize',14)
caxis([-0.06, 0])
print(gcf,'-dpdf','velocityChangeAnalysis.pdf');

% Plotting surface of energy changes
figure
surf(v_infs,altitudes,energyChanges_tan)
title([secondary.name,' - \DeltaEnergy from flyby, km^2/s^2'])
xlabel('Spacecraft v-infinity, km/s')
ylabel('Flyby altitude, km')
view(0,90)
c = colorbar; ylabel(c, 'km^2/s^2','fontsize',14)
caxis([-9e-5, 0])
print(gcf,'-dpdf','energyChangeAnalysis.pdf');

% Combining PDFs and deleting individuals
outputPDF_name = sprintf('output_%s.pdf',secondary.name);
delete(outputPDF_name);
append_pdfs(outputPDF_name,'turningAngleAnalysis.pdf','velocityChangeAnalysis.pdf','energyChangeAnalysis.pdf')
delete('turningAngleAnalysis.pdf','velocityChangeAnalysis.pdf','energyChangeAnalysis.pdf')

end % run_tangentAnalysis
% -------------------------------------------------
% Running Non-Tangent Analyses: Alpha angle
% -------------------------------------------------
if run_nonTangentAnalysis == 1
%%% Choosing constant v_inf & altitude for Alpha analysis
if secondary.name == 'Deimos'
    v_inf_const = .11; % km/s ... Deimos
elseif secondary.name == 'Phobos'
    v_inf_const = .52; % km/s ... Phobos
end
alt_const = .5; % km/s
turningAngle_const = calcTurningAngle( alt_const, v_inf_const, secondary.u ); % rad

%%% Preallocating results
velocityChanges_nontan = zeros(size(alphas));
energyChanges_nontan = zeros(size(alphas));

kk = 0;
for alpha_i = alphas
    kk = kk+1;
    %%% Inertial, Mars-centric velocity change analysis
    v1 = sqrt(v_inf_const^2 + v_secondary^2 - 2*v_inf_const*v_secondary*cos(pi-alpha_i));
    v2 = sqrt(v_inf_const^2 + v_secondary^2 - 2*v_inf_const*v_secondary*cos(pi-alpha_i-turningAngle_const));
    velocityChanges_nontan(kk) = v2-v1; % km/s

    %%% Energy change analysis
    energy1 = v1^2/2 - primary.u/(alt_const + secondary.a);
    energy2 = v2^2/2 - primary.u/(alt_const + secondary.a);
    energyChanges_nontan(kk) = v2^2/2 - v1^2/2; % km^2/s^2
end

figure
subplot(2,1,1)
plot(alphas.*180/pi,velocityChanges_nontan.*1000,'linewidth',2,'color',colors.std.purp)
title('Change in Mars-centric velocity')
PlotBoi2('','\DeltaV, m/s',14)

subplot(2,1,2)
plot(alphas.*180/pi,energyChanges_nontan,'linewidth',2,'color',colors.std.purp)
title('Change in Mars-centric energy')
PlotBoi2('\alpha, deg','\DeltaE, km^2/s^2',14)





end % run_nonTangentAnalysis


if run_loopingFlybyAnalysis == 1
% -------------------------------------------------
% Running Non-Tangent Analyses: Alpha angle
% -------------------------------------------------
%%% Altitude of flybys
alt_flyby = .5; % km/s

%%% interception distance of spacecraft and secondary
r_int = secondary.a;

%%% Values shrinking from target orbit to tangent intercept
ras = linspace(secondary.a,target.ra,100);
rp = target.rp;

%%% preallocating
alphas = zeros(size(ras));
dVs = zeros(size(ras));
dVps = zeros(size(ras)); % periapse velocity change
figure; hold all
%%% Looping through orbits
for kk = 1:length(ras)
%%% Incoming alpha ange (FPA)
ra = ras(kk); % km
a = (ra+rp)./2; % km
e = (ra - rp)/(ra + rp);
E = acos(1/e - r_int/(a*e)); % rad
v = 2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e)); % rad
[ fpa ] = oe2fpa( e, v); % rad
alphas(kk) = fpa; % rad

%%% Spacecraft velocity before flyby
[v_1] = visviva_v( r_int, a, primary.u); % km/s

%%% Calculate incoming vInf
gamma = v_secondary*(pi-alphas(kk))/v_1; % rad
% vInf = sqrt(v_secondary^2 + v_1^2 -2*v_secondary*v_1*cos(pi-alphas(kk)-gamma));
vInf = sqrt(v_secondary^2 + v_1^2 -2*v_secondary*v_1*cos(pi - (pi - alphas(kk))-gamma));
plot(ra,vInf,'ro')
%%% Calculate turning angle
turningAngle = calcTurningAngle( alt_flyby, vInf, secondary.u ); % rad
%%% Inertial, Mars-centric velocity change analysis
v_2 = sqrt(vInf^2 + v_secondary^2 - 2*vInf*v_secondary*cos(pi-alphas(kk)-turningAngle));
dVs(kk) = v_2-v_1; % km/s
fprintf('-------- %3.0f\n%2.3f - vInf (km/s)\n%2.3f - alpha (deg)\n',ras(kk),vInf,alphas(kk)*180/pi)
test_vInf = alphas(kk)*v_1/(pi-alphas(kk)) - v_secondary;
fprintf('%2.3f - new vInf\n\n',test_vInf)
end
figure
plot(ras,alphas.*180./pi,'linewidth',2,'color',colors.std.purp)
PlotBoi2('R_a, km','\alpha, deg',14)

figure
plot(ras,dVs.*1000,'linewidth',2,'color',colors.std.purp)
PlotBoi2('R_a, km','dV, m/s',14)

warning('Something wrong with vInf calculation')
% target.ra = 33793 + primary.R; % 37183 km
% target.rp = 250 + primary.R; % 3640 km
% target.a = (target.ra + target.rp)/2;
% target.e = (target.ra - target.rp)/(target.ra + target.rp);
% % target.ra = 9376.1; % 37183 km
% % target.rp = 250 + primary.R; % 3640 km
% % target.a = (target.ra + target.rp)/2;
% % target.e = (target.ra - target.rp)/(target.ra + target.rp);
% 
% e = target.e;
% E = acos(1/target.e - r_int/(target.a*target.e));
% v = 2*atan2(sqrt(1+e)*tan(E/2),sqrt(1-e));




% 
% [ fpa ] = oe2fpa( e, v);
% fpa*180/pi

end % run_loopingFlybyAnalysis






























