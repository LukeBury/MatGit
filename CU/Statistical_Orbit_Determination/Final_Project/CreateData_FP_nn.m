clear
clc
close all
addpath('../../Research/Europa_Hopper_Analysis/ProjectBin')
addpath('../../bin')
% ------------------------------------------------------------------------
%%% Plot/Run Options; 0 = no, 1 = yes
% ------------------------------------------------------------------------
plotEurCI         = 0;
plotEurCR         = 0;
plotAltitudeEcc   = 0;
plotVisibility    = 0;
plotRangeData     = 0;
% ------------------------------------------------------------------------
%%% System Parameters
% ------------------------------------------------------------------------
d2r = pi/180;
aEur = 671100; % semimajor axis 
wEur = [0,0,2.047200349303344e-05]; % rad/s

%%% Gravitational Parameters
u1 = 126672520; % Jupiter, km^3 / s^2
u2 = 3203.413216; % Europa, km^3 / s^2

% Body Radii
rad1 = 71492; % radius of Jupiter
rad2 = 1560.8; % radius of Europa

%%% Initial Julian Date
[JD0] = JD(2025,3,5,0,0);

% Earth
rE = 6378.1363; % km
wDay = [0,0,2*pi/86400]; % rad/s
% ------------------------------------------------------------------------
%%% Defining Satellite State
% ------------------------------------------------------------------------
%%% Initial Particle Position
% Position of satellite wrt Europa
a = rad2+25;
% e = 0.02845; %w/ i = 90 crashes near 19 hrs
e = .005;
% e = 0.000;
i = 90 *d2r;
raan = 0 *d2r;
w = 0 *d2r;
ta = 90 *d2r;
[rH0_EurCI, vH0_EurCI] = OE2ECI(a, e, i, raan, w, ta, u2);
rH0_EurCI = rH0_EurCI';
vH0_EurCI = vH0_EurCI';

save('OrbitalInfo_nn.mat','a','e','i','raan','w','ta')

% ------------------------------------------------------------------------
%%% Propagating State
% ------------------------------------------------------------------------
%%% Europa period (days)
pEur = 3.551181;
%%% Setting time vector
day = 86400; % sec
ti = 0;
dt = 60;
tf = 1*day;
time = (ti:dt:tf);

%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('Events',@impact_FP_nn,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (ECEF)
X0 = [rH0_EurCI, vH0_EurCI]'; % km, km/s

%%% Propagating the State 
[Times,states_EurCI_truth] = ode45(@integrator_FPdata_nn,time,X0,options,rad2,u1,u2,aEur,wEur);

% ------------------------------------------------------------------------
%%% Other Important Planetary States
% ------------------------------------------------------------------------
%%% Creating Other Important States
states_EJ = zeros(length(Times),6);
states_EaSun = zeros(length(Times),6);
states_JSun =zeros(length(Times),6);
eccentricities = zeros(length(Times),1);
inclinations = zeros(length(Times),1);

%%% Calculating
for k = 1:length(Times)
    %%% Current Rotation Angle
    th = Times(k)*wEur(3); % km
    
    %%% Europa relative to Jupiter Inertial
    states_EJ(k,1:3) = R3([aEur, 0, 0],th); % km
    states_EJ(k,4:6) = R3(cross(wEur,states_EJ(1,1:3)),th); % km/s
    
    %%% Jupiter relative to Sun inertial
    [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,5,'EME2000'); % km, km/s, km^3/s^2
    states_JSun(k,1:6) = [r',v'];
    
    %%% Earth relative to Sun inertial
    [r, v, mu_p] = Ephem(JD0 + Times(k)/86400,3,'EME2000'); % km, km/s, km^3/s^2
    states_EaSun(k,1:6) = [r',v'];
    
    %%% Eccentricities
    [a e i raan w ta] = ECI2OE(states_EurCI_truth(k,1:3)',states_EurCI_truth(k,4:6)',u2);
    eccentricities(k) = e;
    inclinations(k) = i*180/pi;
end

% ------------------------------------------------------------------------
%%% Station States
% ------------------------------------------------------------------------
%%% Initial lat/lon/alts of Stations
latLonAlt = [-35.398333, 148.981944, 0.691750;
    40.427222, -4.250556, 0.834539;
    35.247164, 243.205, 1.07114904]; % deg | deg | km

%%% Calculating Intial Positions and Velocities of Stations in ECI
% statStates [6x3xTimes] = [[6x1]-Station 1 State | [6x1]-Station 2 State | [6x1]-Station 3 State] ... 3rd dimension is time
stations_ECI = zeros(6,3,length(Times));

for k = 1:3
    %%% Acquiring ECEF/ECI station position
    [r] = latlon2surfECEF(latLonAlt(k,1), latLonAlt(k,2), latLonAlt(k,3) + rE); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Correct way of handling altitude?

    %%% Storing initial station positions and velocities
    stations_ECI(1:3,k,1) = r'; % km
    stations_ECI(4:6,k,1) = cross(wDay',r'); % km/s
end; clear r

%%% Propagating Station States in ECI
for k = 2:length(Times)
    %%% Calculating current rotation angle
    th = wDay(3)*Times(k); % rad 
    for s = 1:3
        %%% ECI Position
        stations_ECI(1:3,s,k) = R3(stations_ECI(1:3,s,1),th); % km
        %%% ECI Velocity
        stations_ECI(4:6,s,k) = R3(stations_ECI(4:6,s,1),th); % km/s
    end
end

% ------------------------------------------------------------------------
%%% ECI & ECEF States
% ------------------------------------------------------------------------
states_ECI = zeros(length(Times),6);
states_ECEF = zeros(length(Times),6);
jupiter_ECI = zeros(length(Times),3);
europa_ECI = zeros(length(Times),3);

for k = 1:length(Times)
    %%% Calculating current rotation angle
    th = wDay(3)*Times(k); % rad 
    
    states_ECI(k,1:3) = states_EurCI_truth(k,1:3) + states_EJ(k,1:3) + states_JSun(k,1:3) - states_EaSun(k,1:3); % km
    states_ECI(k,4:6) = states_EurCI_truth(k,4:6) + states_EJ(k,4:6) + states_JSun(k,4:6) - states_EaSun(k,4:6); % km/s
    
    states_ECEF(k,1:3) = R3(states_ECI(k,1:3),-th); % km
    states_ECEF(k,4:6) = R3(states_ECI(k,4:6),-th) - cross(wDay,states_ECI(k,1:3)); % km/s
    
    jupiter_ECI(k,:) = states_JSun(k,1:3) - states_EaSun(k,1:3); % km
    europa_ECI(k,:) = states_EJ(k,1:3) + states_JSun(k,1:3) - states_EaSun(k,1:3); % km
end

% ------------------------------------------------------------------------
%%% Other Data (including ECI & ECEF)
% ------------------------------------------------------------------------
europaAltitudes = zeros(length(Times),1);
states_EurCR = zeros(length(Times),6);

%%% Altitude from Europa Surface
for k = 1:length(Times)
    europaAltitudes(k) = norm(states_EurCI_truth(k,1:3)) - rad2; % km
end

for k = 1:length(Times)
    %%% Calculating current rotation angle
    th = wEur(3)*Times(k); % rad 
    
    states_EurCR(k,1:3) = R3(states_EurCI_truth(k,1:3),-th); % km
    states_EurCR(k,4:6) = R3(states_EurCI_truth(k,4:6),-th) - cross(wEur,states_EurCI_truth(k,1:3)); % km/s
end

% ------------------------------------------------------------------------
%%% Preparing for Data Generation
% ------------------------------------------------------------------------
%%% Setting Noise
pError = (5e-3)*2; % km
dpError = 5e-7; % km/s

%%% Setting Station Elevation Mask
elevMask = 10*d2r; % rad

%%% Initializing Data Vector
data = zeros(1,7);

% ------------------------------------------------------------------------
%%% Checking Visibility, Creating Data, and Adding Noise
% ------------------------------------------------------------------------
shadowAngles_J = zeros(length(Times),1);
shadowAngles_E = zeros(length(Times),1);
satAngles_J = zeros(length(Times),1);
satAngles_E = zeros(length(Times),1);
behindJ = NaN(length(Times),1);
behindE = NaN(length(Times),1);
ranges_truth = NaN(length(Times),3); % [range, time, station]
rangeRates_truth = NaN(length(Times),3);


dk = 1;       % Overall measurement number
for k = 1:length(Times)
    %%% Checking Visibility
    % Angle between center and rim of Jupiter/Europa in ECI
    shadowAngles_J(k) = atan(rad1/norm([jupiter_ECI(k,1:2), 0])); % rad
    shadowAngles_E(k) = atan(rad2/norm(europa_ECI(k,1:3))); % rad
    % Angle between center of Jupiter and satellite vector in ECI
    satAngles_J(k) = atan2(norm(cross([jupiter_ECI(k,1:2), 0],[states_ECI(k,1:2), 0])),dot([jupiter_ECI(k,1:2), 0],[states_ECI(k,1:2), 0])); % rad
    satAngles_E(k) = atan2(norm(cross(europa_ECI(k,1:3),states_ECI(k,1:3))),dot(europa_ECI(k,1:3),states_ECI(k,1:3))); % rad

    if satAngles_J(k) > shadowAngles_J(k) || norm([states_ECI(k,1:2), 0]) < norm([jupiter_ECI(k,1:2), 0]) % If not blocked by Jupiter
        if satAngles_E(k) > shadowAngles_E(k) || norm(states_ECI(k,1:3)) < norm(europa_ECI(k,1:3)) % If not blocked by Europa
            %%% Calculating elevations and storing ranges and range rates
            [az, el1, p] = ECEF2az_el_p(states_ECEF(k,1:3)', latLonAlt(1,1)*d2r, latLonAlt(1,2)*d2r, latLonAlt(1,3));
            [az, el2, p] = ECEF2az_el_p(states_ECEF(k,1:3)', latLonAlt(2,1)*d2r, latLonAlt(2,2)*d2r, latLonAlt(2,3));
            [az, el3, p] = ECEF2az_el_p(states_ECEF(k,1:3)', latLonAlt(3,1)*d2r, latLonAlt(3,2)*d2r, latLonAlt(3,3));
            if el1 > elevMask
                data(dk,1) = Times(k); % sec
%                 ranges_truth(k,1) = sqrt(dot(States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI(1:3,1,k)',States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI(1:3,1,k)')) + mvnrnd(zeros(1,1),pError^2,1); % km ;
%                 rangeRates_truth(k,1) = dot(States_EurCI_truth(k,1:3) + europa_JCI(k,1:3) + jupiter_HC(k,1:3) - earth_HC(k,1:3) - stations_ECI(1:3,1,k)',States_EurCI_truth(k,4:6) + europa_JCI(k,4:6) + jupiter_HC(k,4:6) - earth_HC(k,4:6) - stations_ECI(4:6,1,k)')/ranges_truth(k,1) + mvnrnd(zeros(1,1),dpError^2,1); % km/s ;
                ranges_truth(k,1) = sqrt(dot(states_ECI(k,1:3) - stations_ECI(1:3,1,k)',states_ECI(k,1:3) - stations_ECI(1:3,1,k)')) + mvnrnd(zeros(1,1),pError^2,1); % km ;
                rangeRates_truth(k,1) = dot(states_ECI(k,1:3) - stations_ECI(1:3,1,k)',states_ECI(k,4:6) - stations_ECI(4:6,1,k)')/ranges_truth(k,1) + mvnrnd(zeros(1,1),dpError^2,1); % km/s ;
                ranges_truth(k,2) = Times(k);
                rangeRates_truth(k,2) = Times(k);
                ranges_truth(k,3) = 1;
                rangeRates_truth(k,3) = 1;
                data(dk,2) = ranges_truth(k,1); % km
                data(dk,5) = rangeRates_truth(k,1); % km/s
                dk = dk + 1;

            elseif el2 > elevMask
                data(dk,1) = Times(k);
                ranges_truth(k,1) = sqrt(dot(states_ECI(k,1:3) - stations_ECI(1:3,2,k)',states_ECI(k,1:3) - stations_ECI(1:3,2,k)'))  + mvnrnd(zeros(1,1),pError^2,1); % km ;
                rangeRates_truth(k,1) = dot(states_ECI(k,1:3) - stations_ECI(1:3,2,k)',states_ECI(k,4:6) - stations_ECI(4:6,2,k)')/ranges_truth(k,1) + mvnrnd(zeros(1,1),dpError^2,1); % km/s;
                ranges_truth(k,2) = Times(k);
                rangeRates_truth(k,2) = Times(k);
                ranges_truth(k,3) = 2;
                rangeRates_truth(k,3) = 2;
                data(dk,3) = ranges_truth(k,1); % km
                data(dk,6) = rangeRates_truth(k,1); % km/s
                dk = dk + 1;
                
            elseif el3 > elevMask
                data(dk,1) = Times(k);
                ranges_truth(k,1) = sqrt(dot(states_ECI(k,1:3) - stations_ECI(1:3,3,k)',states_ECI(k,1:3) - stations_ECI(1:3,3,k)')) + mvnrnd(zeros(1,1),pError^2,1); % km ;
                rangeRates_truth(k,1) = dot(states_ECI(k,1:3) - stations_ECI(1:3,3,k)',states_ECI(k,4:6) - stations_ECI(4:6,3,k)')/ranges_truth(k,1) + mvnrnd(zeros(1,1),dpError^2,1); % km/s;
                ranges_truth(k,2) = Times(k);
                rangeRates_truth(k,2) = Times(k);
                ranges_truth(k,3) = 3;
                rangeRates_truth(k,3) = 3;
                data(dk,4) = ranges_truth(k,1); % km
                data(dk,7) = rangeRates_truth(k,1); % km/s
                dk = dk + 1;

            end
        end
    end
end

for k = 1:length(Times)
    if satAngles_J(k) < shadowAngles_J(k) && norm([states_ECI(k,1:2), 0]) > norm([jupiter_ECI(k,1:2), 0]) % If blocked by Jupiter
        behindJ(k) = satAngles_J(k);
    end
    
    if satAngles_E(k) < shadowAngles_E(k) && norm(states_ECI(k,1:3)) > norm(europa_ECI(k,1:3)) % If blocked by Europa
        behindE(k) = satAngles_E(k);
    end
end

%%% Correcting/Configuring Data
firstTime = data(1,1);
data(data == 0) = NaN;
data(1,1) = firstTime;

% ------------------------------------------------------------------------
%%% Writing/Saving Data
% ------------------------------------------------------------------------
%%% Correcting/Configuring Data
firstTime = data(1,1);
data(data == 0) = NaN;
data(1,1) = firstTime;

%%% Saving EurCI Truth
save('states_EurCI_truth.mat','states_EurCI_truth')
%%% Writing Data
dlmwrite('finalProjectData_ECI2.txt',data,'precision',15);

% Old data
dataSet1 = csvread('finalProjectData_ECI.txt');

%%% Collapsing data matrix
data1 = zeros(size(dataSet1,1),4);
for k = 1:size(dataSet1,1)
    if isnan(dataSet1(k,2)) == 0
        data1(k,:) = [1,dataSet1(k,1),dataSet1(k,2),dataSet1(k,5)];
    elseif isnan(dataSet1(k,3)) == 0
        data1(k,:) = [2,dataSet1(k,1),dataSet1(k,3),dataSet1(k,6)];
    elseif isnan(dataSet1(k,4)) == 0
        data1(k,:) = [3,dataSet1(k,1),dataSet1(k,4),dataSet1(k,7)];
    end
end

% New data
dataSet2 = csvread('finalProjectData_ECI2.txt');

%%% Collapsing data matrix
data2 = zeros(size(dataSet2,1),4);
for k = 1:size(dataSet2,1)
    if isnan(dataSet2(k,2)) == 0
        data2(k,:) = [1,dataSet2(k,1),dataSet2(k,2),dataSet2(k,5)];
    elseif isnan(dataSet2(k,3)) == 0
        data2(k,:) = [2,dataSet2(k,1),dataSet2(k,3),dataSet2(k,6)];
    elseif isnan(dataSet2(k,4)) == 0
        data2(k,:) = [3,dataSet2(k,1),dataSet2(k,4),dataSet2(k,7)];
    end
end



% ------------------------------------------------------------------------
%%% Jacobi Constant
% ------------------------------------------------------------------------
% JC = zeros(length(Times),1);
% 
% for k = 1:length(Times)    
%     th = Times(k)*wEur(3);
% %     rBC = R3(rE0_JCI.*(uE/uJ), th); % Barycenter position (JCI)
% %     vBC = cross(wEur,rBC); % Barycenter velocity (JCI) 
% %     x = rH_ECEF_I(k,1)+rE0_JCI(1)-rBC(1);
% %     y = rH_ECEF_I(k,2);
% %     r = [x, y, 0];
% %     r1 = norm(StatesI(k,1:3));
% %     r2 = norm(rH_ECEF_I(k,:));
% %     
% %     vJC = StatesI(k,4:6) - vBC - cross(wE',r);
% %     
% %     JC_I(k) = (nE^2)*(x^2 + y^2) + 2*(uJ/r1 + uE/r2) - norm(vJC)^2;
% 
%     rBC = R3([aEur*(u2/(u1+u2)),0,0], th); % Barycenter position (JCI)
%     vBC = cross(wEur,rBC); % Barycenter velocity (JCI) 
%     x = states_EurCR(k,1)+aEur-rBC(1);
%     y = states_EurCR(k,2);
%     r = [x, y, 0];
%     r1 = norm([aEur,0,0] + states_EurCR(k,1:3));
%     r2 = norm(states_EurCR(k,1:3));
%     
%     vJC = states_EurCI_truth(k,4:6) - vBC - cross(wEur,r);
%     
%     JC(k) = (wEur(3)^2)*(x^2 + y^2) + 2*(u1/r1 + u2/r2) - norm(vJC)^2;
% 
% end
% figure
% plot(Times,JC)
% title('JC')
% 
% clear x y rBC vBC vJC r1 r2 




% ------------------------------------------------------------------------
%%% Plots
% ------------------------------------------------------------------------
if plotEurCI == 1
    figure; hold all
    plot3(states_EurCI_truth(:,1),states_EurCI_truth(:,2),states_EurCI_truth(:,3),'m','linewidth',1.5)
    bodySurface3(rad2,[0,0,0],[0,1,1])
    PlotBoi3('X','Y','Z',16)
    axis equal
end

if plotEurCR == 1
    figure; hold all
    plot3(states_EurCR(:,1),states_EurCR(:,2),states_EurCR(:,3),'m','linewidth',1.5)
    bodySurface3(rad2,[0,0,0],[0,1,1])
    PlotBoi3('X','Y','Z',16)
    axis equal
end

if plotAltitudeEcc == 1
    figure; hold all
    plot(Times./3600,europaAltitudes,'m','linewidth',1.5)
    PlotBoi2('Time, hr','Europa Altitude, km',16)
    
    figure; hold all
    plot(Times./3600, eccentricities,'b','linewidth',1.5)
    PlotBoi2('Time, hr','Eccentricity',16)
    
    figure; hold all
    plot(Times./3600, inclinations,'b','linewidth',1.5)
    PlotBoi2('Time, hr','Inclination, °',16)
    
end

%%% Plotting Visibility 
if plotVisibility == 1
    figure; 
    subplot(2,1,1); hold all
    plot(Times./3600,shadowAngles_J,'k')
    plot(Times./3600,satAngles_J,'m')
    plot(Times./3600,behindJ,'r','linewidth',2)
    PlotBoi2('Times, hr','Visibility Angles, rads',16)
    title('Jupiter Visibility Angles in X-Y axis ONLY')
    legend('shadowAngle','satAngle','Behind Jupiter')
    
    subplot(2,1,2); hold all
    plot(Times./3600,shadowAngles_E,'k')
    plot(Times./3600,satAngles_E,'m')
    plot(Times./3600,behindE,'r','linewidth',2)
    PlotBoi2('Times, hr','Visibility Angles, rads',16)
    title('Europa Visibility Angles')
    legend('shadowAngle','satAngle','Behind Europa')
end

%%% Plotting Range Data
if plotRangeData == 1
    figure
    subplot(2,1,1); hold all
    plot(data(:,1)./3600,data(:,2),'.m')
    plot(data(:,1)./3600,data(:,3),'.r')
    plot(data(:,1)./3600,data(:,4),'.b')
    legend('Canberra','Madrid','Goldstone')
    PlotBoi2('','Range Data, km',16)
    
    subplot(2,1,2); hold all
    plot(data(:,1)./3600,data(:,5),'.r')
    plot(data(:,1)./3600,data(:,6),'.b')
    plot(data(:,1)./3600,data(:,7),'.m')
    legend('Canberra','Madrid','Goldstone')
    PlotBoi2('Time, hr','Range-Rate Data, km/s',16)
end





