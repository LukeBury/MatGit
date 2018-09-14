clear
clc
close all
% ========================================================================
%%% Problem 1
% ========================================================================
%%% Defining Orbit
Re = 6378.1363; % km
alt = 800; % km
a =  Re + alt; % km
uE = 398600.4415; % km^3/s^2
n = sqrt(uE/(a^3)); % rad/s
e = 0;
p = a*(1-e^2);
J2 = 0.0010826269;

% Sun Sync Node Rate
ssnr = 360/365.2421897; % deg/day

%%% Initializing vectors for RAANdots and corresponding degrees/day values
RAANdot = zeros(3601,1);
DegPerDay = zeros(3601,1);
kk = 1;

%%% Looping through inclinations
for i = 0:.05:180
    RAANdot(kk,1) = -3*n*(Re^2)*J2*cosd(i)/(2*(p^2)); % rad/s
    DegPerDay(kk,1) = RAANdot(kk) * (180/pi) * (3600*24); % °/day
    kk = kk + 1;
end

%%% Plotting Results
hold all
plot(0:.05:180, DegPerDay,'linewidth',2)
plot(linspace(0,180,100),linspace(ssnr,ssnr,100),'--r','linewidth',2)
plot(linspace(98.60,98.60,100),linspace(-8,8,100),':k','linewidth',2)
PlotBoi2('Inclination, [°]','Secular Drift of Node, [°/day]',18)
legend('Inclination vs Rate','Sun Synchronous Rate (0.9856)','Ideal Inclination'...
    ,'location','northwest')

dim = [.57 0 1 .5];
str = 'SS Inclination ~ 98.60°, Rate = 0.9853 °/Day';
annotation('textbox',dim,'String',str,'FitBoxToText','on');

% ========================================================================
%%% Problem 2
% ========================================================================
% ------------------------------------------------------------------------
%%% Setting Constants
% ------------------------------------------------------------------------
planets = strings(8,1);
planets(1) = 'Mercury';
planets(2) = 'Venus';
planets(3) = 'Moon';
planets(4) = 'Mars';
planets(5) = 'Jupiter';
planets(6) = 'Saturn';
planets(7) = 'Uranus';
planets(8) = 'Neptune';

%%% J2 Values
J2_Mercury = 0.00006;
J2_Venus = 0.000027;
J2_Moon = 0.0002027;
J2_Mars = 0.001964;
J2_Jupiter = 0.01475;
J2_Saturn = 0.01645;
J2_Uranus = 0.012;
J2_Neptune = 0.004;

J2s = [J2_Mercury;J2_Venus;J2_Moon;J2_Mars;J2_Jupiter;J2_Saturn;...
    J2_Uranus;J2_Neptune];

%%% Equatorial Radius
R_Mercury = 2439.0; % km
R_Venus = 6052.0; % km
R_Moon = 1738.0; % km
R_Mars = 3397.2; % km
R_Jupiter = 71492.0; % km
R_Saturn = 60268.0; % km
R_Uranus = 25559.0; % km
R_Neptune = 24764.0; % km

Rs = [R_Mercury;R_Venus;R_Moon;R_Mars;R_Jupiter;R_Saturn;R_Uranus;...
    R_Neptune];

%%% Orbital Periods (Years per earth year)
Yr_Mercury = 0.24084445; 
Yr_Venus = 0.61518257;
Yr_Moon = 1;%0.0748; 
Yr_Mars = 1.88071105; 
Yr_Jupiter = 11.856525; 
Yr_Saturn = 29.423519; 
Yr_Uranus = 83.747406; 
Yr_Neptune = 163.7232045; 

Ps = [Yr_Mercury;Yr_Venus;Yr_Moon;Yr_Mars;Yr_Jupiter;Yr_Saturn;...
    Yr_Uranus;Yr_Neptune];

%%% Gravitational Parameters
u_Mercury = 2.2032e4; % km^3/s^2 
u_Venus = 3.257e5; % km^3/s^2 
u_Moon = 4902.799; % km^3/s^2 
u_Mars = 4.305e4; % km^3/s^2 
u_Jupiter = 1.268e8; % km^3/s^2 
u_Saturn = 3.794e7; % km^3/s^2 
u_Uranus = 5.794e6; % km^3/s^2 
u_Neptune = 6.809e6; % km^3/s^2 

us = [u_Mercury;u_Venus;u_Moon;u_Mars;u_Jupiter;u_Saturn;u_Uranus;...
    u_Neptune];

%Earth day
dayE = 86400; % sec

% ------------------------------------------------------------------------
%%% Determining Average Drift Rates Necessary For SS Orbit
% ------------------------------------------------------------------------
ssnr_Mercury = ssnr / Yr_Mercury; % deg/EarthDay
ssnr_Venus = ssnr / Yr_Venus; % deg/EarthDay
ssnr_Moon = ssnr / Yr_Moon; % deg/EarthDay
ssnr_Mars = ssnr / Yr_Mars; % deg/EarthDay
ssnr_Jupiter = ssnr / Yr_Jupiter; % deg/EarthDay
ssnr_Saturn = ssnr / Yr_Saturn; % deg/EarthDay
ssnr_Uranus = ssnr / Yr_Uranus; % deg/EarthDay
ssnr_Neptune = ssnr / Yr_Neptune; % deg/EarthDay

ssnrs = [ssnr_Mercury;ssnr_Venus;ssnr_Moon;ssnr_Mars;ssnr_Jupiter;...
    ssnr_Saturn;ssnr_Uranus;ssnr_Neptune];


fprintf('2a:\n')
table(planets,ssnrs)
% table(ssnr_Mercury,ssnr_Venus,ssnr_Moon,ssnr_Mars,ssnr_Jupiter,...
% ssnr_Saturn,ssnr_Uranus,ssnr_Neptune)

% ------------------------------------------------------------------------
%%% Defining Orbits and Iterating over inclinations
% ------------------------------------------------------------------------
inclinations = 0:.05:180;
%%% Initializing vectors for RAANdots and corresponding degrees/day values
RAANdot = zeros(3601,8);
DegPerDay = zeros(3601,8);
%%% Initializing vector to store necessary inclinations for SS orbit
closestInc = zeros(8,1);

%%% Looping through planets
for planet = 1:8
    %%% Setting orbital parameters
    R = Rs(planet); % km
    a = alt + R; % km
    u = us(planet); % km^3/s^2 
    n = sqrt(u/(a^3)); % rad/s
    p = a*(1-e^2); % km
    J2 = J2s(planet);
    kk = 1; % For use with indexing
    
    %%% Looping through inclinations
    for i = 0:.05:180
        RAANdot(kk,planet) = -3*n*(R^2)*J2*cosd(i)/(2*(p^2)); % rad/s
        DegPerDay(kk,planet) = RAANdot(kk,planet) * (180/pi) * (86400); % °/Eday
        kk = kk + 1;
    end
    
    %%% Finding the closest inclination to achieving SS orbit
    [diff, index] = min(abs(DegPerDay(:,planet)-ssnrs(planet)));
    closestInc(planet) = inclinations(index);
    
    %%% If SS orbit is not possible, change value to NaN
    if closestInc(planet) == 180
        closestInc(planet) = NaN;
    end
    
    %%% Plotting Results
    figure
    hold all
    plot(0:.05:180, DegPerDay(:,planet),'linewidth',2)
    plot(linspace(0,180,100),linspace(ssnrs(planet),ssnrs(planet),100),...
        '--r','linewidth',2)
    PlotBoi2('Inclination, [°]','Secular Drift of Node, [°/EarthDay]',18)
    title(planets(planet))
    legend('Inclination','°/Day Necessary for SS Orbit',...
        'location','northwest')
    
    c = closestInc(planet); % shorthand
    
    if isnan(c) == 0 % If SS Orbit is achievable
        plot(linspace(c,c,100),linspace(min(DegPerDay(:,planet)),...
            max(DegPerDay(:,planet)),100),':k','linewidth',2)
        legend('Inclination','°/Day Necessary for SS Orbit',...
        'SS Inclination','location','northwest')
    
        dim = [.57 0 1 .5];
        str = ['SS Inclination ~ ',num2str(c),'°, Rate = ',...
            num2str(ssnrs(planet)),' °/Day'];
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        
    elseif planet ~= 3 % If it's not the moon
        dim = [.45 0 1 .55];
        str = 'SS Orbit Not Achievable';
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        
    elseif planet == 3 % If it is the moon
        dim = [.35 0 1 .55];
        str = 'Earth-Synchronous Orbit Not Achievable';
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
end

%%% Plotting all on one plot
figure 
hold all
for i = 1:8
    plot(0:.05:180, DegPerDay(:,i),'linewidth',2) 
end
PlotBoi2('Inclination, [°]','Secular Drift of Node, [°/EarthDay]',18)
legend('Me','Ve','Mo','Ma','Ju','Sa','Ur','Ne','location','southeast')
    
fprintf('\n\n2b:\n')
fprintf('Inclinations for SS orbit:\n')
table(planets,closestInc)





