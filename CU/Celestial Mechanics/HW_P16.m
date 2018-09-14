clear
clc
close all
addpath('/Users/lukebury/Documents/MATLAB/mbin/plot')

% -------------------------------------------------------------
%%% Setup
% -------------------------------------------------------------
%%% Integrator options
tol = 1e-7;
options = odeset('RelTol',tol,'AbsTol',tol);

% Setting timing options
t0 = 0;
dt = .1;
tf = 50;
tVec = t0:dt:tf; % (s)

%%% Setting parameters
G = 6.67408 * 10-11; % (m^3 kg^-1 s^-2)
dens = 1; % (kg m^-3)
R1R2 = 0:.001:1; % R1/R2
R2 = .5; % m
R1 = R1R2.*R2; % m

%%% Preallocating fission energy vector
Ef = zeros(size(R1R2));

%%% Looping through ratios
for rr = 1:length(R1R2)
    % Seperating distance
    r12 = R1(rr) + R2; % (m)
    
    % Masses
    m1 = 4*dens*pi*(R1(rr)^3)/3; % (kg)
    m2 = 4*dens*pi*(R2^3)/3; % (kg)
    
    % Distance from m1 to COM
    r1c = (R1(rr)+R2)*m2/(m1+m2); % (m)
    
    % Fission spin rate
    w = sqrt((4/3)*pi*G*dens*(R1(rr)^3+R2^3)/((R1(rr)+R2)^3)); % rad/s
    
    % Tangential velocity of m1 wrt COM
    v1 = w*r1c; % (m s^-1)
    
    % Calculate I
    I = m1*m2*r12*r12/(m1+m2) + (2/5)*(m1*R1(rr)*R1(rr) + m2*R2*R2); % (kg m^2)
    
    % Calculate total H
    r2c = R1(rr)+R2-r1c; % m2 to COM
    v2 = w*r2c;
    H = m1*r1c*v1 + m2*r2c*v2 + (2/5)*(m1*R1(rr)*R1(rr) + m2*R2*R2)*w; % (kg m^2 s^-1)

    % Calculate U
    U = -G*m1*m2/r12; % (kg m^2 s^-2)

    % Calculate free energy
    Ef(rr) = H*H/(2*I) + U; % (kg m^2 s^-2)
    
    if Ef(rr) < 0 && Ef(rr-1) > 0
        crossPoint = R1R2(rr);
    end

end

figure; hold all
plot(R1R2,Ef,'-','linewidth',2)
plot([0 1],[0 0], 'r--','linewidth',1.5)
plot([crossPoint crossPoint],[min(Ef) max(Ef)],'r--','linewidth',1.5)
PlotBoi2('R1/R2','Free Energy at Fission (kg m^2/s^2)',14)
str = sprintf('Crossing at R1/R2 = %0.2f',crossPoint);
t = annotation('textbox',[.62 .3 .3 .26],'String',str,'FitBoxToText','on');










