clear
clc
close all
addpath('/Users/lukebury/Documents/MATLAB/mbin/plot')

% -------------------------------------------------------------
%%% Setup
% -------------------------------------------------------------
%%% Variable definitions for ellipsoid shape
% a = alpha (x-axis)
% b = beta  (y-axis), here b = a
% g = gamma (z-axis)

%%% Setting alpha and gamma values based on a=b and a*b*g=1
a = linspace(.1,10,100);
g = 1./(a.^2);
ag = a./g; % vector of ratios for legend purposes

%%% Setting angular momentum values to be used as independant variable
H = linspace(.1,3,length(a));

% -------------------------------------------------------------
%%% Computations and Plotting
% -------------------------------------------------------------
%%% Opening new figure
figure(1); hold all
PlotBoi2('H^2','Amended Potential',14)
plot([H(1)^2 H(end)^2],[0 0],'-k')
colors = colorScale([0 1 1; 1 0 1],length(a));

%%% Preallocating matrix to store H values that make E=0 for given a/g
H2_zero = zeros(size(ag));

%%% Looping through a/g ratios
for rr = 1:length(a)
    %%% Calculating I
    if a(rr) > g(rr)
        I = 2*acos(g(rr)/a(rr)) / sqrt(a(rr)^2-g(rr)^2);
    elseif a(rr) == g(rr)
        I = 2/a(rr);
    elseif a(rr) < g(rr)
        I = log((g(rr) + sqrt(g(rr)^2-a(rr)^2))/(g(rr) - sqrt(g(rr)^2-a(rr)^2))) / sqrt(g(rr)^2-a(rr)^2);
    end
    
    %%% Calculating self-potential
    Us = -3*I/10;
    
    %%% Calculating maximum moment of inertia
    if ag(rr) < 1      % Prolate, a < g
        Ih = (a(rr)^2 + g(rr)^2)/5;
    elseif ag(rr) > 1  % Oblate,  a > g
        
        Ih = (2/5)*a(rr)^2;
    end
    
    %%% --------------------------------------------------
    %%% For plotting E as function of H^2
    %%% --------------------------------------------------
    %%% Reinitializing amended potential
    E = zeros(size(a));
    
    %%% Calculating E values at current a/g ratio for each H value
    for hh = 1:length(H)
        E(hh) = H(hh)*H(hh)/(2*Ih) + Us;
%         if E(hh) == 0
%             
%             break
%         end
    end

    %%% Plotting E values for current a/g ratio
    figure(1)
    plot(H.^2,E,'linewidth',1.5,'color',colors(rr,:))
    
    %%% --------------------------------------------------
    %%% For plotting E as function of a/g
    %%% --------------------------------------------------
    H2_zero(rr) = -2*Ih*Us;
    
end

%%% Creating legend for figure(1)
ag_leg = cell(1,length(ag)+1);
ag_leg{1} = '0';
for kk = 1:length(ag)
    ag_leg{kk+1} = sprintf('a/g - %0.2f',ag(kk));
end
figure(1)
legend(ag_leg,'location','northwest')

figure(2)
plot(ag,H2_zero,'linewidth',2)
PlotBoi2('\alpha/\gamma','H^2 that makes E=0',14)

