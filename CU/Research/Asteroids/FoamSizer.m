clear
clc
close all

PlotsOn = 1;
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
%%% Outer Radius of sphere
d = 35; % cm

%%% Thickness of foam
t = 3.81; % cm

if rem(d,t) ~= 0
    warning('Diameter not divisible by foam thickness')
end

% ------------------------------------------------------------------------
%%% Calculations
% ------------------------------------------------------------------------
r = d/2; % Radius of sphere, cm

V = (4/3)*pi*(r^3); % Volume of sphere, cm^3

nP = ceil(2*r/t); % Number of foam pieces

foams = zeros(nP,3);
dBot = 0; % Bottom of current piece
dTop = t; % Top of current piece
totalFoamArea = 0;
perfectMiddle = 0;
for k = 1:nP % For each foam piece
    if dTop+1e-6 < 2*r
        ny_Top = r - dTop;
        nx_Top = sqrt(abs(r^2 - ny_Top^2));
        nLength_Top = 2*nx_Top;
    else
        nLength_Top = 0;
    end
    
    ny_Bot = r - dBot;
    nx_Bot = sqrt(abs(r^2 - ny_Bot^2));
    nLength_Bot = 2*nx_Bot;
    
    %%% Assigning top/bottom lengths of foam piece
    foams(k,1) = nLength_Bot;
    foams(k,2) = nLength_Top;
    
    %%% Assigning cut-to length
    foams(k,3) = max(foams(k,1:2));
    
    %%% Treating each circle as a square
    totalFoamArea = totalFoamArea + (max(foams(k,1:2))/2)^2;
    
    dBot =  dTop;
    dTop = dTop + t;
    
    if dBot+1e-6 > 2*r
        break
    end
    
end

%%% If odd # of pieces or if diameter not divisible by thickness, cut middle piece to sphere diameter
if rem(nP,2) ~= 0 || rem(d,t) ~= 0
    foams(round(nP/2),3) = d;
end
% ------------------------------------------------------------------------
%%% Printing Data
% ------------------------------------------------------------------------
fprintf('Foam Sheet Lengths and Cut Sizes:\n\n')
fprintf('For Diameter = %3.1f cm\n',d)
fprintf('Sheet # | Bottom D | Top D | Cut To\n')
for k = 1:nP
   fprintf('   %02d   |   %05.2f  | %05.2f | %05.2f\n',k,foams(k,1),foams(k,2),foams(k,3))
end
fprintf('\n**NOTE** All sizes in cm\n')
fprintf('\nPerfect Foam Area Needed: %3.0f cm^2\n',totalFoamArea)

% ------------------------------------------------------------------------
%%% Plotting
% ------------------------------------------------------------------------
if PlotsOn == 1
%%% Preparing to Plot
theta = 0:.01:2*pi;

%%% Plot
figure; hold all
plot(r*cos(theta),r*sin(theta)+r,'k','linewidth',2)
PlotBoi2('[cm]','[cm]',16)
axis equal

for k = 1:nP-1
    plot([-foams(k,2)/2 foams(k,2)/2], [t*k, t*k],'r','linewidth',2);
end
end