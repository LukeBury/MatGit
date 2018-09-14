clear
clc
close all

%% ------------------------------------------------------------------------
%%%                           Acquiring Data
%%%------------------------------------------------------------------------

t1   = xlsread('thin');
t2    = xlsread('medium');
t3  = xlsread('thick');

%% ------------------------------------------------------------------------
%%%                        Parsing Acquired Data
%%%------------------------------------------------------------------------
close all

t1L   = t1(1:end,1); % Load
t1d   = t1(1:end,2); % Extension
t1t   = t1(1:end,3); % Time
t1max = max(t1L);    % Max Load

t2L = t2(1:end,1); % Load
t2d = t2(1:end,2); % Extension
t2t = t2(1:end,3); % Time
t2max = max(t2L);    % Max Load

t3L = t3(1:end,1); % Load
t3d = t3(1:end,2); % Extension
t3t = t3(1:end,3); % Time
t3max = max(t3L);    % Max Load

%%%------------------------------------------------------------------------
%%%                           Measured Data
%%%------------------------------------------------------------------------

b1 = 0.5015; % Thickness, Inches
b2 = 0.2510; % Thickness, Inches
b3 = 0.1255; % Thickness, Inches

w1 = 2.5035; % Length, Inches
w2 = 2.5030; % Length, Inches
w3 = 2.5010; % Length, Inches

a1 = 0.988;  % Inches
a2 = 0.935;  % Inches
a3 = 1.133;  % Inches

x1 = a1/w1;
x2 = a2/w2;
x3 = a3/w3;

kc1 = (t1max/(b1*sqrt(w1)))*(29.6*(x1^(1/2))-185.5*(x1^(3/2))+655.7*(x1^(5/2))-1017*(x1^(7/2))+639*(x1^(9/2)));
kc2 = (t2max/(b2*sqrt(w2)))*(29.6*(x2^(1/2))-185.5*(x2^(3/2))+655.7*(x2^(5/2))-1017*(x2^(7/2))+639*(x2^(9/2)));
kc3 = (t3max/(b3*sqrt(w3)))*(29.6*(x3^(1/2))-185.5*(x3^(3/2))+655.7*(x3^(5/2))-1017*(x3^(7/2))+639*(x3^(9/2)));

ans1 = [kc1 kc2 kc3]

% Theoretical kc line

b  = [.1:.005:.55];
Tkc1 = (t1max./(b.*sqrt(w1))).*(29.6*(x1^(1/2))-185.5*(x1^(3/2))+655.7*(x1^(5/2))-1017*(x1^(7/2))+639*(x1^(9/2)));
Tkc2 = (t2max./(b.*sqrt(w2))).*(29.6*(x2^(1/2))-185.5*(x2^(3/2))+655.7*(x2^(5/2))-1017*(x2^(7/2))+639*(x2^(9/2)));
Tkc3 = (t3max./(b.*sqrt(w3))).*(29.6*(x3^(1/2))-185.5*(x3^(3/2))+655.7*(x3^(5/2))-1017*(x3^(7/2))+639*(x3^(9/2)));
%%%------------------------------------------------------------------------
%%%                             Plotting
%%%------------------------------------------------------------------------

figure(1)                      % Displacement1 vs Loading1
plot(t1d,t1L,'linewidth',2)
PlotBoi('Delta','Loading')

figure(2)                      % Displacement2 vs Loading2
plot(t2d,t2L,'linewidth',2)
PlotBoi('Delta','Loading')

figure(3)                      % Displacement3 vs Loading3
plot(t3d,t3L,'linewidth',2)
PlotBoi('Delta','Loading')

%%
figure(4)                      % Kc vs Thickness
hold all
plot(b1,kc1,'x','markersize',10,'linewidth',3)
plot(b2,kc2,'x','markersize',10,'linewidth',3)
plot(b3,kc3,'x','markersize',10,'linewidth',3)
plot(b,Tkc1)
plot(b,Tkc2)
plot(b,Tkc3)
PlotBoi('Thickness, in','Kc')

















