clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

run_P1  = 0;
run_P2  = 0;
run_P3  = 1;
run_P4  = 0;
run_P57 = 0;

% ------------------------------------------------------------------------
%%% Problem 1
% ------------------------------------------------------------------------
if run_P1 == 1
syms m a l om real
A = -(om^2)*((.5*m*a^2 - (1/12)*m*(3*a^2 + l^2))/((1/12)*m*(3*a^2 + l^2)))
A = simplify(A);
rhs = simplify(.5*A)
    
    



end % run_P1

% ------------------------------------------------------------------------
%%% Problem 2
% ------------------------------------------------------------------------
if run_P2 == 1
syms I w c H real

cs = (H^2-4*I^2*w^2)/(I^2);
eqn = acos((2*I*w^2 + I*cs) / ((4*I^2*w^2+I^2*cs)*(w^2 + cs))^(1/2))
d_eqn = diff(eqn,w)
answer = solve(d_eqn==0,w)

w = (6^(1/2)*H)/(6*I);
final = subs(eqn)
simplify(final)

end % run_P1


% ------------------------------------------------------------------------
%%% Problem 3
% ------------------------------------------------------------------------
if run_P3 == 1
%%% Setting inertia tensor
I1 = 125; I2 = 100; I3 = 75; % kg * m^2
I = blkdiag(125, 100, 75);   % kg * m^2

%%% Setting integration time
t0 = 0; dt = .01; tf = 60;
time = dt:dt:tf;

%%% Angular velocity magnitude
w = 1; % rad/s
w_n = 0.1;

%%% Preallocating
titles = {'Nominal Axis-1 case';...
          'Perturbed Axis-1 case';...
          'Nominal Axis-2 case';...
          'Perturbed Axis-2 case';...
          'Nominal Axis-3 case';...
          'Perturbed Axis-3 case'};
      
sigs = {zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1)};

ws = {zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1)};

%%% creating cases
w0cases = {[w;0;0];...
           [w;0;w_n];...
           [0;w;0];...
           [0;w;w_n];...
           [0;0;w];...
           [w_n;0;w]};

for ci = 1:size(w0cases,1)
%%% Choosing ode45 tolerance
tol = 1e-9;

%%% Setting integrator options
options = odeset('Events',@int_mrpW_event,'RelTol',tol,'AbsTol',tol);

%%% IC
X0 = [0;0;0;w0cases{ci}];

%%% Propagating the States
[time_out, mrpW] = ode45(@int_mrpW, time, X0, options, I);
c = 0;
while time_out(end) ~= time(end) % for switching to shadow sets and restarting integration
    c = c + 1;
    timeIndex = find(round(time,3)==round(time_out(end-1),3));
    mrpW(timeIndex,1:3) = -mrpW(timeIndex,1:3)./(norm(mrpW(timeIndex,1:3))^2);
    [time_n, mrpW_n] = ode45(@int_mrpW, [time(timeIndex):dt:tf], mrpW(timeIndex,:)', options, I);
    mrpW = [mrpW; mrpW_n];
    time_out = [time_out;time_n];
end

%%% Removing duplicate entries
[~,t_idx] = unique(time_out);
time_out = time_out(t_idx);
mrpW = mrpW(t_idx,:);

%%% Storing results
sigs_out = mrpW(:,1:3)';
ws_out   = mrpW(:,4:6)';

sigs{ci} = sigs_out;
ws{ci} = ws_out;

figure; hold all
plot(time_out,ws{ci}(1,:),'.','linewidth',3)
plot(time_out,ws{ci}(2,:),'.','linewidth',3)
plot(time_out,ws{ci}(3,:),'.','linewidth',3)
PlotBoi2('Time, sec','Angular Velocity, rad/s',14)
[legh,objh] = legend('\omega1', '\omega2', '\omega3');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-');
title([titles{ci},' - Angular velocity'])

figure; hold all
plot(time_out,sigs{ci}(1,:),'.','linewidth',2)
plot(time_out,sigs{ci}(2,:),'.','linewidth',2)
plot(time_out,sigs{ci}(3,:),'.','linewidth',2)
PlotBoi2('Time, sec','MRP Value',14)
[legh,objh] = legend('\sigma1', '\sigma2', '\sigma3');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-');
title([titles{ci},' - MRPs'])
ylim([-1 1])

end % for ci

distFig('Screen','Secondary') % 'Main' or 'Secondary'

end % run_P3

% ------------------------------------------------------------------------
%%% Problem 4
% ------------------------------------------------------------------------
if run_P4 == 1
%%% Setting inertia tensor
I1 = 125; I2 = 100; I3 = 75; % kg * m^2
I = blkdiag(125, 100, 75);   % kg * m^2

%%% Setting integration time
t0 = 0; dt = .001; tf = 60;
time = dt:dt:tf;

%%% Angular velocity magnitude
w = 1; % rad/s
w_n = 0.1;

%%% Preallocating
sigs = {zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1)};
    
ws = {zeros(3,length(time)+1);...
      zeros(3,length(time)+1);...
      zeros(3,length(time)+1)};

Ks = {zeros(3,length(time)+1);...
      zeros(3,length(time)+1);...
      zeros(3,length(time)+1)};
    
KsVf = {zeros(3,length(time)+1);...
        zeros(3,length(time)+1);...
        zeros(3,length(time)+1)};

titles = {'Perturbed Axis-1 case';...
          'Perturbed Axis-2 case';...
          'Perturbed Axis-3 case'};

%%% creating cases
w0cases = {[w;0;w_n];...
         [0;w;w_n];...
         [w_n;0;w]};
     
% ---------------------------------------------------------
% Integrating
% ---------------------------------------------------------
%%% Choosing ode45 tolerance
tol = 1e-9;

%%% Setting integrator options
options = odeset('Events',@int_mrpW_event,'RelTol',tol,'AbsTol',tol);

for ci = 1:3
%%% IC
X0 = [0;0;0;w0cases{ci}];

%%% Propagating the States
[time_out, mrpW] = ode45(@int_mrpW, time, X0, options, I);
c = 0;
while time_out(end) ~= time(end) % for switching to shadow sets and restarting integration
    c = c + 1;
    timeIndex = find(round(time,3)==round(time_out(end-1),3));
    mrpW(timeIndex,1:3) = -mrpW(timeIndex,1:3)./(norm(mrpW(timeIndex,1:3))^2);
    [time_n, mrpW_n] = ode45(@int_mrpW, [time(timeIndex):dt:tf], mrpW(timeIndex,:)', options, I);
    mrpW = [mrpW; mrpW_n];
    time_out = [time_out;time_n];
end
[~,t_idx] = unique(time_out);
time_out = time_out(t_idx);
mrpW = mrpW(t_idx,:);

sigs_out = mrpW(:,1:3)';
ws_out   = mrpW(:,4:6)';

sigs{ci,1} = sigs_out;
ws{ci,1} = ws_out;

% ---------------------------------------------------------
% Duffing Equation work
% ---------------------------------------------------------    
for kk = 1:size(mrpW,1)
    % ---------------------------------------------------------
    % Duffing Equation work
    % ---------------------------------------------------------
    
    % ---------------------
    % Computing T and H
    % ---------------------
    T = (1/2) * ws_out(:,kk)' * I * ws_out(:,kk);
    H = I * ws_out(:,kk);
    H = sqrt(H'*H); % taking norm

    % ---------------------
    % Computing As and Bs (slide 42)
    % ---------------------
    A1 = ((I1-I2)*(2*I3*T-H^2) + (I1-I3)*(2*I2*T-H^2))/(I1*I2*I3);
    A2 = ((I2-I3)*(2*I1*T-H^2) + (I2-I1)*(2*I3*T-H^2))/(I1*I2*I3);
    A3 = ((I3-I1)*(2*I2*T-H^2) + (I3-I2)*(2*I1*T-H^2))/(I1*I2*I3);

    B1 = 2*(I1-I2)*(I1-I3)/(I2*I3);
    B2 = 2*(I2-I1)*(I2-I3)/(I1*I3);
    B3 = 2*(I3-I1)*(I3-I2)/(I1*I2);

    % ---------------------
    % Computing K from duffing eqn
    % ---------------------
    dw1 = -(I3-I2)*ws_out(2,kk)*ws_out(3,kk)/I1;
    dw2 = -(I1-I3)*ws_out(1,kk)*ws_out(3,kk)/I2;
    dw3 = -(I2-I1)*ws_out(1,kk)*ws_out(2,kk)/I3;

    K1 = dw1^2 + A1*ws_out(1,kk)^2 + (B1/2)*ws_out(1,kk)^4;
    K2 = dw2^2 + A2*ws_out(2,kk)^2 + (B2/2)*ws_out(2,kk)^4;
    K3 = dw3^2 + A3*ws_out(3,kk)^2 + (B3/2)*ws_out(3,kk)^4;

    Ks{ci,1}(:,kk) = [K1;K2;K3];

    % ---------------------
    % Verifying constant K (slide 43)
    % ---------------------
    K1 = (2*I2*T-H^2)*(H^2-2*I3*T)/(I2*I3*I1^2);
    K2 = (2*I3*T-H^2)*(H^2-2*I1*T)/(I1*I3*I2^2);
    K3 = (2*I1*T-H^2)*(H^2-2*I2*T)/(I1*I2*I3^2);

    KsVf{ci,1}(:,kk) = [K1;K2;K3];
end
figure
plot(time_out,Ks{ci,1}(1,:) - KsVf{ci,1}(1,:))
PlotBoi2('Time, sec','\DeltaK',14);
title([titles{ci},' - Difference between computed K values'])

figure; 
subplot(3,1,1)
plot(time_out, percentchange(KsVf{ci,1}(1,:)))
PlotBoi2('','K1 - %Change',12)
title([titles{ci}])
subplot(3,1,2)
plot(time_out, percentchange(KsVf{ci,1}(2,:)))
PlotBoi2('','K2 - %Change',12)
subplot(3,1,3)
plot(time_out, percentchange(KsVf{ci,1}(3,:)))
PlotBoi2('Time, sec','K3 - %Change',12)

end % cases
% distFig('Screen','Secondary') % 'Main' or 'Secondary'
end % run_P4

% ------------------------------------------------------------------------
%%% Problem 5, 7
% ------------------------------------------------------------------------
if run_P57 == 1
% ----------------------
%%% Problem 5
% ----------------------
syms s1 s2 s3 sd1 sd2 sd3 sdd1 sdd2 sdd3 sSquared n Rc G Me I11 I22 I33 real
s = [s1;s2;s3];
sd = [sd1; sd2; sd3];
st = [0 -s3 s2; s3 0 -s1; -s2 s1 0];

BO = eye(3) + (8*st*st - 4*(1-sSquared)*st)./(1+sSquared)^2

Rc_B = BO * [0;0;Rc];

Rc_B(2)*Rc_B(3)/Rc^2;
Rc_B(1)*Rc_B(3)/Rc^2;
Rc_B(1)*Rc_B(2)/Rc^2;

% ----------------------
%%% Problem 7
% ----------------------
I = blkdiag(I11,I22,I33);
wON_O = [0; n; 0];
% wON_N = -wON_O;

wON_B = BO * wON_O
wBO_B = (4/(1+sSquared)^2)*((1-sSquared)*eye(3) - 2*st + 2*(s*s'))*sd

t = simplify(wBO_B.*((sSquared + 1)^2)./4)
% wBN_B = wON_B + wBO_B

wBN_B_lin = [4*sd1 + 4*n*s3; 4*sd2 + n; 4*sd3 - 4*n*s1];
wt_lin = [0, -wBN_B_lin(3), wBN_B_lin(2); wBN_B_lin(3), 0, -wBN_B_lin(1); -wBN_B_lin(2), wBN_B_lin(1), 0];
wdBN_B_lin = [4*sdd1 + 4*n*sd3; 4*sdd2; 4*sdd3 - 4*n*sd1];

LG_B_lin = 3*(n^2).* [4*s1*(I33-I22); -4*s2*(I11-I33); 0]

I*wdBN_B_lin
-wt_lin*I*wBN_B_lin + LG_B_lin

end % run_P57



% ========================================================================
%%% Functions
% ========================================================================
function [dX] = int_mrpW(t,X,I)
%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack Inertia
I1 = I(1,1);
I2 = I(2,2);
I3 = I(3,3);

%%% Unpack the state vector
s1 = X(1); s2 = X(2); s3 = X(3);
sigs = [s1;s2;s3];
s = norm(sigs);

if s > 1
    sigs = -sigs./(s^2);
    s1 = sigs(1); s2 = sigs(2); s3 = sigs(3);
    s = norm(sigs);
end

w1 = X(4); w2 = X(5); w3 = X(6);

%%% Equations of Motion
A = [1-s^2+2*s1^2, 2*(s1*s2-s3), 2*(s1*s3+s2);...
     2*(s2*s1+s3), 1-s^2+2*s2^2, 2*(s2*s3-s1);...
     2*(s3*s1-s2), 2*(s3*s2+s1), 1-s^2+2*s3^2];
ds = .25*A*[w1;w2;w3];

dw1 = -(I3-I2)*w2*w3/I1;
dw2 = -(I1-I3)*w1*w3/I2;
dw3 = -(I2-I1)*w1*w2/I3;
        
%%% Output the derivative of the state
dX(1:3) = [ds(1); ds(2); ds(3)];
dX(4:6) = [dw1; dw2; dw3];
end % function


function [value, isterminal, direction] = int_mrpW_event(t,X,I)

s1 = X(1); s2 = X(2); s3 = X(3);
sigs = [s1;s2;s3];
s = norm(sigs);

value = s - 1; % When the surface is impacted
isterminal = 1; % stops the integration
direction = 1; % negative direction only
end





