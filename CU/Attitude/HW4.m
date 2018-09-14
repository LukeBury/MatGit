clear
clc
% close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))


run_P3 = 0;
run_P4 = 1;
run_P5 = 0;

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
% *** note on structure of sigs and w0cases *** 
% {case 1 nominal  ,  case 1 nudged;
%  case 2 nominal  ,  case 2 nudged;
%  case 3 nominal  ,  case 3 nudged}
sigs = {zeros(3,length(time)+1),zeros(3,length(time)+1), 'w1 case - nominal','w1 case - nudged';...
        zeros(3,length(time)+1),zeros(3,length(time)+1), 'w2 case - nominal','w2 case - nudged';...
        zeros(3,length(time)+1),zeros(3,length(time)+1), 'w3 case - nominal','w3 case - nudged'};
    
ws = {zeros(3,length(time)+1),zeros(3,length(time)+1), 'w1 case - nominal','w1 case - nudged';...
        zeros(3,length(time)+1),zeros(3,length(time)+1), 'w2 case - nominal','w2 case - nudged';...
        zeros(3,length(time)+1),zeros(3,length(time)+1), 'w3 case - nominal','w3 case - nudged'};

% % Switch times and indexes
% switchTime  = {[],[];[],[];[],[]};
% switchIndex = {[],[];[],[];[],[]};

%%% creating cases
w0cases = {[w,0,0], [w,0,w_n];...
         [0,w,0], [0,w,w_n];...
         [0,0,w], [w_n,0,w]};

%%% Running
for ci = 1:size(w0cases,1)
    w_nom = w0cases{ci,1}'; % nominal
    w_nud = w0cases{ci,2}'; % nudged
    
    w1_nom_new = w_nom(1); w2_nom_new = w_nom(2); w3_nom_new = w_nom(3); % rad/s
    w1_nud_new = w_nud(1); w2_nud_new = w_nud(2); w3_nud_new = w_nud(3); % rad/s
    
    ws{ci,1}(:,1) = w_nom';
    ws{ci,2}(:,1) = w_nud';
    
    for ti = 1:length(time)
        % ---------------------------------------------------------
        % Nominal case
        % ---------------------------------------------------------
        % -------------------
        % Integrating w
        % -------------------
        %%% Assigning friendly names
        w1_nom_old = w1_nom_new; w2_nom_old = w2_nom_new; w3_nom_old = w3_nom_new; % rad/s
        
        %%% Calculating wdots
        w1dot = -(I3-I2)*w2_nom_old*w3_nom_old/I1;
        w2dot = -(I1-I3)*w1_nom_old*w3_nom_old/I2;
        w3dot = -(I2-I1)*w1_nom_old*w2_nom_old/I3;
        
        %%% Updating ws
        w1_nom_new = w1_nom_old + w1dot*dt;
        w2_nom_new = w2_nom_old + w2dot*dt;
        w3_nom_new = w3_nom_old + w3dot*dt;
        
        w_nom_old = [w1_nom_old; w2_nom_old; w3_nom_old]; % rad/s
        w_nom_new = [w1_nom_new; w2_nom_new; w3_nom_new]; % rad/s
        
        ws{ci,1}(:,ti+1) = w_nom_new;
        % -------------------
        % Integrating attitude
        % -------------------
        %%% Assigning variables for swiftness!
        s1_nom = sigs{ci,1}(1,ti);
        s2_nom = sigs{ci,1}(2,ti);
        s3_nom = sigs{ci,1}(3,ti); % so swift!
        s_nom  = norm(sigs{ci,1}(:,ti));
        
        %%% Creating A matrix and updating MRP values
        A_nom = [1-s_nom^2+2*s1_nom^2, 2*(s1_nom*s2_nom-s3_nom), 2*(s1_nom*s3_nom+s2_nom);...
                   2*(s2_nom*s1_nom+s3_nom), 1-s_nom^2+2*s2_nom^2, 2*(s2_nom*s3_nom-s1_nom);...
                   2*(s3_nom*s1_nom-s2_nom), 2*(s3_nom*s2_nom+s1_nom), 1-s_nom^2+s_nom*s3_nom^2];
        rate_nom = .25*A_nom*w_nom_old;
        
        %%% Updating Sigma & Shadow values
        sigs{ci,1}(:,ti+1) = sigs{ci,1}(:,ti) + rate_nom.*dt;
        
        %%% Check if need to switch to shadow set
        if norm(sigs{ci,1}(:,ti+1)) > 1
%             switchIndex{ci,1} = [switchIndex{ci,1}, ti];
%             switchTime{ci,1}  = [switchTime{ci,1},time(ti)];
            sigs{ci,1}(:,ti+1) = -sigs{ci,1}(:,ti+1)./(norm(sigs{ci,1}(:,ti+1))^2);
        end
        
        % ---------------------------------------------------------
        % Nudged case
        % ---------------------------------------------------------
        % -------------------
        % Integrating w
        % -------------------
        %%% Assigning friendly names
        w1_nud_old = w1_nud_new; w2_nud_old = w2_nud_new; w3_nud_old = w3_nud_new; % rad/s
        
        %%% Calculating wdots
        w1dot = -(I3-I2)*w2_nud_old*w3_nud_old/I1;
        w2dot = -(I1-I3)*w1_nud_old*w3_nud_old/I2;
        w3dot = -(I2-I1)*w1_nud_old*w2_nud_old/I3;
        
        %%% Updating ws
        w1_nud_new = w1_nud_old + w1dot*dt;
        w2_nud_new = w2_nud_old + w2dot*dt;
        w3_nud_new = w3_nud_old + w3dot*dt;
        
        w_nud_new = [w1_nud_new; w2_nud_new; w3_nud_new]; % rad/s
        
        ws{ci,2}(:,ti+1) = w_nud_new;
        % -------------------
        % Integrating attitude
        % -------------------
        %%% Assigning variables for swiftness!
        s1_nud = sigs{ci,2}(1,ti); 
        s2_nud = sigs{ci,2}(2,ti);
        s3_nud = sigs{ci,2}(3,ti); % so swift!
        s_nud  = norm(sigs{ci,2}(:,ti));
        
        %%% Creating A matrix and updating MRP values
        A_nud = [1-s_nud^2+2*s1_nud^2, 2*(s1_nud*s2_nud-s3_nud), 2*(s1_nud*s3_nud+s2_nud);...
                   2*(s2_nud*s1_nud+s3_nud), 1-s_nud^2+2*s2_nud^2, 2*(s2_nud*s3_nud-s1_nud);...
                   2*(s3_nud*s1_nud-s2_nud), 2*(s3_nud*s2_nud+s1_nud), 1-s_nud^2+s_nud*s3_nud^2];
        rate_nud = .25*A_nud*w_nud_new;
        
        %%% Updating Sigma & Shadow values
        sigs{ci,2}(:,ti+1) = sigs{ci,2}(:,ti) + rate_nud.*dt;
        
        %%% Check if need to switch to shadow set
        if norm(sigs{ci,2}(:,ti+1)) > 1
%             switchIndex{ci,2} = [switchIndex{ci,2}, ti];
%             switchTime{ci,2}  = [switchTime{ci,2},time(ti)];
            sigs{ci,2}(:,ti+1) = -sigs{ci,2}(:,ti+1)./(norm(sigs{ci,2}(:,ti+1))^2);
        end
    end
    
    figure; hold all
    plot(time,ws{ci,1}(1,1:end-1),'.','linewidth',3)
    plot(time,ws{ci,1}(2,1:end-1),'.','linewidth',3)
    plot(time,ws{ci,1}(3,1:end-1),'.','linewidth',3)
    PlotBoi2('Time, sec','Angular Velocity, rad/s',14)
    legend('\omega1', '\omega2', '\omega3')
    title(ws{ci,3})
    
    figure; hold all
    plot(time,sigs{ci,1}(1,1:end-1),'.','linewidth',3)
    plot(time,sigs{ci,1}(2,1:end-1),'.','linewidth',3)
    plot(time,sigs{ci,1}(3,1:end-1),'.','linewidth',3)
    PlotBoi2('Time, sec','MRP Value',14)
    legend('\sigma1', '\sigma2', '\sigma3')
    title(sigs{ci,3})
    
     figure; hold all
    plot(time,ws{ci,2}(1,1:end-1),'.','linewidth',3)
    plot(time,ws{ci,2}(2,1:end-1),'.','linewidth',3)
    plot(time,ws{ci,2}(3,1:end-1),'.','linewidth',3)
    PlotBoi2('Time, sec','Angular Velocity, rad/s',14)
    legend('\omega1', '\omega2', '\omega3')
    title(ws{ci,4})
    
    figure; hold all
    plot(time,sigs{ci,2}(1,1:end-1),'.','linewidth',3)
    plot(time,sigs{ci,2}(2,1:end-1),'.','linewidth',3)
    plot(time,sigs{ci,2}(3,1:end-1),'.','linewidth',3)
    PlotBoi2('Time, sec','MRP Value',14)
    legend('\sigma1', '\sigma2', '\sigma3')
    title(sigs{ci,4})

end


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
sigs = {zeros(3,length(time)+1),'w1 case - nudged';...
        zeros(3,length(time)+1),'w2 case - nudged';...
        zeros(3,length(time)+1),'w3 case - nudged'};
    
ws = {zeros(3,length(time)+1),'w1 case - nudged';...
        zeros(3,length(time)+1),'w2 case - nudged';...
        zeros(3,length(time)+1),'w3 case - nudged'};

Ks = {zeros(3,length(time)+1),'w1 case - nudged';...
      zeros(3,length(time)+1),'w2 case - nudged';...
      zeros(3,length(time)+1),'w3 case - nudged'};
    
KsVf = {zeros(3,length(time)+1),'w1 case - nudged';...
        zeros(3,length(time)+1),'w2 case - nudged';...
        zeros(3,length(time)+1),'w3 case - nudged'};

%%% creating cases
w0cases = {[w,0,w_n];...
         [0,w,w_n];...
         [w_n,0,w]};
     
%%% Running
for ci = 1:1
    w_nud = w0cases{ci,1}'; % nudged
    
    w1_nud_new = w_nud(1); w2_nud_new = w_nud(2); w3_nud_new = w_nud(3); % rad/s
    
    ws{ci,1}(:,1) = w_nud';
    
    for ti = 1:length(time)
        % ---------------------------------------------------------
        % Integrating w
        % ---------------------------------------------------------
        %%% Assigning friendly names
        w1_nud_old = w1_nud_new; w2_nud_old = w2_nud_new; w3_nud_old = w3_nud_new; % rad/s
        
        %%% Calculating wdots
        w1dot_old = -(I3-I2)*w2_nud_old*w3_nud_old/I1;
        w2dot_old = -(I1-I3)*w1_nud_old*w3_nud_old/I2;
        w3dot_old = -(I2-I1)*w1_nud_old*w2_nud_old/I3;
        
        %%% Updating ws
        w1_nud_new = w1_nud_old + w1dot_old*dt;
        w2_nud_new = w2_nud_old + w2dot_old*dt;
        w3_nud_new = w3_nud_old + w3dot_old*dt;
        
        %%% Calculating new wdots for later use
        w1dot_new = -(I3-I2)*w2_nud_new*w3_nud_new/I1;
        w2dot_new = -(I1-I3)*w1_nud_new*w3_nud_new/I2;
        w3dot_new = -(I2-I1)*w1_nud_new*w2_nud_new/I3;
        
        w_nud_old = [w1_nud_old; w2_nud_old; w3_nud_old]; % rad/s
        w_nud_new = [w1_nud_new; w2_nud_new; w3_nud_new]; % rad/s
        
        ws{ci,1}(:,ti+1) = w_nud_new;
        % ---------------------------------------------------------
        % Integrating attitude
        % ---------------------------------------------------------
        %%% Assigning variables for swiftness!
        s1_nud = sigs{ci,1}(1,ti); 
        s2_nud = sigs{ci,1}(2,ti);
        s3_nud = sigs{ci,1}(3,ti); % so swift!
        s_nud  = norm(sigs{ci,1}(:,ti));
        
        %%% Creating A matrix and updating MRP values
        A_nud = [1-s_nud^2+2*s1_nud^2, 2*(s1_nud*s2_nud-s3_nud), 2*(s1_nud*s3_nud+s2_nud);...
                   2*(s2_nud*s1_nud+s3_nud), 1-s_nud^2+2*s2_nud^2, 2*(s2_nud*s3_nud-s1_nud);...
                   2*(s3_nud*s1_nud-s2_nud), 2*(s3_nud*s2_nud+s1_nud), 1-s_nud^2+s_nud*s3_nud^2];
%         rate_nud = .25*A_nud*[w1_nud_old;w2_nud_old;w3_nud_old];
        rate_nud = .25*A_nud*w_nud_old;
        
        %%% Updating Sigma & Shadow values
        sigs{ci,1}(:,ti+1) = sigs{ci,1}(:,ti) + rate_nud.*dt;
        
        %%% Check if need to switch to shadow set
        if norm(sigs{ci,1}(:,ti+1)) > 1
            sigs{ci,1}(:,ti+1) = -sigs{ci,1}(:,ti+1)./(norm(sigs{ci,1}(:,ti+1))^2);
        end
        
        % ---------------------------------------------------------
        % Duffing Equation work
        % ---------------------------------------------------------
        % ---------------------
        % Computing T and H
        % ---------------------
        T = (1/2) * ws{ci,1}(:,ti+1)' * I * ws{ci,1}(:,ti+1);
        H = I * ws{ci,1}(:,ti+1);
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
        K1 = w1dot_new^2 + A1*ws{ci,1}(1,ti+1)^2 + (B1/2)*ws{ci,1}(1,ti+1)^4;
        K2 = w2dot_new^2 + A2*ws{ci,1}(2,ti+1)^2 + (B2/2)*ws{ci,1}(2,ti+1)^4;
        K3 = w3dot_new^2 + A3*ws{ci,1}(3,ti+1)^2 + (B3/2)*ws{ci,1}(3,ti+1)^4;
%         K1 = w1dot^2 + A1*ws{ci,1}(1,ti)^2 + (B1/2)*ws{ci,1}(1,ti)^4;
%         K2 = w2dot^2 + A2*ws{ci,1}(2,ti)^2 + (B2/2)*ws{ci,1}(2,ti)^4;
%         K3 = w3dot^2 + A3*ws{ci,1}(3,ti)^2 + (B3/2)*ws{ci,1}(3,ti)^4;
        
        Ks{ci,1}(:,ti+1) = [K1;K2;K3];
        
        % ---------------------
        % Verifying constant K (slide 43)
        % ---------------------
        K1 = (2*I2*T-H^2)*(H^2-2*I3*T)/(I2*I3*I1^2);
        K2 = (2*I3*T-H^2)*(H^2-2*I1*T)/(I1*I3*I2^2);
        K3 = (2*I1*T-H^2)*(H^2-2*I2*T)/(I1*I2*I3^2);
        
        KsVf{ci,1}(:,ti+1) = [K1;K2;K3];
    end
end
figure
plot(time(2:end),Ks{ci,1}(1,2:end-1) - KsVf{ci,1}(1,2:end-1))

figure; 
subplot(3,1,1)
plot(time(2:end), KsVf{1,1}(1,2:end-1))
title(KsVf{1,2})
subplot(3,1,2)
plot(time(2:end), KsVf{1,1}(2,2:end-1))
subplot(3,1,3)
plot(time(2:end), KsVf{1,1}(3,2:end-1))

figure; 
subplot(3,1,1)
plot(time(2:end), percentchange(KsVf{1,1}(1,2:end-1)))
title(KsVf{1,2})
subplot(3,1,2)
plot(time(2:end), percentchange(KsVf{1,1}(2,2:end-1)))
subplot(3,1,3)
plot(time(2:end), percentchange(KsVf{1,1}(3,2:end-1)))

% figure; hold all
% subplot(3,1,1)
% plot(time(2:end), percentchange(KsVf{2,1}(1,2:end-1)))
% title(KsVf{2,2})
% subplot(3,1,2)
% plot(time(2:end), percentchange(KsVf{2,1}(2,2:end-1)))
% subplot(3,1,3)
% plot(time(2:end), percentchange(KsVf{2,1}(3,2:end-1)))
% 
% figure; hold all
% subplot(3,1,1)
% plot(time(2:end), percentchange(KsVf{3,1}(1,2:end-1)))
% title(KsVf{3,2})
% subplot(3,1,2)
% plot(time(2:end), percentchange(KsVf{3,1}(2,2:end-1)))
% subplot(3,1,3)
% plot(time(2:end), percentchange(KsVf{3,1}(3,2:end-1)))



% =========================================
%%% Choosing ode45 tolerance
tol = 1e-6;

%%% Setting integrator options
options = odeset('Events',@int_mrpW_event,'RelTol',tol,'AbsTol',tol);

%%% IC
X0 = [0;0;0;w;0;w_n];


%%% Propagating the States without J2
[time_o, mrpW] = ode45(@int_mrpW, time, X0, options, I);
c = 0;
while time_o(end) ~= time(end)
    c = c + 1;
    timeIndex = find(round(time,3)==round(time_o(end-1),3));
    mrpW(timeIndex,1:3) = -mrpW(timeIndex,1:3)./(norm(mrpW(timeIndex,1:3))^2);
    [time_n, mrpW_n] = ode45(@int_mrpW, [time(timeIndex):dt:tf], mrpW(timeIndex,:)', options, I);
    mrpW = [mrpW; mrpW_n];
    time_o = [time_o;time_n];
end
[~,t_idx] = unique(time_o);
time_o = time_o(t_idx);
mrpW = mrpW(t_idx,:);
% =========================================


end % run_P4

% ------------------------------------------------------------------------
%%% Problem 5
% ------------------------------------------------------------------------
if run_P5 == 1
syms s1 s2 s3 sSquared n Rc G Me real
s = [s1;s2;s3];
st = [0 -s3 s2; s3 0 -s1; -s2 s1 0];

BO = eye(3) + (8*st*st - 4*(1-sSquared)*st)./(1+sSquared)^2

Rc_B = BO * [0;0;Rc]

Rc_B(2)*Rc_B(3)/Rc^2
Rc_B(1)*Rc_B(3)/Rc^2
Rc_B(1)*Rc_B(2)/Rc^2

% 
% wON_O = [0; n; 0];
% wON_N = -wON_O;
% 
% wON_B = BO * wON_O;
% 
% wBN = w


end % run_P4



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

%%% Unpack the barycentric state vector
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
     2*(s3*s1-s2), 2*(s3*s2+s1), 1-s^2+s*s3^2];
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





