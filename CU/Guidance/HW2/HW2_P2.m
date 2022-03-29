clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Work
% ========================================================================
A = [5 0 3 0 1;...
    3 0 0 -2 0;...
    0 -2 4 1 0;...
    1 3 -4 1 3;...
    0 2 2 0 -1];

B = [0 1;...
    0 2;...
    0 0;...
    1 3;...
    1 1];

K = [-2.5605   -2.5308    6.9961    3.7510    1.1419;...
   21.8698   -6.9988   26.2313    2.4302    2.3425];


H = A-B*K;
L = eye(5);
R_wi = blkdiag(0.1, 0.2, 0.1, 0.3, 0.2);

% w = rand(5,1);

y_tf = 0;
x_tf = [0, 0, 0, 0, 0]';
tf   = 50; % 10
p_tf = [0 0 1 0 0]';

time0_bkwd = tf:-.001:0;

%%% Choosing ode45 tolerance
tol = 1e-6;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Integrating
[t, p_bkwd] = ode45(@Int_GuidanceHW2P2, time0_bkwd, p_tf, options,H);

t = flipud(t);
p_fwd = flipud(p_bkwd);

figure
subplot(2,1,1)
plot(t,p_fwd,'b','linewidth',1.5)
PlotBoi2('Time, sec','p (adjoint components)',16,'LaTex')
subplot(2,1,2)
indx = ceil(0.95 * size(p_fwd,1));
plot(t(indx:end),p_fwd(indx:end,:),'b','linewidth',1.5)
PlotBoi2('Time, sec','p (adjoint components)',16,'LaTex')

E = zeros(size(p_fwd,1),1);
s2s_almost = zeros(size(p_fwd));
for kk = 1:size(p_fwd,1)
    s2s_almost(kk,1) = (p_fwd(kk,1)*L(1,1)) * R_wi(1,1) * (p_fwd(kk,1)*L(1,1))';
    s2s_almost(kk,2) = (p_fwd(kk,2)*L(2,2)) * R_wi(2,2) * (p_fwd(kk,2)*L(2,2))';
    s2s_almost(kk,3) = (p_fwd(kk,3)*L(3,3)) * R_wi(3,3) * (p_fwd(kk,3)*L(3,3))';
    s2s_almost(kk,4) = (p_fwd(kk,4)*L(4,4)) * R_wi(4,4) * (p_fwd(kk,4)*L(4,4))';
    s2s_almost(kk,5) = (p_fwd(kk,5)*L(5,5)) * R_wi(5,5) * (p_fwd(kk,5)*L(5,5))';
    
    
    
%     E(kk) = p(kk,:) * L * R_wi * L' * p(kk,:)';
end

% s2_y = trapz(t,E);

s2s = trapz(t,s2s_almost)
s2_total = sum(s2s)


% -------------------------------------------------
% Part B
% -------------------------------------------------
time0_b_backwards = time0_bkwd;
time0_b_forewards = fliplr(time0_b_backwards);

x0 = [1 1 1 1 1]';


gamma = 1;
a = 1;
b = 1;

R = [a 0; 0 b];

%  Q is so that x' Q x = gamma*x3^2
Q = [0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 gamma 0 0;...
    0 0 0 0 0;
    0 0 0 0 0];

%%% What I call 'S' is called K on slide 34 of "4 Optimal Control" slide
%%% package. Slide 32 has the info for why S_tf is zero. In slide notation,
%%% K_tf (my S_tf) is equal to S, which is shown in the cost function, but
%%% our cost function has that S = 0
S_tf = zeros(25,1);

%%% Integrating
[t_bkwds, S_out_bkwds] = ode45(@Int_GuidanceHW2P2_S, time0_b_backwards, S_tf, options,B,R,Q,A);
% S_out = flipud(S_out);

S = zeros(5,5,size(S_out_bkwds,1));
for kk = 1:size(S_out_bkwds,1)
    S(:,:,kk) = reshape(S_out_bkwds(kk,:),5,5);
end


time_fwds = flipud(t_bkwds);
S = flipud(S);
S_out_fwds = flipud(S_out_bkwds);


[t_out, X_out] = ode45(@Int_GuidanceHW2P2_State, time0_b_forewards,x0,options,A,B,R,S_out_fwds,time0_b_forewards);

load('GuidanceHw2P2Vars.mat') % y and us
us_hw1 = us;
y_hw1 = y;

y_hw2 = X_out(:,3);

if length(y_hw1) == length(y_hw2)
    figure(10); 
    subplot(2,1,1); hold all
    plot(time0_b_forewards,y_hw2,'linewidth',2)
    plot(time0_b_forewards,y_hw1,'linewidth',2)
    PlotBoi2('Time, sec','Output, $x_3$',16,'LaTex')
    subplot(2,1,2); hold all
    plot(time0_b_forewards,abs(y_hw1-y_hw2),'r','linewidth',2)
    PlotBoi2('Time, sec','Output Difference',16,'LaTex')

% -------------------------------------------------
% Part C (continued)
% -------------------------------------------------
%%% Determing controls from this problem
us_hw2b = zeros(2,size(X_out,1));
for kk = 1:size(X_out,1)
    S = reshape(S_out_fwds(kk,:),5,5);
    us_hw2b(:,kk) = -inv(R)*B'*S*X_out(kk,:)';
end

%%% Evaluating cost functions
J1 = zeros(size(us_hw1,2),1);
J2 = zeros(size(us_hw2b,2),1);

for kk = 1:length(J1)
    J1(kk) = gamma * y_hw1(kk)^2 + us_hw1(:,kk)'*[a 0; 0 b]*us_hw1(:,kk);
    J2(kk) = gamma * y_hw2(kk)^2 + us_hw2b(:,kk)'*[a 0; 0 b]*us_hw2b(:,kk);
end

J1 = trapz(time0_b_forewards,J1);
J2 = trapz(time0_b_forewards,J2);
Jdiff = J2-J1

fprintf('Costs are very very similar\n')
else
    warning('Skipping Hw1 comparison')
end

% -------------------------------------------------
% Part D
% -------------------------------------------------
% time0_d_backwards = 30:-.001:0;
% time0_d_forewards = fliplr(time0_d_backwards);
%% Options for reference tracking
gamma = 100;
a = .0001;
b = .0001;
tol = 1e-6;

R = [a 0; 0 b];
%  Q is so that x' Q x = gamma*x3^2
Q = [0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 gamma 0 0;...
    0 0 0 0 0;
    0 0 0 0 0];

%%% Choosing ode45 tolerance
options = odeset('RelTol',tol,'AbsTol',tol);

%%% What I call 'S' is called K on slide 34 of "4 Optimal Control" slide
%%% package. Slide 32 has the info for why S_tf is zero. In slide notation,
%%% K_tf (my S_tf) is equal to S, which is shown in the cost function, but
%%% our cost function has that S = 0
S_tf = zeros(25,1);

%%% Integrating
[t_bkwds, S_out_bkwds] = ode45(@Int_GuidanceHW2P2_S, time0_b_backwards, S_tf, options,B,R,Q,A);
% S_out = flipud(S_out);

S = zeros(5,5,size(S_out_bkwds,1));
for kk = 1:size(S_out_bkwds,1)
    S(:,:,kk) = reshape(S_out_bkwds(kk,:),5,5);
end


time_fwds = flipud(t_bkwds);
S = flipud(S);
S_out_fwds = flipud(S_out_bkwds);



xref = sin(pi*time_fwds./5);
xref = [zeros(length(xref),1), zeros(length(xref),1), xref, zeros(length(xref),1), zeros(length(xref),1)];

x0_d = [0 0 0 0 0]';

K_tf = zeros(5,5);
S = K_tf;
s_tf = -S*xref(1,:)';

[t_out_bkwd, s_bkwd] = ode45(@Int_GuidanceHW2P2_LittleS, time0_b_backwards,s_tf,options,A,B,Q,R,S_out_bkwds,time0_b_backwards);

s_fwd = flipud(s_bkwd);

[t_out, X_out_d] = ode45(@Int_GuidanceHW2P2_State_d, time0_b_forewards,x0_d,options,A,B,R,S_out_fwds,s_fwd,time0_b_forewards);

y_hw3 = X_out_d(:,3);

figure
subplot(2,1,1); hold all
p1 = plot(time0_b_forewards,y_hw3,'r','linewidth',2);
p2 = plot(time0_b_forewards,xref(:,3),'b','linewidth',2);
legend([p1 p2],'Output','Reference')
PlotBoi2('','Output, $x_3$, and reference',16,'LaTex')
title(['$\gamma$',sprintf(' = %1.0f, a = %1.0e, b = %1.0e',gamma,a,b)],'interpreter','LaTeX')
subplot(2,1,2); hold all
plot(time0_b_forewards,y_hw3 - xref(:,3),'linewidth',2)
PlotBoi2('Time, sec','Output Difference',16,'LaTex')

%%% Determing controls from this problem
us_hw2d = zeros(2,size(X_out_d,1));
for kk = 1:size(X_out_d,1)
    S = reshape(S_out_fwds(kk,:),5,5);
    s = s_fwd(kk,:)';
    us_hw2d(:,kk) = -inv(R)*B'*S*X_out_d(kk,:)' - inv(R)*B'*s;
end


%%% Evaluating cost functions
J3 = zeros(size(us_hw2d,2),1);

for kk = 1:length(J1)
    J3(kk) = gamma * y_hw3(kk)^2 + us_hw2d(:,kk)'*[a 0; 0 b]*us_hw2d(:,kk);
end

J3 = trapz(time0_b_forewards,J3)


toc















function [dP] = Int_GuidanceHW2P2(t,P,A)
%%% Preallocate state output
dP = -A' * P;

end




function [dS] = Int_GuidanceHW2P2_S(t,S,B,R,Q,A)
S = reshape(S,5,5);

%%% Preallocate state output
dS = -S*A + S*B*inv(R)*B'*S - Q - A'*S;
dS = reshape(dS,25,1);

end




function [dX] = Int_GuidanceHW2P2_State(t,x,A,B,R,S,fullTime)

S_now = interp1(fullTime,S,t);
S = reshape(S_now,5,5);

u = -inv(R)*B'*S*x;
dX = A*x + B*u;

end




function [ds] = Int_GuidanceHW2P2_LittleS(t,s,A,B,Q,R,S,fullTime)
S_now = interp1(fullTime,S,t);
S = reshape(S_now,5,5);

ref = [0; 0; sin(pi*t/5); 0; 0];

ds = -[A' - S * B * inv(R) *B']*s + Q*ref;

end





function [dX] = Int_GuidanceHW2P2_State_d(t,x,A,B,R,S,s,fullTime)

S_now = interp1(fullTime,S,t);
S = reshape(S_now,5,5);

s_now = interp1(fullTime,s,t);
s_now = s_now';

u = -inv(R)*B'*S*x - inv(R)*B'*s_now;
dX = A*x + B*u;

end



























