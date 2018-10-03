clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

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


H = A-B*K
L = eye(5);
R_wi = blkdiag(0.1, 0.2, 0.1, 0.3, 0.2);

% w = rand(5,1);

y_tf = -0.0035;
x_tf = [0.0094, 0.0069, -0.0035, 0.0388, -0.0808]';
tf   = 4;
p_tf = [0 0 1 0 0]';

time0 = tf:-.0001:0;

%%% Choosing ode45 tolerance
tol = 1e-10;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Integrating
[t, p] = ode45(@Int_GuidanceHW2P2, time0, p_tf, options,H);

t = flipud(t);
p = flipud(p);

plot(t,p,'b','linewidth',1.5)
PlotBoi2('Time, sec','p (adjoint components)',16,'LaTex')

E = zeros(size(p,1),1);
s2s_almost = zeros(size(p));
for kk = 1:size(p,1)
    s2s_almost(kk,1) = (p(kk,1)*L(1,1)) * R_wi(1,1) * (p(kk,1)*L(1,1))';
    s2s_almost(kk,2) = (p(kk,2)*L(2,2)) * R_wi(2,2) * (p(kk,2)*L(2,2))';
    s2s_almost(kk,3) = (p(kk,3)*L(3,3)) * R_wi(3,3) * (p(kk,3)*L(3,3))';
    s2s_almost(kk,4) = (p(kk,4)*L(4,4)) * R_wi(4,4) * (p(kk,4)*L(4,4))';
    s2s_almost(kk,5) = (p(kk,5)*L(5,5)) * R_wi(5,5) * (p(kk,5)*L(5,5))';
    
    
    
%     E(kk) = p(kk,:) * L * R_wi * L' * p(kk,:)';
end

% s2_y = trapz(t,E);

s2s = trapz(t,s2s_almost)















function [dP] = Int_GuidanceHW2P2(t,P,A)
%%% Preallocate state output
dP = -A' * P;

end
































