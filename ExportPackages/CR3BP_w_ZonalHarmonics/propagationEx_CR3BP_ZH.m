% ========================================================================
%%% Description
% ========================================================================
% Example trajectory integrated in the CR3BP with J2, and J4 of the
% primary and J2 of the secondary. Example takes place in Jupiter-Europa
% system

% Created: 12/11/20
% Author : Luke Bury, luke.bury@colorado.edu
%
% ========================================================================
%%% References
% ========================================================================
% 1) Luke Bury and Jay McMahon. The Effect of Zonal Harmonics on Dynamical
% Structures in the Circular Restricted Three-Body Problem Near the
% Secondary Body. Celestrial Mechanics and Dynamical Astronomy, 132(45),
% 2020.
%
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
ticWhole = tic;

% ========================================================================
%%% Run Switches
% ========================================================================
%%% Choose which method to use for calculating normalized mean motion [1]
use_ephemeris_method = true;
use_theory_method    = false;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Set initial conditions
% -------------------------------------------------
%%% Normalized state [6x1]
X0 = [1.009;
      0;
      0;
      0;
      0.06364003647575;
      0];

%%% Normalized final time
Tf = 3.61769928;

% -------------------------------------------------
%%% Necessary parameters
% -------------------------------------------------
%%% Mass ratio of Jupiter-Europa system
m1 = 1.898190000000000e+27;
m2 = 4.799000000000000e+22;
mu = m2 / (m1 + m2); 

%%% Semi-major axis
a = 671100; % km

%%% Body radii (dimensionless)
R1n = 69911  / a; 
R2n = 1560.8 / a; 

%%% Set desired zonal harmonics terms 
%%% *note: the 'p' and 's' denote terms for the 'primary' or 'secondary'
%%% body, so here we have J2, and J4 set for the primary, and J2 set
%%% for the secondary
J2p = 0.014695624770700;
J4p = -5.913138887463000e-04;
J2s = 4.333968885716003e-04;

%%% Measured mean motion value of system
n_real = 2.047827248637295e-05; % rad/s

%%% Gravitational constant
G = 6.6726e-20; % km^3 * kg^-1 * s^-2

% -------------------------------------------------
%%% Calculate normalized mean motion
% -------------------------------------------------
if use_ephemeris_method
    n = n_real*sqrt((a^3)/(G*(m1+m2)));
elseif use_theory_method
    n = ((15*J4p*R1n^4*(mu - 1))/8 - (3*J2p*R1n^2*(mu - 1))/2 + (3*J2s*R2n^2*mu)/2 + 1)^(1/2);
end

if use_ephemeris_method && use_theory_method
    warning('Must only choose one method for normalized mean motion')
    return
end
% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol     = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Integrate and Plot
% ========================================================================
% -------------------------------------------------
%%% Integrate trajectory
% -------------------------------------------------
% [~, X] = ode113(@Int_CR3Bn_ZH, [0, Tf], X0, options, mu, n, R1n, R2n, J2p, J4p, J2s);
[~, X] = ode113(@Int_CR3Bn_ZH, linspace(0, Tf, 10000), X0, options, mu, n, R1n, R2n, J2p, J4p, J2s);

% -------------------------------------------------
%%% Plot trajectory
% -------------------------------------------------
figure; hold all
% axis equal
xlabel('$x_n$','FontName','Times New Roman','Fontsize',20,'Interpreter','LaTex')
ylabel('$y_n$','FontName','Times New Roman','Fontsize',20,'Interpreter','LaTex')
zlabel('$z_n$','FontName','Times New Roman','Fontsize',20,'Interpreter','LaTex')
grid on
set(gcf,'color','white')

plot3(X(:,1), X(:,2), X(:,3),'b','linewidth',1.5)

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\n(Elapsed time: %1.4f seconds)\n',tocWhole)

% ========================================================================
%%% Integration function
% ========================================================================
function [dX] = Int_CR3Bn_ZH(t, X, mu, n, R1n, R2n, J2p, J4p, J2s)
%%% For numerical integration in the normalized CR3BP with capabilities for
%%% zonal harmonics terms from either body (p-primary, s-secondary)
%%% Inputs:
%          t   - normalized time vector
%          X   - initial state [6x1]
%          mu  - mass ratio of CR3BP system
%          n   - normalized mean motion of system
%          R1n - normalized radius of primary body
%          R2n - normalized radius of secondary body
%          J2p - J2 of primary body
%          J4p - J4 of primary body
%          J2s - J2 of secondary body
% =======================================================================
%%% Preallocate output as column vector
dX = zeros(6,1);

%%% Unpack state vector
x = X(1); y = X(2); z = X(3); 
xd = X(4); yd = X(5); zd = X(6);

%%% Calculate acceleration terms
xdd = 2*n*yd + n^2*x + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2)) - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2)) - (J2p*R1n^2*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + (J2s*R2n^2*mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) + (J4p*R1n^4*(6*(2*mu + 2*x)*((mu + x)^2 + y^2 + z^2) - 30*z^2*(2*mu + 2*x))*(mu - 1))/(8*((mu + x)^2 + y^2 + z^2)^(9/2)) - (9*J4p*R1n^4*(2*mu + 2*x)*(mu - 1)*(3*((mu + x)^2 + y^2 + z^2)^2 + 35*z^4 - 30*z^2*((mu + x)^2 + y^2 + z^2)))/(16*((mu + x)^2 + y^2 + z^2)^(11/2)) + (5*J2p*R1n^2*(2*mu + 2*x)*(mu - 1)*((mu + x)^2 + y^2 - 2*z^2))/(4*((mu + x)^2 + y^2 + z^2)^(7/2)) - (5*J2s*R2n^2*mu*(2*mu + 2*x - 2)*((mu + x - 1)^2 + y^2 - 2*z^2))/(4*((mu + x - 1)^2 + y^2 + z^2)^(7/2));
ydd = n^2*y - 2*n*xd - (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) + (J2s*R2n^2*mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (J2p*R1n^2*y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (J4p*R1n^4*(12*y*((mu + x)^2 + y^2 + z^2) - 60*y*z^2)*(mu - 1))/(8*((mu + x)^2 + y^2 + z^2)^(9/2)) - (5*J2s*R2n^2*mu*y*((mu + x - 1)^2 + y^2 - 2*z^2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(7/2)) - (9*J4p*R1n^4*y*(mu - 1)*(3*((mu + x)^2 + y^2 + z^2)^2 + 35*z^4 - 30*z^2*((mu + x)^2 + y^2 + z^2)))/(8*((mu + x)^2 + y^2 + z^2)^(11/2)) + (5*J2p*R1n^2*y*(mu - 1)*((mu + x)^2 + y^2 - 2*z^2))/(2*((mu + x)^2 + y^2 + z^2)^(7/2));
zdd = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (2*J2s*R2n^2*mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (J4p*R1n^4*(48*z*((mu + x)^2 + y^2 + z^2) - 80*z^3)*(mu - 1))/(8*((mu + x)^2 + y^2 + z^2)^(9/2)) + (2*J2p*R1n^2*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) - (5*J2s*R2n^2*mu*z*((mu + x - 1)^2 + y^2 - 2*z^2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(7/2)) - (9*J4p*R1n^4*z*(mu - 1)*(3*((mu + x)^2 + y^2 + z^2)^2 + 35*z^4 - 30*z^2*((mu + x)^2 + y^2 + z^2)))/(8*((mu + x)^2 + y^2 + z^2)^(11/2)) + (5*J2p*R1n^2*z*(mu - 1)*((mu + x)^2 + y^2 - 2*z^2))/(2*((mu + x)^2 + y^2 + z^2)^(7/2));

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; 
dX(4:6) = [xdd; ydd; zdd]; 
end
















