clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

% ========================================================================
%%% Run switches
% ========================================================================
run_symbollicWork = 0;
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Symbolic Work
% ========================================================================
if run_symbollicWork == 1
syms x y z R1 R2 J21 J22 u dx dy dz real

r1 = sqrt((-u-x)^2+y^2+z^2);
r2 = sqrt((1-u-x)^2+y^2+z^2);

gamma1 = 3*(1-u)*R1^2*J21*(5*z^2-r1^2)/(2*r1^7);
gamma2 = 3*u*R2^2*J22*(5*z^2-r2^2)/(2*r2^7);
gamma1_z = 3*(1-u)*R1^2*J21*(5*z^2-3*r1^2)/(2*r1^7);
gamma2_z = 3*u*R2^2*J22*(5*z^2-3*r2^2)/(2*r2^7);

%%% Defining Equations of Motion, State Variables, and A Matrix
x1 = -u;
x2 = 1-u;

ddx = 2*dy + x + ((1-u)/r1^3-gamma1)*(x1-x) + (u/r2^3 - gamma2)*(x2-x);
ddy = -2*dx + y*(-(1-u)/r1^3 - u/r2^3 + gamma1 + gamma2 + 1);
ddz = z*(-(1-u)/r1^3 - u/r2^3 + gamma1_z + gamma2_z);

EOM = [dx; dy; dz; ddx; ddy; ddz];
state = [x; y; z; dx; dy; dz];
A = jacobian(EOM, state);

989
return
end % run_symbollic work
% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------------------- 
%%% Choosing Bodies
% -------------------------------------------------
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Mass ratio of system
u = secondary.MR;

% ------------------------------------------------- 
%%% Periodic Orbit options
% -------------------------------------------------
n = 300;
L_Point = 2;
% dS_PO = .0001;
dS_PO = .0001;
dynamic_dS_PO = 0;
PO_plot_skip = 50; %5

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-9;
options = odeset('RelTol',tol,'AbsTol',tol);

% ------------------------------------------------- 
%%% Finding Equilibrium Points
% -------------------------------------------------
% Ls_n = EquilibriumPoints(secondary.MR);
Ls_J2_n = EquilibriumPoints_J2(secondary.MR,primary.J2,0,primary.R/rNorm,secondary.R_n);
% ========================================================================
%%% Evaluating A matrix at equilibrium point
% ========================================================================
% ------------------------------------------------- 
%%% Choosing equilibrium point and setting initial state
% -------------------------------------------------
%%% Coordinates of Lagrange Point to be evaluated
xL = Ls_J2_n(L_Point,1);
yL = Ls_J2_n(L_Point,2);
zL = Ls_J2_n(L_Point,3);

%%% Radii of bodies
R1 = primary.R/rNorm;
R2 = secondary.R_n;
J21 = primary.J2;
989
J22 = 0;

% ------------------------------------------------- 
%%% Evaluating A at equilibrium point
% -------------------------------------------------
A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4,5) = 2;
A(5,4) = -2;
A(4,1) = (u + xL - 1)*((3*u*(2*u + 2*xL - 2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(5/2)) - (3*J22*R2^2*u*(2*u + 2*xL - 2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(7/2)) + (21*J22*R2^2*u*(2*u + 2*xL - 2)*((u + xL - 1)^2 + yL^2 - 4*zL^2))/(4*((u + xL - 1)^2 + yL^2 + zL^2)^(9/2))) + (u - 1)/((u + xL)^2 + yL^2 + zL^2)^(3/2) - (u + xL)*((3*(2*u + 2*xL)*(u - 1))/(2*((u + xL)^2 + yL^2 + zL^2)^(5/2)) - (J21*R1^2*(2*u + 2*xL)*(3*u - 3))/(2*((u + xL)^2 + yL^2 + zL^2)^(7/2)) + (7*J21*R1^2*(2*u + 2*xL)*(3*u - 3)*((u + xL)^2 + yL^2 - 4*zL^2))/(4*((u + xL)^2 + yL^2 + zL^2)^(9/2))) - u/((u + xL - 1)^2 + yL^2 + zL^2)^(3/2) - (3*J22*R2^2*u*((u + xL - 1)^2 + yL^2 - 4*zL^2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(7/2)) + (J21*R1^2*(3*u - 3)*((u + xL)^2 + yL^2 - 4*zL^2))/(2*((u + xL)^2 + yL^2 + zL^2)^(7/2)) + 1;
A(4,2) = (u + xL - 1)*((3*u*yL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*J22*R2^2*u*yL)/((u + xL - 1)^2 + yL^2 + zL^2)^(7/2) + (21*J22*R2^2*u*yL*((u + xL - 1)^2 + yL^2 - 4*zL^2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(9/2))) - (u + xL)*((3*yL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2) - (J21*R1^2*yL*(3*u - 3))/((u + xL)^2 + yL^2 + zL^2)^(7/2) + (7*J21*R1^2*yL*(3*u - 3)*((u + xL)^2 + yL^2 - 4*zL^2))/(2*((u + xL)^2 + yL^2 + zL^2)^(9/2)));
A(4,3) = (u + xL - 1)*((3*u*zL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) + (12*J22*R2^2*u*zL)/((u + xL - 1)^2 + yL^2 + zL^2)^(7/2) + (21*J22*R2^2*u*zL*((u + xL - 1)^2 + yL^2 - 4*zL^2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(9/2))) - (u + xL)*((3*zL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2) + (4*J21*R1^2*zL*(3*u - 3))/((u + xL)^2 + yL^2 + zL^2)^(7/2) + (7*J21*R1^2*zL*(3*u - 3)*((u + xL)^2 + yL^2 - 4*zL^2))/(2*((u + xL)^2 + yL^2 + zL^2)^(9/2)));
A(5,1) = -yL*((3*(2*u + 2*xL)*(u - 1))/(2*((u + xL)^2 + yL^2 + zL^2)^(5/2)) - (3*u*(2*u + 2*xL - 2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(5/2)) - (J21*R1^2*(2*u + 2*xL)*(3*u - 3))/(2*((u + xL)^2 + yL^2 + zL^2)^(7/2)) + (3*J22*R2^2*u*(2*u + 2*xL - 2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(7/2)) - (21*J22*R2^2*u*(2*u + 2*xL - 2)*((u + xL - 1)^2 + yL^2 - 4*zL^2))/(4*((u + xL - 1)^2 + yL^2 + zL^2)^(9/2)) + (7*J21*R1^2*(2*u + 2*xL)*(3*u - 3)*((u + xL)^2 + yL^2 - 4*zL^2))/(4*((u + xL)^2 + yL^2 + zL^2)^(9/2)));
A(5,2) = (u - 1)/((u + xL)^2 + yL^2 + zL^2)^(3/2) + yL*((3*u*yL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*yL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2) + (J21*R1^2*yL*(3*u - 3))/((u + xL)^2 + yL^2 + zL^2)^(7/2) - (3*J22*R2^2*u*yL)/((u + xL - 1)^2 + yL^2 + zL^2)^(7/2) - (7*J21*R1^2*yL*(3*u - 3)*((u + xL)^2 + yL^2 - 4*zL^2))/(2*((u + xL)^2 + yL^2 + zL^2)^(9/2)) + (21*J22*R2^2*u*yL*((u + xL - 1)^2 + yL^2 - 4*zL^2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(9/2))) - u/((u + xL - 1)^2 + yL^2 + zL^2)^(3/2) - (3*J22*R2^2*u*((u + xL - 1)^2 + yL^2 - 4*zL^2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(7/2)) + (J21*R1^2*(3*u - 3)*((u + xL)^2 + yL^2 - 4*zL^2))/(2*((u + xL)^2 + yL^2 + zL^2)^(7/2)) + 1;
A(5,3) = yL*((3*u*zL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*zL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2) - (4*J21*R1^2*zL*(3*u - 3))/((u + xL)^2 + yL^2 + zL^2)^(7/2) + (12*J22*R2^2*u*zL)/((u + xL - 1)^2 + yL^2 + zL^2)^(7/2) - (7*J21*R1^2*zL*(3*u - 3)*((u + xL)^2 + yL^2 - 4*zL^2))/(2*((u + xL)^2 + yL^2 + zL^2)^(9/2)) + (21*J22*R2^2*u*zL*((u + xL - 1)^2 + yL^2 - 4*zL^2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(9/2)));
A(6,1) = -zL*((3*(2*u + 2*xL)*(u - 1))/(2*((u + xL)^2 + yL^2 + zL^2)^(5/2)) - (3*u*(2*u + 2*xL - 2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(5/2)) - (J21*R1^2*(6*u + 6*xL)*(3*u - 3))/(2*((u + xL)^2 + yL^2 + zL^2)^(7/2)) + (3*J22*R2^2*u*(6*u + 6*xL - 6))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(7/2)) + (7*J21*R1^2*(2*u + 2*xL)*(3*u - 3)*(3*(u + xL)^2 + 3*yL^2 - 2*zL^2))/(4*((u + xL)^2 + yL^2 + zL^2)^(9/2)) - (21*J22*R2^2*u*(2*u + 2*xL - 2)*(3*(u + xL - 1)^2 + 3*yL^2 - 2*zL^2))/(4*((u + xL - 1)^2 + yL^2 + zL^2)^(9/2)));
A(6,2) = zL*((3*u*yL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*yL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2) + (3*J21*R1^2*yL*(3*u - 3))/((u + xL)^2 + yL^2 + zL^2)^(7/2) - (9*J22*R2^2*u*yL)/((u + xL - 1)^2 + yL^2 + zL^2)^(7/2) + (21*J22*R2^2*u*yL*(3*(u + xL - 1)^2 + 3*yL^2 - 2*zL^2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(9/2)) - (7*J21*R1^2*yL*(3*u - 3)*(3*(u + xL)^2 + 3*yL^2 - 2*zL^2))/(2*((u + xL)^2 + yL^2 + zL^2)^(9/2)));
A(6,3) = zL*((3*u*zL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*zL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2) - (2*J21*R1^2*zL*(3*u - 3))/((u + xL)^2 + yL^2 + zL^2)^(7/2) + (6*J22*R2^2*u*zL)/((u + xL - 1)^2 + yL^2 + zL^2)^(7/2) - (7*J21*R1^2*zL*(3*u - 3)*(3*(u + xL)^2 + 3*yL^2 - 2*zL^2))/(2*((u + xL)^2 + yL^2 + zL^2)^(9/2)) + (21*J22*R2^2*u*zL*(3*(u + xL - 1)^2 + 3*yL^2 - 2*zL^2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(9/2))) + (u - 1)/((u + xL)^2 + yL^2 + zL^2)^(3/2) - u/((u + xL - 1)^2 + yL^2 + zL^2)^(3/2) + (J21*R1^2*(3*u - 3)*(3*(u + xL)^2 + 3*yL^2 - 2*zL^2))/(2*((u + xL)^2 + yL^2 + zL^2)^(7/2)) - (3*J22*R2^2*u*(3*(u + xL - 1)^2 + 3*yL^2 - 2*zL^2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(7/2));

% ------------------------------------------------- 
%%% Eigenvectors and Eigenvalues
% -------------------------------------------------
%%% Getting eigenvectors and eigenvalues
[eVecs, eVals] = eig(A);

%%% Finding center manifold (purely imaginary eigenvalue) ... or maybe not
989
for ev_i = 1:size(eVals,1)
    if abs(real(eVals(ev_i,ev_i))) < 5e-16 && abs(imag(eVals(ev_i,ev_i))) > 5e-16 % If real() == 0 and imag() ~= 0
        ev_i = 3
        warning('eigenvalue not dynamic')
        eVal_CM = eVals(ev_i,ev_i);
        eVec_CM = eVecs(:,ev_i);
%         eVal_CM = eVals(3,3);
%         eVec_CM = eVecs(:,3);

        break
    end
end
ev_i

% ========================================================================
%%% Periodic Orbits
% ========================================================================
% ------------------------------------------------- 
%%% Perturbing the Lagrange-point state in the direction of the center manifold
% -------------------------------------------------
%%% Perturbing the state to find first periodic orbit and estimating period
X_pert = [xL; yL; zL; 0; 0; 0] + real(eVec_CM).*1e-4;
T_pert = 2*pi/imag(eVal_CM);

dS0 = 1e-5; % Tuning parameter
% dS0 = 1e-8; % Tuning parameter
%%% First Tangent Orbit
X_tang = real(eVec_CM)/norm(real(eVec_CM));
T_tang = 0;

%%% Stepping in direction of first tangent orbit
X_guess = X_pert + dS0*X_tang;
T_guess = T_pert + dS0*T_tang;

%%% Storing first periodic orbit
POs = zeros(7,n);
POs(:,1) = [X_guess; T_guess];

% ------------------------------------------------- 
%%% Iterating through orbits
% -------------------------------------------------
%%% Initializating
error_tol = 1e-10;
for PO_i = 1:n
%     PO_i
    
    error = 1;
    counter = 0;
    %%% Finding new periodic orbit
    while error > error_tol
        %%% Preparing for integration
        % Initialize STM
        stm0 = eye(6);
        % Store full initial state
        X = [X_guess; reshape(stm0,36,1)];
        % Time vector for integration
        time = [0,T_guess];
        
        %%% Integrating
        [T, X] = ode45(@int_CR3BnSTM_J2, time, X, options, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);
        
        %%% Retreiving updated STM (monodromy matrix guess)
        stm = reshape(X(end,7:end),6,6);
        
        %%% Correcting the Initial Guess (crazy shit)
        F = [X(end,1:6)-X(1,1:6), (X(1, 1:6)' - X_pert)'*CR3Bn_EOM_J2(X_pert,secondary.MR,primary.R/rNorm,secondary.R_n,primary.J2,0), (X(1, 1:6)' - X_pert)'*X_tang + (T_guess - T_pert)*T_tang - dS_PO]';
        D = [stm-eye(6), CR3Bn_EOM_J2(X_guess,secondary.MR,primary.R/rNorm,secondary.R_n,primary.J2,0);...
             CR3Bn_EOM_J2(X_pert,secondary.MR,primary.R/rNorm,secondary.R_n,primary.J2,0)', 0;...
             X_tang', T_tang];
        delta = (D'*D)\D'*(-F);
        
        X_guess = X_guess + delta(1:6);
        T_guess = T_guess + delta(7);
        
        %%% Storing error
        error = norm(F);
        
        %%% Counting error tries
        counter = counter + 1;
    end
        
    %%% Storing new periodic orbit
    POs(:,PO_i+1) = [X_guess; T_guess];

    %%% Update perturbed orbit
    X_pert = X_guess;
    T_pert = T_guess;

    %%% Compute new tangent orbit
    tangentOrbit = (POs(:,PO_i+1) - POs(:,PO_i))/norm(POs(:,PO_i+1)-POs(:,PO_i));

    %%% Update the tangent orbit
    X_tang = tangentOrbit(1:6);
    T_tang = tangentOrbit(7);

    %%% Update guess for first tangent orbit
    X_guess = X_pert + dS_PO*X_tang;
    T_guess = T_pert + dS_PO*T_tang;

    %%% Occasionally plotting periodic orbit
    if mod(PO_i,PO_plot_skip) == 0
        PO_i
        figure(1); hold all; grid on;
        p1 = plot3(X(:,1),X(:,2),X(:,3),'m');
        
    end
    
    %%% Determine if tangent orbit step size should change (or some weird shit like that)
    if dynamic_dS_PO == 1
        if counter > 5
            dS_PO = dS_PO/counter;
            if dS_PO < 1e-5
                dS_PO = 1e-5;
            elseif dS_PO > 1e-2
                dS_PO = 1e-2;
            end

        elseif counter <= 5
            dS_PO = 1.2*dS_PO;
            if dS_PO < 1e-5
                dS_PO = 1e-5;
            elseif dS_PO > 1e-2
                dS_PO = 1e-2;
            end 
        end
    end
    
end

%%% Adding detail to plot
figure(1); hold all; grid on;

if isfield(secondary,'img') == 1
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
else
    plotBody3(secondary.R_n,[1-secondary.MR,0,0],secondary.color)
end
% if isfield(primary,'img') == 1
%     plotBodyTexture3(primary.R/rNorm,[-secondary.MR,0,0],primary.img)
% else
%     plotBody3(primary.R/rNorm,[-secondary.MR,0,0],primary.color)
% end
plot3(xL,yL,zL,'^','markeredgecolor',colors.std.black,'markerfacecolor','m')
% legend([p1],'With J2')
PlotBoi3('X_n','Y_n','Z_n',14)
view(0,90)
axis square



% [JC_L1]   = JacobiConstantCalculator_J2(secondary.MR,Ls_n(1,:),[0,0,0],primary.R/rNorm,secondary.R_n,primary.J2,0);
% [JC_L2]   = JacobiConstantCalculator_J2(secondary.MR,Ls_n(2,:),[0,0,0],primary.R/rNorm,secondary.R_n,primary.J2,0);
% [JC_Surf] = JacobiConstantCalculator_J2(secondary.MR,[1-secondary.MR-secondary.R_n,0,0],[0,0,0],primary.R/rNorm,secondary.R_n,primary.J2,0);
% 
% JCs_PO = zeros(1,size(POs,2));
% for kk = 1:size(POs,2)
%     [JC] = JacobiConstantCalculator_J2(secondary.MR,POs(1:3,kk)',POs(4:6,kk)',primary.R/rNorm,secondary.R_n,primary.J2,0);
%     JCs_PO(kk) = JC;
% end
%     
% figure(2); hold all
% plot(1:size(POs,2),JCs_PO,'o','markersize',5)


















% ========================================================================
%%% Functions
% ========================================================================

function [dX] = CR3Bn_EOM_J2(X,u,R1,R2,J21,J22)

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Define position of bodies
x1 = -u;
x2 = 1-u;

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x+u)^2 + y^2 + z^2);
r2 = sqrt((x+u-1)^2 + y^2 + z^2);


%%% Define Gamma values for shortening J2 terms
gamma1 = 3*(1-u)*R1^2*J21*(5*z^2-r1^2)/(2*r1^7);
gamma2 = 3*u*R2^2*J22*(5*z^2-r2^2)/(2*r2^7);
gamma1_z = 3*(1-u)*R1^2*J21*(5*z^2-3*r1^2)/(2*r1^7);
gamma2_z = 3*u*R2^2*J22*(5*z^2-3*r2^2)/(2*r2^7);

%%% Equations of Motion
ddx = 2*dy + x + ((1-u)/r1^3-gamma1)*(x1-x) + (u/r2^3 - gamma2)*(x2-x);
ddy = -2*dx + y*(-(1-u)/r1^3 - u/r2^3 + gamma1 + gamma2 + 1);
ddz = z*(-(1-u)/r1^3 - u/r2^3 + gamma1_z + gamma2_z);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end





