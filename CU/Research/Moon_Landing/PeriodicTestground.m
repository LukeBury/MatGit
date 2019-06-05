clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Run switches
% ========================================================================
run_symbollicWork = 0;
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Symbolic Work
% ========================================================================
if run_symbollicWork == 1
syms x y z R1 R2 J21 J22 u dx dy dz real

r1 = sqrt((-u-xL)^2+yL^2+zL^2);
r2 = sqrt((1-u-xL)^2+yL^2+zL^2);

%%% Defining Equations of Motion, State Variables, and A Matrix
ddx = 2*dy + xL - (1-u)*(xL+u)/(r1^3) - u*(xL+u-1)/(r2^3);
ddy = -2*dx + yL -((1-u)/(r1^3) + u/(r2^3))*yL;
ddz = -((1-u)/(r1^3) + u/(r2^3))*zL;

EOM = [dx; dy; dz; ddx; ddy; ddz];
state = [xL; yL; zL; dx; dy; dz];
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
n = 1000;
L_Point = 2;
dS_PO = .0001;
% dS_PO = .000001;
dynamic_dS_PO = 0;
PO_plot_skip = 50;

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-9;
options = odeset('RelTol',tol,'AbsTol',tol);
prms.u = secondary.MR;

% ------------------------------------------------- 
%%% Finding Equilibrium Points
% -------------------------------------------------
Ls_n = EquilibriumPoints(secondary.MR);

% ========================================================================
%%% Evaluating A matrix at equilibrium point
% ========================================================================
% ------------------------------------------------- 
%%% Choosing equilibrium point and setting initial state
% -------------------------------------------------
%%% Coordinates of Lagrange Point to be evaluated
xL = Ls_n(L_Point,1);
yL = Ls_n(L_Point,2);
zL = Ls_n(L_Point,3);

% ------------------------------------------------- 
%%% Evaluating A at equilibrium point
% -------------------------------------------------
A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4,5) = 2;
A(5,4) = -2;
A(4,1) = (u - 1)/((u + xL)^2 + yL^2 + zL^2)^(3/2) - u/((u + xL - 1)^2 + yL^2 + zL^2)^(3/2) + (3*u*(2*u + 2*xL - 2)*(u + xL - 1))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(5/2)) - (3*(2*u + 2*xL)*(u + xL)*(u - 1))/(2*((u + xL)^2 + yL^2 + zL^2)^(5/2)) + 1;
A(4,2) = (3*u*yL*(u + xL - 1))/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*yL*(u + xL)*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2);
A(4,3) = (3*u*zL*(u + xL - 1))/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*zL*(u + xL)*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2);
A(5,1) = -yL*((3*(2*u + 2*xL)*(u - 1))/(2*((u + xL)^2 + yL^2 + zL^2)^(5/2)) - (3*u*(2*u + 2*xL - 2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(5/2)));
A(5,2) = (u - 1)/((u + xL)^2 + yL^2 + zL^2)^(3/2) + yL*((3*u*yL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*yL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2)) - u/((u + xL - 1)^2 + yL^2 + zL^2)^(3/2) + 1;
A(5,3) = yL*((3*u*zL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*zL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2));
A(6,1) = -zL*((3*(2*u + 2*xL)*(u - 1))/(2*((u + xL)^2 + yL^2 + zL^2)^(5/2)) - (3*u*(2*u + 2*xL - 2))/(2*((u + xL - 1)^2 + yL^2 + zL^2)^(5/2)));
A(6,2) = zL*((3*u*yL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*yL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2));
A(6,3) = (u - 1)/((u + xL)^2 + yL^2 + zL^2)^(3/2) + zL*((3*u*zL)/((u + xL - 1)^2 + yL^2 + zL^2)^(5/2) - (3*zL*(u - 1))/((u + xL)^2 + yL^2 + zL^2)^(5/2)) - u/((u + xL - 1)^2 + yL^2 + zL^2)^(3/2);

% ------------------------------------------------- 
%%% Eigenvectors and Eigenvalues
% -------------------------------------------------
%%% Getting eigenvectors and eigenvalues
[eVecs, eVals] = eig(A);

%%% Finding center manifold (purely imaginary eigenvalue) ... or maybe not
for ev_i = 1:size(eVals,1)
    if abs(real(eVals(ev_i,ev_i))) < 5e-16 && abs(imag(eVals(ev_i,ev_i))) > 5e-16 % If real() == 0 and imag() ~= 0
        ev_i = 5
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
    PO_i
    
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
        [T, X] = ode45(@Int_CR3BnSTM, time, X, options, prms);
        
        %%% Retreiving updated STM
        stm = reshape(X(end,7:end),6,6);
        
        %%% Correcting the Initial Guess (crazy shit)
        F = [X(end,1:6)-X(1,1:6), (X(1, 1:6)' - X_pert)'*CR3Bn_EOM(X_pert,secondary.MR), (X(1, 1:6)' - X_pert)'*X_tang + (T_guess - T_pert)*T_tang - dS_PO]';
        D = [stm-eye(6), CR3Bn_EOM(X_guess,secondary.MR);...
             CR3Bn_EOM(X_pert,secondary.MR)', 0;...
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

    %%% Should this orbit be plotted?
    if mod(PO_i,PO_plot_skip) == 0
        figure(1); hold all; grid on;
        p2 = plot3(X(:,1),X(:,2),X(:,3),'k');
        
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
plot3(xL,yL,zL,'^','markeredgecolor',colors.std.black,'markerfacecolor','k')
% legend([p2],'Basic')
PlotBoi3('X_n','Y_n','Z_n',14)
view(0,90)
axis square




[JC_L1]   = JacobiConstantCalculator(secondary.MR,Ls_n(1,:),[0,0,0])
[JC_L2]   = JacobiConstantCalculator(secondary.MR,Ls_n(2,:),[0,0,0])
[JC_sNX] = JacobiConstantCalculator(secondary.MR,[1-secondary.MR-secondary.R_n,0,0],[0,0,0])
% 
JCs_PO = zeros(1,size(POs,2));
for kk = 1:size(POs,2)
    [JC] = JacobiConstantCalculator(secondary.MR,POs(1:3,kk)',POs(4:6,kk)');
    JCs_PO(kk) = JC;
end
    
figure(2); hold all
plot(1:size(POs,2),JCs_PO,'o','markersize',5)





















% ========================================================================
%%% Functions
% ========================================================================

function [dX] = CR3Bn_EOM(X,u)

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x+u)^2 + y^2 + z^2);
r2 = sqrt((x+u-1)^2 + y^2 + z^2);

%%% Equations of Motion
ddx = 2*dy + x - (1-u)*(x+u)/(r1^3) - u*(x+u-1)/(r2^3);
ddy = -2*dx + y -((1-u)/(r1^3) + u/(r2^3))*y;
ddz = -((1-u)/(r1^3) + u/(r2^3))*z;

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz];
dX(4:6) = [ddx; ddy; ddz];

end





