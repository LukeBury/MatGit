clear
clc
% close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();



%% =======================================================================
%%% PO With J2
% ========================================================================



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
PO_n = 2000;
L_Point = 2;
dS_PO = .0001;
dynamic_dS_PO = 0;
PO_plot_skip = 399; 
ev_force = 3;

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
J22 = 0;

% ------------------------------------------------- 
%%% Evaluating A at equilibrium point
% -------------------------------------------------
[A] = get_Amat_CR3BP_J2(u, [xL, yL, zL], J21, J22, R1, R2);

% ------------------------------------------------- 
%%% Eigenvectors and Eigenvalues
% -------------------------------------------------
%%% Getting eigenvectors and eigenvalues
[eVecs, eVals] = eig(A);

%%% Finding center manifold (purely imaginary eigenvalue) ... or maybe not

for ev_i = 1:size(eVals,1)
    if abs(real(eVals(ev_i,ev_i))) < 5e-16 && abs(imag(eVals(ev_i,ev_i))) > 5e-16 % If real() == 0 and imag() ~= 0
%         ev_i = 3
%         warning('eigenvalue not dynamic')
        eVal_CM = eVals(ev_force,ev_force);
        eVec_CM = eVecs(:,ev_force);
%         eVal_CM = eVals(3,3);
%         eVec_CM = eVecs(:,3);

        break
    end
end


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
POs_J2 = zeros(7,PO_n);
POs_J2(:,1) = [X_guess; T_guess];

% ------------------------------------------------- 
%%% Iterating through orbits
% -------------------------------------------------
%%% Initializating
error_tol = 1e-10;
for PO_i = 1:(PO_n-1)
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
        [T, X] = ode45(@Int_CR3BnSTM_J2, time, X, options, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);
        
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
    POs_J2(:,PO_i+1) = [X_guess; T_guess];

    %%% Update perturbed orbit
    X_pert = X_guess;
    T_pert = T_guess;

    %%% Compute new tangent orbit
    tangentOrbit = (POs_J2(:,PO_i+1) - POs_J2(:,PO_i))/norm(POs_J2(:,PO_i+1)-POs_J2(:,PO_i));

    %%% Update the tangent orbit
    X_tang = tangentOrbit(1:6);
    T_tang = tangentOrbit(7);

    %%% Update guess for first tangent orbit
    X_guess = X_pert + dS_PO*X_tang;
    T_guess = T_pert + dS_PO*T_tang;

    %%% Occasionally plotting periodic orbit
    if mod(PO_i,PO_plot_skip) == 0
        figure(1); hold all; grid on; % [245 463 805 270], [245 494 609 239]
        p1 = plot3(X(:,1),X(:,2),X(:,3),'m','linewidth',1.5);
        
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
PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
view(0,90)
axis square




%% =======================================================================
%%% No J2
% ========================================================================

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
A = get_Amat_CR3BP(u, [xL, yL, zL]);

% ------------------------------------------------- 
%%% Eigenvectors and Eigenvalues
% -------------------------------------------------
%%% Getting eigenvectors and eigenvalues
[eVecs, eVals] = eig(A);

%%% Finding center manifold (purely imaginary eigenvalue) ... or maybe not
for ev_i = 1:size(eVals,1)
    if abs(real(eVals(ev_i,ev_i))) < 5e-16 && abs(imag(eVals(ev_i,ev_i))) > 5e-16 % If real() == 0 and imag() ~= 0
%         ev_i = 3
%         warning('eigenvalue not dynamic')
        eVal_CM = eVals(ev_force,ev_force);
        eVec_CM = eVecs(:,ev_force);
%         eVal_CM = eVals(3,3);
%         eVec_CM = eVecs(:,3);

        break
    end
end

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
POs = zeros(7,PO_n);
POs(:,1) = [X_guess; T_guess];

% ------------------------------------------------- 
%%% Iterating through orbits
% -------------------------------------------------
%%% Initializating
error_tol = 1e-10;
for PO_i = 1:(PO_n-1)
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
        prms.u = secondary.MR;
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
        p2 = plot3(X(:,1),X(:,2),X(:,3),'k','linewidth',1.5);
        
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
PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
view(0,90)
axis equal



legend([p1,p2],'With J_{2p}','Without J_{2p}')



JCs_PO = zeros(1,size(POs,2));
JCs_PO_J2 = zeros(1,size(POs_J2,2));
for kk = 1:PO_n
    [JC] = JacobiConstantCalculator(secondary.MR,POs(1:3,kk)',POs(4:6,kk)');
    [JC_J2] = JacobiConstantCalculator_J2(secondary.MR,POs_J2(1:3,kk)',POs_J2(4:6,kk)',primary.R/rNorm,secondary.R_n,primary.J2,0);
    
    JCs_PO(kk) = JC;
    JCs_PO_J2(kk) = JC_J2;
end
    
figure(2); hold all
plot([1:PO_n],JCs_PO,'mo','markersize',5)
plot([1:PO_n],JCs_PO_J2,'ko','markersize',5)

% JacobiConstantCalculator_J2(u,rBCR_n,vBCR_n, R1_n, R2_n, J21, J22)



















%% =======================================================================
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