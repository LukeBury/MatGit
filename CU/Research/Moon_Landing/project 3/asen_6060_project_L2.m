%% ASEN 6060 Project L2
% Don Kuettel

clc; clear; close all

tic

% ODE45 Tolerances
myoptions = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Constants
G = 6.673e-20;
R = 384400;
dS = 1e-5;

M1 = 5.972e24;      % kg (Earth)
M2 = 7.35e22;       % kg (Moon)

n = sqrt((G*(M1+M2))/R^3);

v = M2/(M1+M2);
R1 = -v;
R2 = 1 - v;

%% Finding Lagrange points

x1 = 0.8;
y1 = 0;

x2 = 1.15;
y2 = 0;

x3 = -1;
y3 = 0;

x4 = 0.5;
y4 = 0.8;

x5 = 0.5;
y5 = -0.8;

[x_L1, y_L1] = CR3BP_newton_raph(x1, y1, v);
[x_L2, y_L2] = CR3BP_newton_raph(x2, y2, v);
[x_L3, y_L3] = CR3BP_newton_raph(x3, y3, v);
[x_L4, y_L4] = CR3BP_newton_raph(x4, y4, v);
[x_L5, y_L5] = CR3BP_newton_raph(x5, y5, v);

LPs_norm = [x_L1 y_L1; x_L2 y_L2; x_L3 y_L3; x_L4 y_L4; x_L5 y_L5];
LPs_dim = R*[x_L1 y_L1; x_L2 y_L2; x_L3 y_L3; x_L4 y_L4; x_L5 y_L5];

z_L1 = 0; z_L2 = 0; z_L3 = 0; z_L4 = 0; z_L5 = 0;

%% Finding Initial Orbit around L1

%  Evaluate the partial matrix A(t)
A_L2 = zeros(6,6);

A_L2(1,4) = 1; A_L2(2,5) = 1; A_L2(3,6) = 1; A_L2(4,5) = 2; A_L2(5,4) = -2;

A_L2(4,1) = (v - 1)/((v + x_L2)^2 + y_L2^2 + z_L2^2)^(3/2) - v/((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(3/2) + (3*v*(2*v + 2*x_L2 - 2)^2)/(4*((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(5/2)) - (3*(2*v + 2*x_L2)^2*(v - 1))/(4*((v + x_L2)^2 + y_L2^2 + z_L2^2)^(5/2)) + 1;
A_L2(4,2) = (3*v*y_L2*(2*v + 2*x_L2 - 2))/(2*((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(5/2)) - (3*y_L2*(2*v + 2*x_L2)*(v - 1))/(2*((v + x_L2)^2 + y_L2^2 + z_L2^2)^(5/2));
A_L2(4,3) = (3*v*z_L2*(2*v + 2*x_L2 - 2))/(2*((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(5/2)) - (3*z_L2*(2*v + 2*x_L2)*(v - 1))/(2*((v + x_L2)^2 + y_L2^2 + z_L2^2)^(5/2));
A_L2(5,1) = (3*v*y_L2*(2*v + 2*x_L2 - 2))/(2*((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(5/2)) - (3*y_L2*(2*v + 2*x_L2)*(v - 1))/(2*((v + x_L2)^2 + y_L2^2 + z_L2^2)^(5/2));
A_L2(5,2) = (v - 1)/((v + x_L2)^2 + y_L2^2 + z_L2^2)^(3/2) - v/((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(3/2) - (3*y_L2^2*(v - 1))/((v + x_L2)^2 + y_L2^2 + z_L2^2)^(5/2) + (3*v*y_L2^2)/((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(5/2) + 1;
A_L2(5,3) = (3*v*y_L2*z_L2)/((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(5/2) - (3*y_L2*z_L2*(v - 1))/((v + x_L2)^2 + y_L2^2 + z_L2^2)^(5/2);
A_L2(6,1) = (3*v*z_L2*(2*v + 2*x_L2 - 2))/(2*((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(5/2)) - (3*z_L2*(2*v + 2*x_L2)*(v - 1))/(2*((v + x_L2)^2 + y_L2^2 + z_L2^2)^(5/2));
A_L2(6,2) = (3*v*y_L2*z_L2)/((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(5/2) - (3*y_L2*z_L2*(v - 1))/((v + x_L2)^2 + y_L2^2 + z_L2^2)^(5/2);
A_L2(6,3) = (v - 1)/((v + x_L2)^2 + y_L2^2 + z_L2^2)^(3/2) - v/((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(3/2) - (3*z_L2^2*(v - 1))/((v + x_L2)^2 + y_L2^2 + z_L2^2)^(5/2) + (3*v*z_L2^2)/((v + x_L2 - 1)^2 + y_L2^2 + z_L2^2)^(5/2);

[Evecs, Evals] = eig(A_L2);

% Compute First Periodic Orbit
X_tilda = [x_L2; y_L2; 0; 0; 0; 0] + 1e-4*real(Evecs(:,3));
T_tilda = 2*pi/imag(Evals(3,3));

% Compute First Orbit Tangent
X_tangent = real(Evecs(:,3))/norm(real(Evecs(:,3)));
T_tangent = 0;

% Compute first Orbit guess
X_in = X_tilda + dS*X_tangent;
T_in = T_tilda + dS*T_tangent;

PO_Matrix_IC(:,1) = [X_in; T_in];

%% Iteration of Orbits

index = 1;
counter = 0;
for k = 1:1000;
    k
    if counter > 5
        dS = dS/counter;
        
        if dS < 1e-4
            dS = 1e-4;
            
        elseif dS > 1e-2
            dS = 1e-2;
            
        end
        
    elseif counter <= 5
        dS = 1.1*dS;
        
        if dS < 1e-4
            dS = 1e-4;
            
        elseif dS > 1e-2
            dS = 1e-2;
            
        end
        
    end
    
    err = 1;
    counter = 0;
    % Newton Iteration
    while err > 1e-10
        
        % Initial State
        Phi0 = eye(6);                    %  I.C. for the STM
        
        % Store all initial conditions in a single vector.  To do this
        % we need to reshape the initial STM to a 6^2x1 vector
        IC = [X_in' v reshape(Phi0,1,6*6)];
        
        time = [0 T_in];
        
        % Now, we call ode45() to numerically solve the EoM
        [T, X] = ode45(@CR3BP_norm, time, IC, myoptions);
        
        % Breaking down state matrix
        x = X(:,1);
        y = X(:,2);
        z = X(:,3);
        vx = X(:,4);
        vy = X(:,5);
        vz = X(:,6);
        
        STM = reshape(X(end,8:end),6,6);
        
        % Correct the Initial Guess
        F = [(X(end, 1:6)' - X(1, 1:6)')' (X(1, 1:6)' - X_tilda)'*CR3BPdynamics(X_tilda,v) (X(1, 1:6)' - X_tilda)'*X_tangent + (T_in - T_tilda)*T_tangent - dS]';
        
        D = [STM - eye(6) CR3BPdynamics(X_in,v); CR3BPdynamics(X_tilda,v)' 0; X_tangent' T_tangent];
        
        delta = (D'*D)\D'*-F;
        
        X_in = X_in + delta(1:6);
        T_in = T_in + delta(7);
        
        err = norm(F);
        
        counter = counter + 1;
        
    end
    
    dS_tally(k) = dS;
    
    if mod(k,10) == 0 && k >= 20;
        % Stats
        PO_Matrix_IC_plot(:,index) = [X_in; T_in];
        lamda = abs(eig(STM));
        stab_index(index) = (real(max(lamda)) + real(1/max(lamda)))/2;
        
        figure(1)
        hold on; grid on; box on
        h_earth = plot(R1, 0, 'ok');
        h_moon = plot(R2, 0, 'ok');
        h_L1 = plot(x_L2,y_L2,'or');
        plot(x, y,'-k')
        set(h_earth,'MarkerEdgeColor','k','MarkerFaceColor','b','markersize',10)
        set(h_moon,'MarkerEdgeColor','k','MarkerFaceColor','k','markersize',5)
        set(h_L1,'MarkerEdgeColor','k','MarkerFaceColor','r','markersize',5)
        hold off
        xlabel('x [-]')
        ylabel('y [-]')
        legend('Earth','Moon','L_2','Lyapunov Orbits','location','best')
        title('L2 Planar Lyapunov Orbits')
        axis('square')
        
        index = index + 1;
    end
    
    % PO_Matrix_coords(:,6*k-5:6*k) = X(:,1:6);
    PO_Matrix_IC(:,k+1) = [X_in; T_in];
    
    % Update tilda
    X_tilda = X_in;
    T_tilda = T_in;
    
    % Compute Orbit Tangent
    Tan = (PO_Matrix_IC(:,k+1) - PO_Matrix_IC(:,k))/norm(PO_Matrix_IC(:,k+1) - PO_Matrix_IC(:,k));
    
    % Update Tangent
    X_tangent = Tan(1:6);
    T_tangent = Tan(7);
    
    % Update first Orbit guess
    X_in = X_tilda + dS*X_tangent;
    T_in = T_tilda + dS*T_tangent;
    
end

toc

%%
figure
hold on; grid on; box on
plot(PO_Matrix_IC_plot(1,:),stab_index)
hold off
xlabel('x [-]')
ylabel('Stability Index')
title('L_2 Stability Index')