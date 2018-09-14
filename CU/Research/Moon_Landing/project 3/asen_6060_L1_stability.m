%% ASEN 6060 Project - Stability
% Don Kuettel

clc; clear; close all

tic

% ODE45 Tolerances
myoptions = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Constants
G = 6.673e-20;
R = 384400;

M1 = 5.972e24;      % kg (Earth)
M2 = 7.35e22;       % kg (Moon)

% v = M2/(M1+M2);
v = 2.5281e-05; % europa
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

%% Full Period

% From L1 script
X_in = [0.808769550768705;-2.84553528821764e-05;0;-2.91048806199374e-05;0.284396224284932;0];

P = 3.02385643730079;

% Initial State
Phi0 = eye(6);                    %  I.C. for the STM

% Store all initial conditions in a single vector.  To do this
% we need to reshape the initial STM to a 6^2x1 vector
IC = [X_in' v reshape(Phi0,1,6*6)];

time = [0 P];

% Now, we call ode45() to numerically solve the EoM
[T, X] = ode45(@CR3BP_norm, time, IC, myoptions);

% Breaking down state matrix
x_PO = X(:,1);
y_PO = X(:,2);
z_PO = X(:,3);
vx_PO = X(:,4);
vy_PO = X(:,5);
vz_PO = X(:,6);

STM = reshape(X(end,8:end),6,6);

[Evecs, Evals] = eig(STM);

% Unstable Manifold Conditions
X_in_U = X_in + 5e-8*Evecs(:,2);


%% Integration
n = 300;

for k = 1:n
    
    % Initial State
    Phi0 = eye(6);                    %  I.C. for the STM
    
    % Store all initial conditions in a single vector.  To do this
    % we need to reshape the initial STM to a 6^2x1 vector
    IC = [X_in' v reshape(Phi0,1,6*6)];
    
    time = [0 P/n*k];
    
    % Now, we call ode45() to numerically solve the EoM
    [T, X] = ode45(@CR3BP_norm, time, IC, myoptions);
    
    % Breaking down state matrix
    x = X(:,1);
    y = X(:,2);
    z = X(:,3);
    vx = X(:,4);
    vy = X(:,5);
    vz = X(:,6);
    
    X_partial = [x(end); y(end); z(end); vx(end); vy(end); vz(end)];
    
    STM = reshape(X(end,8:end),6,6);
    
    %% Unstable Exterior Manifold
    X_in_U_partial = X_partial - 5e-8*(STM*Evecs(:,2)/norm(STM*Evecs(:,2)));
    
    % Initial State
    Phi0 = eye(6);                    %  I.C. for the STM
    
    % Store all initial conditions in a single vector.  To do this
    % we need to reshape the initial STM to a 6^2x1 vector
    IC = [X_in_U_partial' v reshape(Phi0,1,6*6)];
            
    time = [P*2.5 0];
%     time = [0 P*2.5];
    
    % Now, we call ode45() to numerically solve the EoM
    [T, X] = ode45(@CR3BP_norm, time, IC, myoptions);
    
    % Breaking down state matrix
    x = X(:,1);
    y = X(:,2);
    z = X(:,3);
    vx = X(:,4);
    vy = X(:,5);
    vz = X(:,6);
    
    figure(1)
    hold on; grid on; box on
    h_earth = plot(R1, 0, 'ok');
    h_moon = plot(R2, 0, 'ok');
    h_L1 = plot(x_L1,y_L1,'or');
    plot(x, y,'-k')
    set(h_earth,'MarkerEdgeColor','k','MarkerFaceColor','b','markersize',10)
    set(h_moon,'MarkerEdgeColor','k','MarkerFaceColor','k','markersize',5)
    set(h_L1,'MarkerEdgeColor','k','MarkerFaceColor','r','markersize',5)
    hold off
    xlabel('x [-]')
    ylabel('y [-]')
    legend('Earth','Moon','L_1','Stable Orbits','location','best')
    title('L1 Stable Manifold')
    axis('square')
    
end

toc
%%
figure
hold on; grid on; box on
h_earth = plot(R1, 0, 'ok');
h_moon = plot(R2, 0, 'ok');
h_L1 = plot(x_L1,y_L1,'or');
plot(x_PO, y_PO,'-k')
set(h_earth,'MarkerEdgeColor','k','MarkerFaceColor','b','markersize',10)
set(h_moon,'MarkerEdgeColor','k','MarkerFaceColor','k','markersize',5)
set(h_L1,'MarkerEdgeColor','k','MarkerFaceColor','r','markersize',5)
hold off
xlabel('x [-]')
ylabel('y [-]')
legend('Earth','Moon','L_1','Orbit','location','best')
title('L1 Periodic Orbit')
axis('square')