clear
clc
close all
addpath('../../bin')

% ------------------------------------------------------------------------
%%% 2b
% ------------------------------------------------------------------------
syms x y z dx dy dz u ae J2 J3

U = -(3*u*ae*ae*J2*z*z)/(2*((x^2+y^2+z^2)^(5/2)))...
    + (u*ae*ae*J2)/(2*((x^2+y^2+z^2)^(3/2)))...
    - (5*u*ae*ae*ae*J3*z*z*z)/(2*((x^2+y^2+z^2)^(7/2)))...
    + (3*u*ae*ae*ae*J3*z)/(2*((x^2+y^2+z^2)^(5/2)));



dudx = diff(U,x); % 1x1

dudy = diff(U,y); % 1x1

dudz = diff(U,z); % 1x1

acc = [dudx; dudy; dudz]; % 3x1

dadx = simplify(diff(acc,x)); % 3x1
dady = simplify(diff(acc,y)); % 3x1
dadz = simplify(diff(acc,z)); % 3x1
dadu = simplify(diff(acc,u)); % 3x1
dadJ2 = simplify(diff(acc,J2)); % 3x1
dadJ3 = simplify(diff(acc,J3)); % 3x1

% fprintf('dadx \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadx(1),dadx(2),dadx(3))
% fprintf('dady \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dady(1),dady(2),dady(3))
% fprintf('dadz \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadz(1),dadz(2),dadz(3))
% fprintf('dadu \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadu(1),dadu(2),dadu(3))
% fprintf('dadJ2 \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadJ2(1),dadJ2(2),dadJ2(3))
% fprintf('dadJ3 \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadJ3(1),dadJ3(2),dadJ3(3))



%%% Don's EOM Method
r = sqrt(x^2 + y^2 + z^2);
U_J2 = -(3*u*ae*ae*J2*z*z)/(2*(r^(5))) + (u*ae*ae*J2)/(2*(r^(3))); % Potential from J2

state = [x; y; z; dx; dy; dz];

Utot_J2 = -u/r + U_J2; % Total potential

EQM = [dx; dy; dz; diff(Utot_J2, x); diff(Utot_J2, y); diff(Utot_J2,z)];
EQM(4)

Asym = jacobian(EQM, state);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
%%% Getting J2 & u accelerations
% ------------------------------------------------------------------------

tl = zeros(3,3);
tm = eye(3,3);
tr = zeros(3,3);
ml = [dadx(1), dady(1), dadz(1)
    dadx(2), dady(2), dadz(2)
    dadx(3), dady(3), dadz(3)];
mm = zeros(3,3);
mr = [simplify(diff(dudx,u)), simplify(diff(dudx,J2)), simplify(diff(dudx,ae))
    simplify(diff(dudy,u)), simplify(diff(dudy,J2)), simplify(diff(dudy,ae))
    simplify(diff(dudz,u)), simplify(diff(dudz,J2)), simplify(diff(dudz,ae))];
bl = zeros(3,3);
bm = zeros(3,3);
br = zeros(3,3);

stm = [tl tm tr; ml mm mr; bl bm br];
stm2 = [tl tm; ml mm];

stm0 = eye(6,6);

% ------------------------------------------------------------------------
%%% Getting J2 & u accelerations
% ------------------------------------------------------------------------
a1 = dudx + (15*J3*ae^3*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*ae^3*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
a2 = dudy + (15*J3*ae^3*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*ae^3*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
a3 = dudz - (3*J3*ae^3*u)/(2*(x^2 + y^2 + z^2)^(5/2)) + (15*J3*ae^3*u*z^2)/(x^2 + y^2 + z^2)^(7/2) - (35*J3*ae^3*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2));
a_3 = [a1; a2; a3];

% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
%%% Orbital Elements
a = 7000; % km
e = 0.001; 
i = 30*pi/180; % rad
RAAN = 80*pi/180; % rad
w = 40*pi/180; % rad
ta = 0; % rad

%%% State deviation vector
dx = [1; 0; 0; 0; .01; 0];% km, km/s

%%% Other
uE = 398600.4415; % km^3/s^2
J2 = 0.0010826269;
g = 9.81;
Re = 6378.1;

% ------------------------------------------------------------------------
%%% Finding Initial Position and Setting State
% ------------------------------------------------------------------------
%%% Converting OE to ECI state
[r0, v0] = OE2ECI(a, e, i, RAAN, w, ta, uE); % km, km/s

%%% Setting initial state vectors (6x1)
X0 = [r0; v0]

X02 = X0 + dx;

% ------------------------------------------------------------------------
%%% Propagating the State with Numerical Integration
% ------------------------------------------------------------------------
%%% Setting time frame
ti = 0; % sec
tf = 24*3600; % sec
dt = 1; % sec
time = ti:dt:tf; % sec

%%% Setting integrator accuracy
tol = 1E-10;
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating the State
[Times,States] = ode45(@StatOD_Hw1_Int,time,X0,options,uE,Re,J2);


%%% Propagating the State (X + dx)
[Times,States_pdx] = ode45(@StatOD_Hw1_Int,time,X02,options,uE,Re,J2);

numDevs = States_pdx - States; % Numerical Deviations
% ------------------------------------------------------------------------
%%% Propagating the State Transition Matrix
% ------------------------------------------------------------------------
%%% Preparing for STM propagation
IC = [X0; reshape(stm0,36,1)];

%%% Propagating the State (STM)
[Times,States_dx] = ode45(@StatOD_Hw1_STMInt,time,IC,options,uE,Re,J2,dadx,dady,dadz,dadJ2);

%%% Unpacking the STMS and adding the new dxs to the computed states
stms = zeros(6,6,length(Times));
stm_states = zeros(length(Times),6);
stm_Devs = zeros(length(Times),6);
for k = 1:length(Times)
    stms(:,:,k) = reshape(States_dx(k,7:end),6,6);
    stm_Devs(k,:) = stms(:,:,k)*dx;
    stm_states(k,:) = States_dx(k,1:6) + (stms(:,:,k)*dx)';
end


% ------------------------------------------------------------------------
%%% Plotting
% ------------------------------------------------------------------------
%%% Plotting nominal trajectory
figure
plot3(States(:,1),States(:,2),States(:,3))
title('Nominal Trajectory')
PlotBoi3('X, km', 'Y, km', 'Z, km', 14)
view(-60,15)

%%% Plotting perturbed trajectory
figure
plot3(States_pdx(:,1),States_pdx(:,2),States_pdx(:,3))
title('Perturbed Trajectory')
PlotBoi3('X, km', 'Y, km', 'Z, km', 14)
view(-60,15)

%%% Results from STM
figure
plot3(stm_states(:,1),stm_states(:,2),stm_states(:,3))
title('STM Trajectory')
PlotBoi3('X, km','Y, km', 'Z, km', 14)
view(-60,15)

%%% Comparing Nominal to numerically integrated perturbed
figure
subplot(3,2,1)
plot(Times./3600,States(:,1)-States_pdx(:,1))
title('Error between nominal and perturbed trajectories')
PlotBoi2('','X Error, km',14)
subplot(3,2,3)
plot(Times./3600,States(:,2)-States_pdx(:,2))
PlotBoi2('','Y Error, km',14)
subplot(3,2,5)
plot(Times./3600,States(:,3)-States_pdx(:,3))
PlotBoi2('Time, hr','Z Error, km',14)
subplot(3,2,2)
plot(Times./3600,States(:,4)-States_pdx(:,4))
PlotBoi2('','X-Dot Error, km',14)
subplot(3,2,4)
plot(Times./3600,States(:,5)-States_pdx(:,5))
PlotBoi2('','Y-Dot Error, km',14)
subplot(3,2,6)
plot(Times./3600,States(:,6)-States_pdx(:,6))
PlotBoi2('Time, hr','Z-Dot Error, km',14)

%%% Comparing STM to numerically integrated perturbed
figure
subplot(3,2,1)
plot(Times./3600,stm_states(:,1)-States_pdx(:,1))
title('Error between stm and perturbed trajectories')
PlotBoi2('','X Error, km',14)
subplot(3,2,3)
plot(Times./3600,stm_states(:,2)-States_pdx(:,2))
PlotBoi2('','Y Error, km',14)
subplot(3,2,5)
plot(Times./3600,stm_states(:,3)-States_pdx(:,3))
PlotBoi2('Time, hr','Z Error, km',14)
subplot(3,2,2)
plot(Times./3600,stm_states(:,4)-States_pdx(:,4))
PlotBoi2('','X-Dot Error, km',14)
subplot(3,2,4)
plot(Times./3600,stm_states(:,5)-States_pdx(:,5))
PlotBoi2('','Y-Dot Error, km',14)
subplot(3,2,6)
plot(Times./3600,stm_states(:,6)-States_pdx(:,6))
PlotBoi2('Time, hr','Z-Dot Error, km',14)

%%% Comparing STM to nominal
figure
subplot(3,2,1)
plot(Times./3600,stm_states(:,1)-States(:,1))
title('Error between stm and nominal trajectories')
PlotBoi2('','X Error, km',14)
subplot(3,2,3)
plot(Times./3600,stm_states(:,2)-States(:,2))
PlotBoi2('','Y Error, km',14)
subplot(3,2,5)
plot(Times./3600,stm_states(:,3)-States(:,3))
PlotBoi2('Time, hr','Z Error, km',14)
subplot(3,2,2)
plot(Times./3600,stm_states(:,4)-States(:,4))
PlotBoi2('','X-Dot Error, km',14)
subplot(3,2,4)
plot(Times./3600,stm_states(:,5)-States(:,5))
PlotBoi2('','Y-Dot Error, km',14)
subplot(3,2,6)
plot(Times./3600,stm_states(:,6)-States(:,6))
PlotBoi2('Time, hr','Z-Dot Error, km',14)

%%% Comparing STM to nominal
figure
subplot(3,1,1)
hold all
plot(Times./3600,numDevs(:,1))
plot(Times./3600,stm_Devs(:,1))
subplot(3,1,2)
hold all
plot(Times./3600,numDevs(:,2))
plot(Times./3600,stm_Devs(:,2))
subplot(3,1,3)
hold all
plot(Times./3600,numDevs(:,3))
plot(Times./3600,stm_Devs(:,3))


devDiff = stm_Devs - numDevs;
%%% Comparing STM to nominal
figure
subplot(3,1,1)
plot(Times./3600,devDiff(:,1))
PlotBoi2('','x %-Change',14)
subplot(3,1,2)
plot(Times./3600,devDiff(:,2))
PlotBoi2('','y %-Change',14)
subplot(3,1,3)
plot(Times./3600,devDiff(:,3))
PlotBoi2('','z %-Change',14)

% (JC_B-JC_B(1))*100./(JC_B(1))

% ------------------------------------------------------------------------
%%% Numerical Integrators
% ------------------------------------------------------------------------
function [ dY ] = StatOD_Hw1_Int(t,Y,u,ae,J2)
dY = zeros(6,1);
J3 = 0;

%%% Unpack the state vector (ECI)
x = Y(1);
y = Y(2);
z = Y(3);
dy = Y(4:6); % Satellite Velocity, km/s

%%% Output the derivative of the state
dY(1:3) = dy; % km/s
% dudx - u*x/r^3
% dY(4) = -u*x/(norm([x,y,z])^3) + (15*J2*ae^2*u*x*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*ae^3*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*J2*ae^2*u*x)/(2*(x^2 + y^2 + z^2)^(5/2)) + (35*J3*ae^3*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)); % km/s^2
% EQM(4)
dY(4) = (3*J2*ae^2*u*x)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*x)/(x^2 + y^2 + z^2)^(3/2) - (15*J2*ae^2*u*x*z^2)/(2*(x^2 + y^2 + z^2)^(7/2));
% dudy - u*y/r^3
% dY(5) = -u*y/(norm([x,y,z])^3) + (15*J2*ae^2*u*y*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*ae^3*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*J2*ae^2*u*y)/(2*(x^2 + y^2 + z^2)^(5/2)) + (35*J3*ae^3*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)); % km/s^2
% EQM(5)
dY(5) = (3*J2*ae^2*u*y)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*y)/(x^2 + y^2 + z^2)^(3/2) - (15*J2*ae^2*u*y*z^2)/(2*(x^2 + y^2 + z^2)^(7/2));
% dudz - u*z/r^3
% dY(6) = -u*z/(norm([x,y,z])^3) + (3*J3*ae^3*u)/(2*(x^2 + y^2 + z^2)^(5/2)) + (15*J2*ae^2*u*z^3)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*ae^3*u*z^2)/(x^2 + y^2 + z^2)^(7/2) + (35*J3*ae^3*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2)) - (9*J2*ae^2*u*z)/(2*(x^2 + y^2 + z^2)^(5/2)); % km/s^2
% EQM(6)
dY(6) = (9*J2*ae^2*u*z)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*ae^2*u*z^3)/(2*(x^2 + y^2 + z^2)^(7/2)) - (u*z)/(x^2 + y^2 + z^2)^(3/2);
end



function [ dY ] = StatOD_Hw1_STMInt(t,Y,u,ae,J2,dadx,dady,dadz,dadJ2)
%%% Size dY to fit state and all of reshaped (n^2,1) STM
dY = zeros(6+6^2,1);

%%% Unpack state
x = Y(1);
y = Y(2);
z = Y(3);
dx = Y(4);
dy = Y(5);
dz = Y(6);

%%% Get rid of J3 terms
J3 = 0;

%%% Reshape (n^2,1) stm to (n,n)
stm = reshape(Y(7:end),6,6);

%%% Build A matrix and evaluate at current state
A = zeros(6,6);
A(1:3,4:6) = eye(3,3);
% dadx(1) ... (dax/dx)
% A(4,1) = -(3*J2*ae^2*u*(x^2 + y^2 + z^2)^3 - 15*J2*ae^2*u*z^2*(x^2 + y^2 + z^2)^2 + 315*J3*ae^3*u*x^2*z^3 + 15*J3*ae^3*u*z*(x^2 + y^2 + z^2)^2 - 35*J3*ae^3*u*z^3*(x^2 + y^2 + z^2) - 15*J2*ae^2*u*x^2*(x^2 + y^2 + z^2)^2 - 105*J3*ae^3*u*x^2*z*(x^2 + y^2 + z^2) + 105*J2*ae^2*u*x^2*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2));
% Asym(4,1)
A(4,1) = (3*u*x^2)/(x^2 + y^2 + z^2)^(5/2) - u/(x^2 + y^2 + z^2)^(3/2) + (3*J2*ae^2*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*ae^2*u*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*ae^2*u*x^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (105*J2*ae^2*u*x^2*z^2)/(2*(x^2 + y^2 + z^2)^(9/2));
% dady(1) ... (dax/dy)
% A(4,2) = (15*ae^2*u*x*y*(J2*x^4 + 2*J2*x^2*y^2 - 5*J2*x^2*z^2 + 7*J3*ae*x^2*z + J2*y^4 - 5*J2*y^2*z^2 + 7*J3*ae*y^2*z - 6*J2*z^4 - 14*J3*ae*z^3))/(2*(x^2 + y^2 + z^2)^(11/2));
% Asym(4,2)
A(4,2) = (3*u*x*y)/(x^2 + y^2 + z^2)^(5/2) - (15*J2*ae^2*u*x*y)/(2*(x^2 + y^2 + z^2)^(7/2)) + (105*J2*ae^2*u*x*y*z^2)/(2*(x^2 + y^2 + z^2)^(9/2));
% dadz(1) ... (dax/dz)
% A(4,3) = -(315*J3*ae^3*u*x*z^4 + 15*J3*ae^3*u*x*(x^2 + y^2 + z^2)^2 - 45*J2*ae^2*u*x*z*(x^2 + y^2 + z^2)^2 + 105*J2*ae^2*u*x*z^3*(x^2 + y^2 + z^2) - 210*J3*ae^3*u*x*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2));
% Asym(4,3)
A(4,3) = (3*u*x*z)/(x^2 + y^2 + z^2)^(5/2) - (45*J2*ae^2*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (105*J2*ae^2*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
% dadx(2) ... (day/dx)
% A(5,1) = (15*ae^2*u*x*y*(J2*x^4 + 2*J2*x^2*y^2 - 5*J2*x^2*z^2 + 7*J3*ae*x^2*z + J2*y^4 - 5*J2*y^2*z^2 + 7*J3*ae*y^2*z - 6*J2*z^4 - 14*J3*ae*z^3))/(2*(x^2 + y^2 + z^2)^(11/2));
% Asym(5,1)
A(5,1) = (3*u*x*y)/(x^2 + y^2 + z^2)^(5/2) - (15*J2*ae^2*u*x*y)/(2*(x^2 + y^2 + z^2)^(7/2)) + (105*J2*ae^2*u*x*y*z^2)/(2*(x^2 + y^2 + z^2)^(9/2));
% dady(2) ... (day/dy)
% A(5,2) = -(3*J2*ae^2*u*(x^2 + y^2 + z^2)^3 - 15*J2*ae^2*u*z^2*(x^2 + y^2 + z^2)^2 + 315*J3*ae^3*u*y^2*z^3 + 15*J3*ae^3*u*z*(x^2 + y^2 + z^2)^2 - 35*J3*ae^3*u*z^3*(x^2 + y^2 + z^2) - 15*J2*ae^2*u*y^2*(x^2 + y^2 + z^2)^2 - 105*J3*ae^3*u*y^2*z*(x^2 + y^2 + z^2) + 105*J2*ae^2*u*y^2*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2));
% Asym(5,2)
A(5,2) = (3*u*y^2)/(x^2 + y^2 + z^2)^(5/2) - u/(x^2 + y^2 + z^2)^(3/2) + (3*J2*ae^2*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*ae^2*u*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*ae^2*u*y^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (105*J2*ae^2*u*y^2*z^2)/(2*(x^2 + y^2 + z^2)^(9/2));
% dadz(2) ... (day/dz)
% A(5,3) = -(315*J3*ae^3*u*y*z^4 + 15*J3*ae^3*u*y*(x^2 + y^2 + z^2)^2 - 45*J2*ae^2*u*y*z*(x^2 + y^2 + z^2)^2 + 105*J2*ae^2*u*y*z^3*(x^2 + y^2 + z^2) - 210*J3*ae^3*u*y*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2));
% Asym(5,3)
A(5,3) = (3*u*y*z)/(x^2 + y^2 + z^2)^(5/2) - (45*J2*ae^2*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (105*J2*ae^2*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
% dadx(3) ... (daz/dx)
% A(6,1) = -(315*J3*ae^3*u*x*z^4 + 15*J3*ae^3*u*x*(x^2 + y^2 + z^2)^2 - 45*J2*ae^2*u*x*z*(x^2 + y^2 + z^2)^2 + 105*J2*ae^2*u*x*z^3*(x^2 + y^2 + z^2) - 210*J3*ae^3*u*x*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2));
% Asym(6,1)
A(6,1) = (3*u*x*z)/(x^2 + y^2 + z^2)^(5/2) - (45*J2*ae^2*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (105*J2*ae^2*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
% dady(3) ... (daz/dy)
% A(6,2) = -(315*J3*ae^3*u*y*z^4 + 15*J3*ae^3*u*y*(x^2 + y^2 + z^2)^2 - 45*J2*ae^2*u*y*z*(x^2 + y^2 + z^2)^2 + 105*J2*ae^2*u*y*z^3*(x^2 + y^2 + z^2) - 210*J3*ae^3*u*y*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2));
% Asym(6,2)
A(6,2) = (3*u*y*z)/(x^2 + y^2 + z^2)^(5/2) - (45*J2*ae^2*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (105*J2*ae^2*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
% dadz(3) ... (daz/dz)
% A(6,3) = -(9*J2*ae^2*u*(x^2 + y^2 + z^2)^3 + 315*J3*ae^3*u*z^5 - 90*J2*ae^2*u*z^2*(x^2 + y^2 + z^2)^2 + 105*J2*ae^2*u*z^4*(x^2 + y^2 + z^2) + 75*J3*ae^3*u*z*(x^2 + y^2 + z^2)^2 - 350*J3*ae^3*u*z^3*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(11/2));
% Asym(6,3)
A(6,3) = (3*u*z^2)/(x^2 + y^2 + z^2)^(5/2) - u/(x^2 + y^2 + z^2)^(3/2) + (9*J2*ae^2*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (45*J2*ae^2*u*z^2)/(x^2 + y^2 + z^2)^(7/2) + (105*J2*ae^2*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2));
%%% Calculate new STM
stm_dot = A*stm;

%%% Creating X-dot
% Spacecraft velocities
dY(1:3) = [dx; dy; dz];

%%% Using u and J2 terms as dynamics
% -u*x/r^3 + dudx
% dY(4) = -u*x/(norm([x,y,z])^3) + (15*J2*ae^2*u*x*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*ae^3*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*J2*ae^2*u*x)/(2*(x^2 + y^2 + z^2)^(5/2)) + (35*J3*ae^3*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
% EQM(4)
dY(4) = (3*J2*ae^2*u*x)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*x)/(x^2 + y^2 + z^2)^(3/2) - (15*J2*ae^2*u*x*z^2)/(2*(x^2 + y^2 + z^2)^(7/2));
% -u*y/r^3 + dudy
% dY(5) = -u*y/(norm([x,y,z])^3) + (15*J2*ae^2*u*y*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*ae^3*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*J2*ae^2*u*y)/(2*(x^2 + y^2 + z^2)^(5/2)) + (35*J3*ae^3*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2));
% EQM(5)
dY(5) = (3*J2*ae^2*u*y)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*y)/(x^2 + y^2 + z^2)^(3/2) - (15*J2*ae^2*u*y*z^2)/(2*(x^2 + y^2 + z^2)^(7/2));
% -u*z/r^3 + dudz
% dY(6) = -u*z/(norm([x,y,z])^3) + (3*J3*ae^3*u)/(2*(x^2 + y^2 + z^2)^(5/2)) + (15*J2*ae^2*u*z^3)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*ae^3*u*z^2)/(x^2 + y^2 + z^2)^(7/2) + (35*J3*ae^3*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2)) - (9*J2*ae^2*u*z)/(2*(x^2 + y^2 + z^2)^(5/2));
% EQM(6)
dY(6) = (9*J2*ae^2*u*z)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*ae^2*u*z^3)/(2*(x^2 + y^2 + z^2)^(7/2)) - (u*z)/(x^2 + y^2 + z^2)^(3/2);

% Filling in reshaped (6^2,1) STM to state
dY(7:end) = reshape(stm_dot,36,1);

end



