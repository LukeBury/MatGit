clear
clc
close all
% ------------------------------------------------------------------------
%%% 3a
% ------------------------------------------------------------------------
syms Om dt s
A = [0 1 0 0; 0 0 0 -Om; 0 0 0 1;0 Om 0 0]

% ------------------------------------------------------------------------
%%% 3bi
% ------------------------------------------------------------------------
fprintf('---------------------------- 3bi ----------------------------\n')
clear Om
dt = 0.5; % sec
Om = 0.045; % rad/s
Omdt = Om*dt; 

F = [1 sin(Omdt)/Om 0 -(1-cos(Omdt))/Om;...
    0 cos(Omdt) 0 -sin(Omdt);...
    0 (1-cos(Omdt))/Om 1 sin(Omdt)/Om;...
    0 sin(Omdt) 0 cos(Omdt)]

% ------------------------------------------------------------------------
%%% 3bii
% ------------------------------------------------------------------------
fprintf('---------------------------- 3bii ----------------------------\n')
%%% Givens
u0 = [0, 85*cos(pi/4), 0, -85*sin(pi/4)]'; % m m/s m m/s
P0 = diag([10, 2, 10, 2],0); % m^2 (m/s)^2 m^2 (m/s)^2

%%% Preparing for iterations
ti = 1;
tf = 101;
time = [ti:tf];
uk = zeros(4,101);
Pk = zeros(4,4,101);
% 2-sigma magnitude vector
sig2 = zeros(4,101);

%%% Iterating
for i = ti:tf
    %%% Updating uk & Pk
    uk(:,i) = F^(i)*u0;
    Pk(:,:,i) = (F^(i)) * P0 * (F^(i))';
    sig2(1,i) = 2*sqrt(Pk(1,1,i));
    sig2(2,i) = 2*sqrt(Pk(2,2,i));
    sig2(3,i) = 2*sqrt(Pk(3,3,i));
    sig2(4,i) = 2*sqrt(Pk(4,4,i));
end

fontsize = 20;
figure
hold all
%%% Plot East position and 2-sigma
p1 = plot(time*dt-.5,uk(1,:),'b','linewidth',2);
e1 = plot(time*dt-.5,uk(1,:)+sig2(1,:),'--g','linewidth',2);
plot(time*dt-.5,uk(1,:)-sig2(1,:),'--g','linewidth',2)
%%% Plot North position and 2-sigma
p2 = plot(time*dt-.5,uk(3,:),'r','linewidth',2);
plot(time*dt-.5,uk(3,:)+sig2(3,:),'--g','linewidth',2)
plot(time*dt-.5,uk(3,:)-sig2(3,:),'--g','linewidth',2)
%%% Formatting plot
PlotBoi2('Time, sec','Position, m',fontsize)
legend([p1 p2 e1],'East Position','North Position','2-Sigma')

figure
hold all
%%% Plot East velocity and 2-sigma
p1 = plot(time*dt-.5,uk(2,:),'b','linewidth',2);
e1 = plot(time*dt-.5,uk(2,:)+sig2(2,:),'--g','linewidth',2);
plot(time*dt-.5,uk(2,:)-sig2(2,:),'--g','linewidth',2)
%%% Plot North velocity and 2-sigma
p2 = plot(time*dt-.5,uk(4,:),'r','linewidth',2);
plot(time*dt-.5,uk(4,:)+sig2(4,:),'--g','linewidth',2)
plot(time*dt-.5,uk(4,:)-sig2(4,:),'--g','linewidth',2)
PlotBoi2('Time, sec','Velocity, m/s',fontsize)
legend([p1 p2 e1],'East Velocity','North Velocity','2-Sigma')

% ------------------------------------------------------------------------
%%% 3ci
% ------------------------------------------------------------------------
fprintf('---------------------------- 3ci ----------------------------\n')
dt = 0.5; % sec
Oma = 0.045; % rad/s
Omb = -0.045; % rad/s
Omat = dt * Oma; % rad
Ombt = dt * Omb; % rad

Fa = [1 sin(Omat)/Oma 0 -(1-cos(Omat))/Oma;...
    0 cos(Omat) 0 -sin(Omat);...
    0 (1-cos(Omat))/Oma 1 sin(Omat)/Oma;...
    0 sin(Omat) 0 cos(Omat)];

Fb = [1 sin(Ombt)/Omb 0 -(1-cos(Ombt))/Omb;...
    0 cos(Ombt) 0 -sin(Ombt);...
    0 (1-cos(Ombt))/Omb 1 sin(Ombt)/Omb;...
    0 sin(Ombt) 0 cos(Ombt)];
fprintf('\n')
% ------------------------------------------------------------------------
%%% 3ciii
% ------------------------------------------------------------------------
fprintf('---------------------------- 3ciii ----------------------------\n')
%%% Givens
ua0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)]; % m m/s m m/s
Pa0 = diag([10, 4, 10, 4],0); % m^2 (m/s)^2 m^2 (m/s)^2
ub0 = [3200; 85*cos(pi/4); 3200; -85*sin(pi/4)]; % m m/s m m/s
Pb0 = diag([11, 3.5, 11, 3.5],0); % m^2 (m/s)^2 m^2 (m/s)^2

%%% Preparing for iterations
ti = 1;
tf = 201;
time = [ti:tf];

uak = zeros(4,tf);
Pak = zeros(4,4,tf);

ubk = zeros(4,tf);
Pbk = zeros(4,4,tf);

uck = zeros(4,tf);
Pck = zeros(4,4,tf);


%%% For resizing
uck2 = zeros(2,tf);
Pck2 = zeros(2,2,tf);

uaktest = zeros(4,2);
ubktest = zeros(4,2);
Paktest = zeros(4,4,2)
Pbktest =zeros(4,4,2)
ucktest = zeros(2,2);
Pcktest = zeros(2,2,2);

ZR = 100; % m
ER = 100; % m

mvns = zeros(tf,1);
%%% Iterating
for i = ti:tf
    %%% Updating uk & Pk
    uak(:,i) = Fa^(i)*ua0;
    Pak(:,:,i) = (Fa^(i)) * Pa0 * (Fa^(i))';
    
    ubk(:,i) = Fb^(i)*ub0;
    Pbk(:,:,i) = (Fb^(i)) * Pb0 * (Fb^(i))';
    
    uck(:,i) = uak(:,i) - ubk(:,i);
    Pck(:,:,i) = Pak(:,:,i) + Pbk(:,:,i);
    
    %%% Resizing
    uck2(:,i) = [uck(1,i);uck(3,i)];
    Pck2(:,:,i) = [Pck(1,1,i), Pck(1,3,i); Pck(3,1,i), Pck(3,3,i)];
        
    %%% Calculating probablity of collision
    mvns(i,1) = mvncdf(-[ZR;ER],[ZR;ER],uck2(:,i),Pck2(:,:,i));
end

figure
hold all
plot(uak(1,:),uak(3,:),'linewidth',2)
plot(ubk(1,:),ubk(3,:),'linewidth',2)
legend('Plane A','Plane B')
PlotBoi2('East-position, m','West-position, m',fontsize)

figure
plot(time*dt-.5,mvns,'linewidth',2)
PlotBoi2('Time,sec','Probability of collision',fontsize)




