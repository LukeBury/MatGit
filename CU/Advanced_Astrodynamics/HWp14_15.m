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
%%% Run Switches
% ========================================================================


% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% System
% ------------------------------------
u = 1;

% ------------------------------------
%%% Initial conditions
% ------------------------------------
% a0 = 2;
% e0 = 0.1;
% i0_deg = 30; % deg
% raan0_deg = 0;
% w0_deg = 45;
% tau0_deg = 0;
a0 = 1;
e0 = 0.1;
i0_deg = 45; % deg
raan0_deg = 0;
w0_deg = 45;
tau0_deg = 0;
[r0, v0] = COE2RV(a0, e0, i0_deg, raan0_deg, w0_deg, tau0_deg, u)

% r0 = [1; 0 ; 0]; % 
% v0 = [0; 1; 0]; % km/s
X0 = [r0; v0];

% ------------------------------------
%%% Perturbing Force
% ------------------------------------
aPerts = [0.001; 0.01; 0.1; 1];

% ------------------------------------
%%% Initial elements
% ------------------------------------
% [a0,e0,i0,raan0,w0,ta0] = RV2COE(r0,v0,u);
Tp0 = 2*pi*sqrt((a0^3)/u);
meanMot0 = 2*pi/Tp0;
% ------------------------------------
%%% Integrator prep
% ------------------------------------
%%% Time info
t0 = 0;
% tf = 10*Tp0;
tf = 5*Tp0;

n_t = 10000;
time0 = linspace(t0, tf, n_t);

time_orbits = time0./Tp0;

%%% ODE Tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Integration and Analyzing
% ========================================================================
% ------------------------------------
%%% Integrating
% ------------------------------------
%%% Integrating each case
[~, X_1] = ode113(@Int_2BI_constPert, time0, X0, options, u, [0; 0; aPerts(1)]);
[~, X_2] = ode113(@Int_2BI_constPert, time0, X0, options, u, [0; 0; aPerts(2)]);
[~, X_3] = ode113(@Int_2BI_constPert, time0, X0, options, u, [0; 0; aPerts(3)]);
[~, X_4] = ode113(@Int_2BI_constPert, time0, X0, options, u, [0; 0; aPerts(4)]);

% ------------------------------------
%%% Calculating OEs
% ------------------------------------
%%% Preallocating for numerical
as    = zeros(n_t,4);
es    = zeros(n_t,4);
is_deg    = zeros(n_t,4);
raans_deg = zeros(n_t,4);
ws_deg    = zeros(n_t,4);
tas_deg   = zeros(n_t,4);

%%% preallocating for analytical
as_analytical = zeros(n_t,4);
es_analytical = zeros(n_t,4);

delta_a = @(n,e,i,w,as,E) (2*as*sin(i)/(n*n))*(sqrt(1-e^2)*cos(w)*sin(E) - sin(w)*(1-cos(E)));
delta_e = @(n,e,i,w,a,as,E) (as*sqrt(1-e^2)*sin(i)/(n*n*a))*((3/2)*cos(w)*E - 2*e*cos(w)*sin(E) + (1/4)*cos(w)*sin(2*E) - (sqrt(1-e^2)/4)*sin(w)*(1-cos(2*E)));
adjust_e1 = 0;
adjust_e2 = 0;
adjust_e3 = 0;
adjust_e4 = 0;

d2r = pi/180;

E1_rad = 0;
E2_rad = 0;
E3_rad = 0;
E4_rad = 0;

E_augmentation = 0;
for kk = 1:n_t
    %%% Calculating OEs from numerical
    [a1,e1,i1_deg,raan1_deg,w1_deg,ta1_deg] = RV2COE(X_1(kk,1:3),X_1(kk,4:6),u);
    [a2,e2,i2_deg,raan2_deg,w2_deg,ta2_deg] = RV2COE(X_2(kk,1:3),X_2(kk,4:6),u);
    [a3,e3,i3_deg,raan3_deg,w3_deg,ta3_deg] = RV2COE(X_3(kk,1:3),X_3(kk,4:6),u);
    [a4,e4,i4_deg,raan4_deg,w4_deg,ta4_deg] = RV2COE(X_4(kk,1:3),X_4(kk,4:6),u);
    
    if kk == 1
        ta1_deg = tau0_deg;
        ta2_deg = tau0_deg;
        ta3_deg = tau0_deg;
        ta4_deg = tau0_deg;
    end
        
    %%% Assigning old Es
    E1_old = E1_rad;
    E2_old = E2_rad;
    E3_old = E3_rad;
    E4_old = E4_rad;
    
    %%% Getting necessary parameters
    [E1_rad] = T2E(ta1_deg*d2r,e1);
    [E2_rad] = T2E(ta2_deg*d2r,e2);
    [E3_rad] = T2E(ta3_deg*d2r,e3);
    [E4_rad] = T2E(ta4_deg*d2r,e4);
    
    %%% Correcting E so it can go beyond 2pi
    E1_new = E1_rad;
    E2_new = E2_rad;
    E3_new = E3_rad;
    E4_new = E4_rad;
    
    if kk > 2
        if (E1_old - E1_new > 0)
            E_augmentation = ceil((E1_old - E1_new)/(2*pi))*2*pi;
        end
    end
    E1_rad = E1_rad + E_augmentation;
    
    %%% Computing analytical changes
    n1 = sqrt(u/(a1^3));
    n2 = sqrt(u/(a2^3));
    n3 = sqrt(u/(a3^3));
    n4 = sqrt(u/(a4^3));
    
    a1_an = delta_a(n1,e0,i0_deg*d2r,w0_deg*d2r,aPerts(1),E1_rad);
    a2_an = delta_a(n2,e0,i0_deg*d2r,w0_deg*d2r,aPerts(2),E2_rad);
    a3_an = delta_a(n3,e0,i0_deg*d2r,w0_deg*d2r,aPerts(3),E3_rad);
    a4_an = delta_a(n4,e0,i0_deg*d2r,w0_deg*d2r,aPerts(4),E4_rad);
    
    e1_an = delta_e(n1,e0,i0_deg*d2r,w0_deg*d2r,a1,aPerts(1),E1_rad);
    e2_an = delta_e(n2,e0,i0_deg*d2r,w0_deg*d2r,a2,aPerts(2),E1_rad);
    e3_an = delta_e(n3,e0,i0_deg*d2r,w0_deg*d2r,a3,aPerts(3),E1_rad);
    e4_an = delta_e(n4,e0,i0_deg*d2r,w0_deg*d2r,a4,aPerts(4),E1_rad);
    
    as_analytical(kk,:) = [a1_an; a2_an; a3_an; a4_an];
    es_analytical(kk,:) = [e1_an; e2_an; e3_an; e4_an];
    
%     if kk > 1
%         if es_analytical(kk-1,1) > e1_an
%             d
%         end
%     end
%     es_analytical(kk,:) = [e1_an+adjust_e1; e2_an+adjust_e2; e3_an+adjust_e3; e4_an+adjust_e4];
%     
%     
    
    %%% Storing OEs
    as(kk,:)        = [   a1,    a2,    a3,    a4];
    es(kk,:)        = [   e1,    e2,    e3,    e4];
    is_deg(kk,:)    = [   i1_deg,    i2_deg,    i3_deg,    i4_deg];
    raans_deg(kk,:) = [raan1_deg, raan2_deg, raan3_deg, raan4_deg];
    ws_deg(kk,:)    = [   w1_deg,    w2_deg,    w3_deg,    w4_deg];
    tas_deg(kk,:)   = [  ta1_deg,   ta2_deg,   ta3_deg,   ta4_deg];
end




% ========================================================================
%%% Plotting
% ========================================================================
% ------------------------------------
%%% Plotting full systems
% ------------------------------------
figure; hold all
p1 = plot(X_1(:,1),X_1(:,2),'b','linewidth',1);
p2 = plot(0,0,'k.','markersize',30);
PlotBoi2('$X$','$Y$',18,'LaTex')
axis equal
legend([p1 p2],'Perturbed Orbit','Central Body')
title(sprintf('a_s = %1.3f',aPerts(1)))

figure; hold all
p1 = plot(X_2(:,1),X_2(:,2),'b','linewidth',1);
p2 = plot(0,0,'k.','markersize',30);
PlotBoi2('$X$','$Y$',18,'LaTex')
axis equal
legend([p1 p2],'Perturbed Orbit','Central Body')
title(sprintf('a_s = %1.3f',aPerts(2)))

% figure; hold all
% p1 = plot(X_3(:,1),X_3(:,2),'b','linewidth',1);
% p2 = plot(0,0,'k.','markersize',30);
% PlotBoi2('$X$','$Y$',18,'LaTex')
% axis equal
% legend([p1 p2],'Perturbed Orbit','Central Body')
% title(sprintf('a_s = %1.3f',aPerts(3)))
% 
% figure; hold all
% p1 = plot(X_4(:,1),X_4(:,2),'b','linewidth',1);
% p2 = plot(0,0,'k.','markersize',30);
% PlotBoi2('$X$','$Y$',18,'LaTex')
% axis equal
% legend([p1 p2],'Perturbed Orbit','Central Body')
% title(sprintf('a_s = %1.3f',aPerts(4)))
% ------------------------------------
%%% Plotting OEs
% ------------------------------------
figure('position',[307 303 888 420])
subplot(2,2,1); hold all
title(sprintf('a_s = %1.3f',aPerts(1)))
plot(time_orbits,as(:,1),'b','linewidth',2);
plot(time_orbits,as_analytical(:,1)+a0,'r','linewidth',2);
PlotBoi2('','$a$',18,'LaTex')

subplot(2,2,3); hold all
plot(time_orbits,as_analytical(:,1)+a0 - as(:,1),'m','linewidth',2);
PlotBoi2('Tp','$a$ discrepancy',18,'LaTex')

subplot(2,2,2); hold all
p1 = plot(time_orbits,es(:,1),'b','linewidth',2);
p2 = plot(time_orbits,es_analytical(:,1)+e0,'r','linewidth',2);
legend([p1, p2],'Numerical','Analytical')
PlotBoi2('','$e$',18,'LaTex')

subplot(2,2,4); hold all
plot(time_orbits,es_analytical(:,1)+e0 - es(:,1),'m','linewidth',2);
PlotBoi2('Tp','$e$ discrepancy',18,'LaTex')




figure('position',[486 182 888 420])
subplot(2,2,1); hold all
title(sprintf('a_s = %1.3f',aPerts(2)))
plot(time_orbits,as(:,2),'b','linewidth',2);
plot(time_orbits,as_analytical(:,2)+a0,'r','linewidth',2);
PlotBoi2('','a',18,'LaTex')

subplot(2,2,3); hold all
plot(time_orbits,as_analytical(:,2)+a0 - as(:,2),'m','linewidth',2);
PlotBoi2('Tp','$a$ discrepancy',18,'LaTex')

subplot(2,2,2); hold all
p1 = plot(time_orbits,es(:,2),'b','linewidth',2);
p2 = plot(time_orbits,es_analytical(:,2)+e0,'r','linewidth',2);
legend([p1, p2],'Numerical','Analytical')
PlotBoi2('Tp','e',18,'LaTex')

subplot(2,2,4); hold all
plot(time_orbits,es_analytical(:,2)+e0 - es(:,2),'m','linewidth',2);
PlotBoi2('Tp','$e$ discrepancy',18,'LaTex')





figure('position',[486 182 888 420])
subplot(2,2,1); hold all
plot(time_orbits, is_deg(:,1),'b','linewidth',2)
title(sprintf('a_s = %1.3f',aPerts(1)))
PlotBoi2('','$i$, $^{\circ}$',18,'LaTex')

subplot(2,2,3); hold all
plot(time_orbits, raans_deg(:,1),'b.','linewidth',2)
PlotBoi2('Tp','$\Omega$, $^{\circ}$',18,'LaTex')

subplot(2,2,2); hold all
plot(time_orbits, ws_deg(:,1),'b','linewidth',2)
PlotBoi2('','$\omega$, $^{\circ}$',18,'LaTex')

subplot(2,2,4); hold all
plot(time_orbits, tas_deg(:,1),'b.','linewidth',2)
PlotBoi2('Tp','$f$, $^{\circ}$',18,'LaTex')



figure('position',[486 182 888 420])
subplot(2,2,1); hold all
plot(time_orbits, is_deg(:,2),'b','linewidth',2)
title(sprintf('a_s = %1.3f',aPerts(2)))
PlotBoi2('','$i$, $^{\circ}$',18,'LaTex')

subplot(2,2,3); hold all
plot(time_orbits, raans_deg(:,2),'b.','linewidth',2)
PlotBoi2('Tp','$\Omega$, $^{\circ}$',18,'LaTex')

subplot(2,2,2); hold all
plot(time_orbits, ws_deg(:,2),'b','linewidth',2)
PlotBoi2('','$\omega$, $^{\circ}$',18,'LaTex')

subplot(2,2,4); hold all
plot(time_orbits, tas_deg(:,2),'b.','linewidth',2)
PlotBoi2('Tp','$f$, $^{\circ}$',18,'LaTex')


% figure
% subplot(2,1,1); hold all
% title(sprintf('a_s = %1.3f',aPerts(3)))
% plot(time_orbits,as(:,3),'b','linewidth',2)
% plot(time_orbits,as_analytical(:,3),'r','linewidth',2)
% PlotBoi2('','a',18,'LaTex')
% subplot(2,1,2); hold all
% plot(time_orbits,es(:,3),'b','linewidth',2)
% plot(time_orbits,es_analytical(:,3),'r','linewidth',2)
% PlotBoi2('Tp','e',18,'LaTex')
% 
% 
% figure
% subplot(2,1,1); hold all
% title(sprintf('a_s = %1.3f',aPerts(4)))
% plot(time_orbits,as(:,4),'b','linewidth',2)
% plot(time_orbits,as_analytical(:,4),'r','linewidth',2)
% PlotBoi2('','a',18,'LaTex')
% subplot(2,1,2); hold all
% plot(time_orbits,es(:,4),'b','linewidth',2)
% plot(time_orbits,es_analytical(:,4),'r','linewidth',2)
% PlotBoi2('Tp','e',18,'LaTex')




% ========================================================================
%%% 15-f
% ========================================================================

dabar1 = -aPerts(1)*sind(i0_deg)*sind(w0_deg)*(e0+2) / (meanMot0^2)
dabar2 = -aPerts(2)*sind(i0_deg)*sind(w0_deg)*(e0+2) / (meanMot0^2)
mean1_an = mean(as_analytical(:,1))
mean2_an = mean(as_analytical(:,2))
mean1_num = mean(as(:,1)) - a0
mean2_num = mean(as(:,2)) - a0





