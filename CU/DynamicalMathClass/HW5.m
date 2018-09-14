clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
colors = get_colors();

run_p4 = 0;
% ========================================================================
%%% P4
% ========================================================================
if run_p4 == 1
% ----------------
%%% a,b,d - computing exponent, plotting attractor, destroying attractor
% ----------------
%%% Setting time vector
t0 = 1;
dt = 1;
tf = 100;
time = t0:dt:tf;

%%% Initial Conditions (given) % (a) and (d) are reduced parameters
x0 = [.7;0;-.5];
v0 = [0;-.1;.1];
w0 = [x0;v0];
b = 1; c = 1; e = 1; f = 1;
a = 0.4;
d = -1;

%%% Initial Conditions (DESTROY) % (a) and (d) are reduced parameters
% x0 = [.7;0;-.5];
% v0 = [0;-.1;.1];
% w0 = [x0;v0];
% b = 1; c = 1; e = 1; f = 1;
% a = 3;
% d = -1;

%%% Computing the scale factor
scalefactor = -log(norm(v0));

%%% Preallocating
T=[]; L=[];states = [];
%%% Integrating
for kk = 1:length(time)
    ti = time(kk);
    
    tol = 1e-9;
    options = odeset('RelTol',tol,'AbsTol',tol);
    [t,state] = ode45(@int_appm5460_hw5,[ti,ti+dt],w0,options,a,b,c,d,e,f);
    states = [states;state];
    
    Lyp = (scalefactor + 0.5*log(sum(state(:,4:6).^2,2)))./t;
    T = [T; t];  % Store integration output
    L=[L; Lyp];
    nm = norm(state(end,4:6));
    scalefactor = scalefactor + log(nm);
    w0 = horzcat(state(end,1:3), state(end,4:6)/nm);
end
plot(T,L,'.','color',colors.std.purp);
PlotBoi2('Time','Lyapunov Exponent',12)
grid on
maxLyapunovExponent = Lyp(end)

figure
plot3(states(:,1),states(:,2),states(:,3),'.','color',colors.std.purp)
PlotBoi3('x','y','z',12)
grid on
title(sprintf('a = %1.2f     d = %1.2f',a,d))


% ----------------
%%% c - Computing the box counting dimension
% ----------------
%%% Finding bounds for box counting
warning('May not work perfectly...')
xmin = min(states(:,1)); xmax = max(states(:,1));
ymin = min(states(:,2)); ymax = max(states(:,2));
zmin = min(states(:,3)); zmax = max(states(:,3));

xrange = xmax-xmin;
yrange = ymax-ymin;
zrange = zmax-zmin;
maxRange = max([xrange,yrange,zrange]); % max range between x, y, and z

side = maxRange;
sides = [maxRange, maxRange/2, maxRange/4, maxRange/8, maxRange/16, maxRange/32, maxRange/64, maxRange/128]; % side lengths for cubes (boxes)
boxCounts = zeros(length(sides),2);
tic
for s_i = 1:length(sides) % for each box length
    side = sides(s_i); % set box length
    boxCount = 0; % set box count back to zero
    xp = xmin;
    while xp <= (xmax) % set x-range for current box
        xm = xp;
        xp = xm + side;
        yp = ymin;
        while yp <= (ymax) % set y-range for current box
            ym = yp; 
            yp = ym + side;
            zp = zmin;
            while zp <= (zmax) % set z-range for current box
                zm = zp;
                zp = zm + side;
                boxCounted = 0; % reset whether this box includes the attractor
                for kk = 1:size(states,1) % for every state
                    if (xm <= states(kk,1) <= xp) && (ym <= states(kk,2) <= yp) && (zm <= states(kk,3) <= zp) && (boxCounted == 0) % is this state in the current box?
                        boxCount = boxCount + 1; % count this box towards the total for this box size
                        boxCounted = 1; % this box includes the attractor
                        break % no need to look through further states
                    end
                end
            end
        end
    end
    boxCounts(s_i,1) = side; % side length of current box search
    boxCounts(s_i,2) = boxCount; % how many boxes captured the attractor
    s_i
    toc
end
plot_x = log(1./boxCounts(:,1));
plot_y = log(boxCounts(:,2));
linearFit = polyfit(plot_x,plot_y,1);
boxCountingDimension = linearFit(1)

figure; hold all
plot(plot_x,plot_y,'linewidth',1.5)
x = [min(plot_x):.1:max(plot_x)];
y = linearFit(1).*x + linearFit(2);
plot(x,y,'--','linewidth',1.5)
legend('Results','Linear fit')
grid on
PlotBoi2('log(1/d)','log(N(d))',14)

distFig('Screen','Main','Scale',0.8) % 'Main' or 'Secondary'

end % run_p4
% ========================================================================
%%% P5
% ========================================================================
% ----------------
%%% b
% ----------------
% u = -3/(2^(2/3))
% x = -3:.1:3;
% xd = 1 + u*x + x.^3;
% 
% figure
% plot(x,xd)
% grid on
% ylim([-3 3])

% ----------------
%%% c
% ----------------
% u1 = 0
% u2 = -1
% u3 = -.5
% u4 = -2
% u5 = 2
% x = -3:.001:3;
% xd1 = u1*x + 2*x.^2 - x.^3;
% xd2 = u2*x + 2*x.^2 - x.^3;
% xd3 = u3*x + 2*x.^2 - x.^3;
% xd4 = u4*x + 2*x.^2 - x.^3;
% xd5 = u5*x + 2*x.^2 - x.^3;
% 
% figure; hold all
% plot(x,xd1,'linewidth',1.5)
% plot(x,xd2,'linewidth',1.5)
% plot(x,xd3,'linewidth',1.5)
% plot(x,xd4,'linewidth',1.5)
% plot(x,xd5,'linewidth',1.5)
% legend('u1','u2','u3','u4','u5')
% plot([-3,3],[0,0],'k')
% plot([0,0],[-3,3],'k')
% grid on
% ylim([-3 3])

% ----------------
%%% d
% ----------------
% xLimits = 10;
% u1 = .2;
% u2 = -.15;
% x = -xLimits:.001:xLimits;
% xd1 = u1*x + sin(x);
% xd2 = u2*x + sin(x);
% 
% figure; hold all
% plot(x,xd1,'.','linewidth',1.5)
% plot(x,xd2,'.','linewidth',1.5)
% plot([-xLimits,xLimits],[0,0],'k')
% plot([0,0],[-xLimits,xLimits],'k')
% legend('u1','u2')
% 
% ylim([min(x) max(x)])
% 
% 
% xLimits = 5;
% x = -xLimits:.00001:xLimits;
% figure; hold all
% y = x - tan(x);
% plot(x,y,'.','linewidth',1.5)
% plot([-xLimits,xLimits],[0,0],'k')
% plot([0,0],[-xLimits,xLimits],'k')
% ylim([min(x) max(x)])


% ========================================================================
%%% P7
% ========================================================================
syms x y lam real
%%% Finding equilibria
eqns = [x*lam + 2*x*y - x^2 == 0,(lam - 1)*y + x^2 == 0];
eq_pts = solve(eqns, [x y]);
pretty(eq_pts.x)
pretty(eq_pts.y)

%%% Jacobian
A = [lam+2*y-2*x, 2*x;...
     2*x, lam-1]

%%% eq 1
x=0, y=0;
lam = 0;
A_1_1 = subs(A)
[v1_1,e1_1] = eig(A_1_1);
e1_1

lam = 1/9;
A_1_2 = subs(A)
[v1_2,e1_2] = eig(A_1_2);
e1_2

lam = 1;
A_1_3 = subs(A)
[v1_3,e1_3] = eig(A_1_3);
e1_3

%%% eq 2
x = 1/4 + lam/4 + sqrt((9*lam-1)*(lam-1)/16);
y = 1/8 - 5*lam/8 + sqrt((9*lam-1)*(lam-1)/16);
lam = 0;
A_2_1 = subs(A);
[v2_1,e2_1] = eig(A_2_1);
e2_1

lam = 1/9;
A_2_2 = subs(A);
[v2_2,e2_2] = eig(A_2_2);
e2_2

lam = 1;
A_2_3 = subs(A);
[v2_3,e2_3] = eig(A_2_3);
e2_3

%%% eq 3
x = 1/4 + lam/4 - sqrt((9*lam-1)*(lam-1)/16);
y = 1/8 - 5*lam/8 - sqrt((9*lam-1)*(lam-1)/16);
lam = 0;
A_3_1 = subs(A);
[v3_1,e3_1] = eig(A_3_1);
e3_1

lam = 1/9;
A_3_2 = subs(A);
[v3_2,e3_2] = eig(A_3_2);
e3_2

lam = 1;
A_3_3 = subs(A);
[v3_3,e3_3] = eig(A_3_3);
e3_3

















function [dX] = int_appm5460_hw5(t,X,a,b,c,d,e,f)
%%% Preallocate state output
dX = zeros(6,1);

%%% My Equations of Motion
v1 = a*X(1) + b*X(3);
v2 = c*X(1)*X(3) + d*X(2);
v3 = -e*X(1) + f*X(2);
vd1 = a*X(4)+b*X(6);
vd2 = c*(X(4)*X(3) + X(1)*X(6)) + d*X(5);
vd3 = -e*X(4) + f*X(5);

%%% Output the derivative of the state
dX(1:3) = [v1;v2;v3];
dX(4:6) = [vd1;vd2;vd3];
end
