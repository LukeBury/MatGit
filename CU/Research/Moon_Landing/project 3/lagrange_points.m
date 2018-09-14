clc; clear; close all

tic

M1 = 5.972e24;  % kg (Earth)
M2 = 7.35e22;  % kg (Moon)

G = 6.673e-20;
R = 384400;

% Mass Ratio
v = M2/(M1 + M2);

% Nondimensional masses
m1 = 1 - v;
m2 = v;

r1 = -v;
r2 = 1 - v;

x = linspace(-2, 2, 1000);
y = linspace(-2, 2, 1000);

levels = [3 3.02 3.05 3.1 3.18 linspace(3.3,6,12)];

for i = 1:length(x)
    for j = 1:length(y)
        R1 = sqrt((x(i) + v)^2 + y(j)^2);
        R2 = sqrt((x(i) -1 + v)^2 + y(j)^2);
        
        Omega = (x(i)^2 + y(j)^2)/2 + m1/R1 + m2/R2;
        
        C(j,i) = 2*Omega;
    end
end

toc

figure
box on; grid on; hold on
contour(y, x, C, levels)
h1 = plot(r1, 0, 'ob'); % Earth
h2 = plot(r2, 0, 'ok'); % Moon
h3 = plot(1 - (v/3)^(1/3), 0, 'or'); % L1
h4 = plot(1 + (v/3)^(1/3), 0, 'or'); % L2
h5 = plot(-1 - 5*v/12, 0, 'or'); % L3
h6 = plot((M1-M2)/(M1+M2)/2, sqrt(3)/2, 'or'); % L4
h7 = plot((M1-M2)/(M1+M2)/2, -sqrt(3)/2, 'or'); % L5
set(h1,'MarkerEdgeColor','k','MarkerFaceColor','b','markersize',10)
set(h2,'MarkerEdgeColor','k','MarkerFaceColor','k','markersize',5)
set(h3,'MarkerEdgeColor','k','MarkerFaceColor','r','markersize',5)
set(h4,'MarkerEdgeColor','k','MarkerFaceColor','r','markersize',5)
set(h5,'MarkerEdgeColor','k','MarkerFaceColor','r','markersize',5)
set(h6,'MarkerEdgeColor','k','MarkerFaceColor','r','markersize',5)
set(h7,'MarkerEdgeColor','k','MarkerFaceColor','r','markersize',5)
colormap(jet)
contourcbar
title('Earth-Moon System: Non-Dimensional Lagrange Points')
xlabel('x [-]')
ylabel('y [-]')
axis('square')