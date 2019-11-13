clear, clc
close all

th = linspace(0, 2*pi, 30)';

% 3D lines
r1 = [ones(30, 1), sin(th), 1.2 * cos(th)];
r2 = [ones(30, 1) * 2, sin(th) * 0.8, cos(th) * 0.6];
r3 = [ones(30, 1) * 3, sin(th) + 0.5, cos(th)];

% Patches between lines
[Xp1, Yp1, Zp1] = patchLines(r1, r2);
[Xp2, Yp2, Zp2] = patchLines(r2, r3);

% Combines patches
Xp = [Xp1, Xp2];
Yp = [Yp1, Yp2];
Zp = [Zp1, Zp2];

% Plot total object
figure
patch(Xp, Yp, Zp,'b')
grid on
axis equal
view(3)