function plotCR3BP_YZNeck( JC, u , Lpoint, on_J21, prms, color, linewidth)
%%% Returns points of certain level of a given contour map
%%% Inputs:
%       1) JC - Corresponding Jacobi Constant
%       2) u  - Mass ration of CR3BP system
%       3) Lpoint - Lagrange point (which neck - 1 or 2)
%       4) on_J21 - account for J2 of primary? (1-yes, 0-no)
%       5) prms - structure, fields: J21, R1_n, R2_n
%       6) color - color for plot [1x3]
%       7) linewidth - linewidth for plot
%
%%% Outputs: 
%       1) 
%=========================================================================
% -------------------------------------------------
% Acquire Collinear Lagrange points
% -------------------------------------------------
if on_J21 == 0
    L123 = EquilibriumPoints(u,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
elseif on_J21 == 1
    L123 = EquilibriumPoints_J2(u,prms.J21,0,prms.R1_n,prms.R2_n,1:3);
end

% -------------------------------------------------
% Finding upper-y-value of neck at L-point
% -------------------------------------------------
%%% Set values
x = L123(Lpoint,1);
z = 0;

%%% JC function equal to zero
r1 = @(x,y,z) sqrt((x+u)^2+y.^2+z^2);
r2 = @(x,y,z) sqrt((x-1+u)^2+y.^2+z^2);
if on_J21 == 0
    f = @(y) x^2 + y.^2 - JC + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z);
elseif on_J21 == 1
    f = @(y) x^2 + y.^2 - JC + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z) + prms.R1_n*prms.R1_n*prms.J21*(1-u)*(r1(x,y,z)^2 - 3*z*z)/(r1(x,y,z)^5);
end


%%% Find the root of this function in the appropriate range
y_neck_upper = fzero(f,[0 10*prms.R2_n]);

%%% Clear variables
clear x f z


% -------------------------------------------------
% Finding upper-z-value of neck at L-point
% -------------------------------------------------
%%% Set values
x = L123(Lpoint,1);
y = 0;

%%% JC function equal to zero
r1 = @(x,y,z) sqrt((x+u)^2+y.^2+z^2);
r2 = @(x,y,z) sqrt((x-1+u)^2+y.^2+z^2);

if on_J21 == 0
    f = @(z) x^2 + y.^2 - JC + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z);
elseif on_J21 == 1
    f = @(z) x^2 + y.^2 - JC + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z) + prms.R1_n*prms.R1_n*prms.J21*(1-u)*(r1(x,y,z)^2 - 3*z*z)/(r1(x,y,z)^5);
end

%%% Find the root of this function in the appropriate range
z_neck_upper = fzero(f,[0 10*prms.R2_n]);

%%% Clear variables
clear x f y

% -------------------------------------------------
% Create grid of starting locations based on y-z neck to find contour
% points
% -------------------------------------------------
ys = linspace(-2*y_neck_upper, 2*y_neck_upper, 400);
zs = linspace(-2*z_neck_upper, 2*z_neck_upper, 400);
[Y_yz,Z_yz] = meshgrid(ys,zs);

% -------------------------------------------------
% Only keep starting positions that are valid for the energy level
% -------------------------------------------------
%%% Calculating JCs across y-z grid
JCs_yz_Lpoint = zeros(size(Y_yz));
for yk = 1:size(Y_yz,1)
    for zk = 1:size(Y_yz,2)
        %%% Zero-Velocity Curve
        if on_J21 == 0
            zv = JacobiConstantCalculator(u,[L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)] ,[0, 0, 0]);
        elseif on_J21 == 1
            zv = JacobiConstantCalculator_J2(u,[L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)],[0,0,0], prms.R1_n, prms.R2_n, prms.J21, 0);
        end
        JCs_yz_Lpoint(yk,zk) = zv;
    end
end

%%% Get points of y-z contour in 3D space
[ yzContourPoints ] = getContourPoints( Y_yz, Z_yz, JCs_yz_Lpoint, JC );

% -------------------------------------------------
% Plotting
% -------------------------------------------------
plot3(ones(size(yzContourPoints,2),1).*L123(Lpoint,1), yzContourPoints(1,:),yzContourPoints(2,:),'linewidth',linewidth,'color',color)

end

