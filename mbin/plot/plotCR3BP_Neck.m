function plotCR3BP_Neck(prms,L123,JC_value,xRes,yRes,contColor,lw)    
%%% Creating ZV contours
% Contour bounds
xCont_min = L123(1,1)-2*prms.R2_n;
xCont_max = L123(2,1)+2*prms.R2_n;
yCont_min = -9*prms.R2_n;
yCont_max = 9*prms.R2_n;
% yCont_min = -5.5*secondary.R_n;
% yCont_max = 5.5*secondary.R_n;

% Creating x-y grid
xs = linspace(xCont_min,xCont_max,xRes); % 600
ys = linspace(yCont_min,yCont_max,yRes); % 200
[X_xy, Y_xy] = meshgrid(xs,ys);
clear xs ys zs

% Calculating JCs across x-y grid
JCs_xy = zeros(size(X_xy));
for xk = 1:size(X_xy,1)
    for yk = 1:size(X_xy,2)
        %%% Zero-Velocity Curve
%         zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0] ,[0, 0, 0]);
        zv = getJacobiConstant_ZH([X_xy(xk,yk), Y_xy(xk,yk), 0, 0, 0, 0],prms);
        JCs_xy(xk,yk) = zv;
    end
end

[xyContourPoints,href] = contour(X_xy,Y_xy,JCs_xy,[JC_value, JC_value],...
    'color',contColor,'linewidth',lw);
end