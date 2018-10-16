clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Running
% ========================================================================
u = 1

%%% 1-impulse escape
dVt_1 = @(r, l, vInf) sqrt(2*u/r + vInf^2) - sqrt(u/r);

%%% 2-impulse escape, l > r (increasing orbit)
dV1_2i = @(r, l, vInf) sqrt(u*(2/r - 2/(l+r))) - sqrt(u/r);
dV2_2i = @(r, l, vInf) sqrt(2*u/l + vInf^2) - sqrt(u*(2/l - 2/(l+r)));

dVt_2i = @(r, l, vInf) dV1_2i(r, l, vInf) + dV2_2i(r, l, vInf);

%%% 2-impulse escape, l < r (decreasing orbit)
dV1_2d = @(r, l, vInf) sqrt(u/r)   - sqrt(u*(2/r - 2/(l+r)));
dV2_2d = @(r, l, vInf) sqrt(2*u/l + vInf^2) - sqrt(u*(2/l - 2/(l+r)));

dVt_2d = @(r, l, vInf) dV1_2d(r, l, vInf) + dV2_2d(r, l, vInf);


% Creating x-y grid
resolution = 100;
rs = linspace(0.01,20,resolution);
vInfs = linspace(0,5,resolution);
[rs, vInfs] = meshgrid(rs,vInfs);

u = 1
l = 25

%%% Calculating dvs across di1-di3 grid
dv_total = zeros(size(rs));
for row = 1:size(rs,1)
    for col = 1:size(rs,2)
        r = rs(row,col);
        vInf = vInfs(row,col);
        
        a = sqrt(2*u/r + vInf^2) - 
        
        
    end
end







