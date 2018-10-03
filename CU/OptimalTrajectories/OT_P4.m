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
%%% Defining local-circular and bi-elliptic velocities (normalized)
r1 = 1;
r2 = 1;
r3 = 10;
r = r2/r1;
l = r3/r1;

vLc = 1;
vBi_p = sqrt(2*l/(1+l));
vBi_a = sqrt(2/(l*(1+l)));

% Creating x-y grid
resolution = 10000;
di_1 = linspace(0,90,resolution);
di_3 = di_1;
[di1, di3] = meshgrid(di_1,di_3);

%%% Calculating dvs across di1-di3 grid
dv_total = zeros(size(di1));
for k_di1 = 1:size(di1,1)
    for k_di3 = 1:size(di1,2)
        %%% Assign various di values
        di1_k = di1(k_di1,k_di3);
        di3_k = di3(k_di1,k_di3);
        di2_k = 90 - di1_k - di3_k;
        
        if di2_k < 0
            dv_total(k_di1,k_di3) = NaN;
        else
            %%% Calculate dvs
            % first burn at periapsis
            dv1 = sqrt(vLc^2 + vBi_p^2 - 2 * vLc * vBi_p * cosd(di1_k));

            % burn at apoapsis
            dv2 = 2 * vBi_a * sind(di2_k/2);

            % final burn at periapsis
            dv3 = sqrt(vLc^2 + vBi_p^2 - 2 * vLc * vBi_p * cosd(di3_k));

            % total dv
            dv_total(k_di1,k_di3) = dv1 + dv2 + dv3;
        end
        
    end
end

min_dv = inf;
for kk = 1:size(dv_total,1)
  for jj = 1:size(dv_total,2)
      if dv_total(kk, jj) < min_dv
          min_dv = dv_total(kk, jj);
          minRow = kk;
          minCol = jj;
      end
  end
end

min_di1 = di1(minRow, minCol);
min_di3 = di3(minRow, minCol);
min_di2 = 90 - min_di3 - min_di1;

min_dv
min_di1
min_di2
min_di3

if resolution <= 1000
    surf(di_1,di_3,dv_total)
end

toc



















