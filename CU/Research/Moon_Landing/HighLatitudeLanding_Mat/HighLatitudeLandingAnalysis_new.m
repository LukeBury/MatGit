% ========================================================================
%%% Description
% ========================================================================
% For analyzing the high-latitude-landing-problem data generated in Julia

% Created: 03/20/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
ticWhole = tic;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

%%% Periodic orbit ICs
PO_ICs = get_PO_ICs();

% ========================================================================
%%% Run Switches
% ========================================================================
loadPlot_1yrData    = false;

plot_allLowTrajs    = false;
plot_allHighTrajs   = false;
plot_allTrajs       = false;
plot_periapseMap    = false;
plot_apoapseMap     = false;
plot_RPeriapse      = false;
plot_RApoapse       = false;
plot_xyCrossing     = false;
plot_xyCrossing_3D  = false;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Load the data generated in Julia
% -------------------------------------------------
%%% Data Path
dataPath = '/Volumes/LB_External_Drive/Research/High_Latitude_Landing_Study/Data/Full_Results/';

%%% Columb specifiers
c_trajID    = 1;
c_x         = 2;
c_y         = 3;
c_z         = 4;
c_xd        = 5;
c_yd        = 6;
c_zd        = 7;
c_t         = 8;
c_JC        = 9;
c_R         = 10;
c_a         = 11;
c_e         = 12;
c_i         = 13;
c_E         = 14;
c_H         = 15;
c_hx        = 16;
c_hy        = 17;
c_hz        = 18;
c_lat       = 19;
c_lon       = 20;
c_isPer     = 21;
c_isApo     = 22;
c_isxyCross = 23;


%%% Data files
if ~loadPlot_1yrData
%     file_lowImpactTrajData          = [dataPath, 'lowTrajs_150mps_neckFull_L2_250km_22v0s_4pi.txt'];
%     file_highImpactTrajData          = [dataPath, 'lowTrajs_150mps_neckFull_L1_250km_22v0s_4pi.txt'];
%     file_lowImpactTrajData          = [dataPath, 'lowTrajs_100mps_neckFull_L2_250km_22v0s_4pi.txt'];
%     file_highImpactTrajData          = [dataPath, 'lowTrajs_100mps_neckFull_L1_250km_22v0s_4pi.txt'];
    file_lowImpactTrajData          = [dataPath, 'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi.txt'];
%     file_highImpactTrajData          = [dataPath, 'lowTrajs_50mps_neckFull_L1_250km_22v0s_4pi.txt'];
%     file_lowImpactTrajData          = [dataPath, 'lowTrajs_1mps_neckFull_L2_2km_22v0s_4pi.txt'];
%     file_highImpactTrajData          = [dataPath, 'lowTrajs_1mps_neckFull_L1_250km_22v0s_4pi.txt'];

    % file_lowImpactTrajData          = [dataPath, 'lowImpactTrajs_150mps_4pi.txt'];
    % file_lowImpactTrajData          = [dataPath, 'lowImpactTrajsnew_L1_150.txt'];
    % file_lowImpactTrajData          = [dataPath, 'lowImpactTrajsnew_L1_fullNeck_150.txt'];

    % file_lowImpactTrajData          = [dataPath, 'lowTrajs_1mps_neckFull_L1_2km_22v0s_1yr.txt'];
    % file_lowImpactTrajData          = [dataPath, 'lowTrajs_1mps_neckFull_L2_2km_22v0s_1yr.txt'];
    % file_lowImpactTrajData          = [dataPath, 'lowTrajs_1mps_nuckFull_L1_5km_22v0s_4pi.txt'];
    % file_lowImpactTrajData          = [dataPath, 'lowTrajs_1mps_neckFull_L1_5km_22v0s_1yr.txt'];


%     file_highImpactTrajData         = [dataPath, 'highImpactTrajs_100mps_1pi.txt'];
    % file_highImpactTrajData         = [dataPath, 'highImpactTrajs_150mps_1pi.txt'];
    % file_highImpactTrajData         = [dataPath, 'highImpactTrajs_stopAtImpact_150mps_4pi.txt'];
    
    %%% Load data
    lowImpactTrajData_BaCR          = dlmread(file_lowImpactTrajData,',',1,0);
    % lowImpactEvents_BaCR_periapses  = lowImpactTrajData_BaCR(lowImpactTrajData_BaCR(:,c_isPer) == true, :);
    % lowImpactEvents_BaCR_apoapses   = lowImpactTrajData_BaCR(lowImpactTrajData_BaCR(:,c_isApo) == true, :);
    % lowImpactEvents_BaCR_zEquals0   = lowImpactTrajData_BaCR(lowImpactTrajData_BaCR(:,c_isxyCross) == true, :);

%     highImpactTrajData_BaCR         = dlmread(file_highImpactTrajData,',',1,0);
    % highImpactEvents_BaCR_periapses = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,c_isPer) == true, :);
    % highImpactEvents_BaCR_apoapses  = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,c_isApo) == true, :);
    % highImpactEvents_BaCR_zEquals0  = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,c_isxyCross) == true, :);

%%% 1-year numerical proofs
elseif loadPlot_1yrData 
%     file_L1TrajData_1yr = [dataPath, 'lowTrajs_1mps_neckFull_L1_250km_22v0s_1yr.txt'];
%     file_L2TrajData_1yr = [dataPath, 'lowTrajs_1mps_neckFull_L2_2km_22v0s_1yr.txt'];
    file_L2TrajData_1yr = [dataPath, 'lowImpactTrajsNew.txt'];

%     lowData_L1_1yr = dlmread(file_L1TrajData_1yr,',',1,0); % min(min(abs(i)-90)): 4.21e-4 ... corresponding to i=90.0004 (index 10754971..trajID 1065) |||| max(abs(lat)): 39.7
    lowData_L2_1yr = dlmread(file_L2TrajData_1yr,',',1,0);
end


% % % % figure; plot(lowData_L2_1yr(:,c_i),'.')
% % % indices = find(abs(abs(lowData_L2_1yr(:,c_i))-90) < 1);
% % % figure; hold all
% % % axis equal
% % % plot3(lowData_L2_1yr(indices,c_x),lowData_L2_1yr(indices,c_y),lowData_L2_1yr(indices,c_z),'k.')
% % % PlotBoi3_CR3Bn(20)
% % % plotSecondary(secondary)
% % % quiver3(lowData_L2_1yr(indices,c_x),lowData_L2_1yr(indices,c_y),lowData_L2_1yr(indices,c_z),lowData_L2_1yr(indices,c_xd),lowData_L2_1yr(indices,c_yd),lowData_L2_1yr(indices,c_zd))



% -------------------------------------------------
%%% For scanning columns of huge files
% -------------------------------------------------
% % % % fid = fopen(file_lowImpactTrajData);
% % % % n = 0;
% % % % maxLat = 0;
% % % % while ischar(tline)
% % % %     n = n+1;
% % % % %     disp(tline)
% % % % 
% % % %     % If past the header, split the line and check the latitude
% % % %     if n > 1
% % % %         x = strsplit(tline,',');
% % % %         lat = str2double(x{c_lat});
% % % %         if abs(lat) > maxLat
% % % %             maxLat = abs(lat);
% % % %         end
% % % % %         disp(lat)
% % % %     end
% % % %     
% % % %     % Print status
% % % %     if mod(n, 100000) == 0
% % % %         fprintf('%1.0d\n',n)
% % % %     end
% % % %     
% % % %     tline = fgetl(fid);
% % % % 
% % % % end
% % % % fclose(fid);
% % % % 
% % % % %%% Number of unique trajectories
% % % % n_lowTrajs  = sum(isnan(lowImpactTrajData_BaCR(:,1)));
% % % % % n_highTrajs = sum(isnan(highImpactTrajData_BaCR(:,1)));

% -------------------------------------------------
%%% System
% -------------------------------------------------
primary = bodies.jupiter;   secondary = bodies.europa;

% -------------------------------------------------
%%% Normalizing factors and equillibrium points
% -------------------------------------------------
%%% Normalizing factors
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Collinear equilibrium points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates




if ~loadPlot_1yrData
    r_L1s = NaN(size(lowImpactTrajData_BaCR,1),1);
    for kk = 1:size(lowImpactTrajData_BaCR,1)
        r_L1s(kk) = norm(lowImpactTrajData_BaCR(kk,c_x:c_z) - L123(1,:));
    end
    close_to_L1 = find(r_L1s == min(r_L1s));
    norm(lowImpactTrajData_BaCR(close_to_L1(1),c_xd:c_zd))*vNorm*1000

    t2 = find(lowImpactTrajData_BaCR(:,c_trajID)==25);
    t2_1 = t2(1);
    norm(lowImpactTrajData_BaCR(1,c_x:c_z) - lowImpactTrajData_BaCR(t2_1,c_x:c_z))*rNorm
end
% % % % -------------------------------------------------
% % % %%% Store trajectory data in separate cell arrays
% % % % -------------------------------------------------
% % % %%% Initialize cell arrays
% % % CA_lowImpactTrajData_BaCR          = cell(n_lowTrajs,1);
% % % CA_lowImpactEvents_BaCR_periapses  = cell(n_lowTrajs,1);
% % % CA_lowImpactEvents_BaCR_apoapses   = cell(n_lowTrajs,1);
% % % CA_lowImpactEvents_BaCR_zEquals0   = cell(n_lowTrajs,1);
% CA_highImpactTrajData_BaCR         = cell(n_highTrajs,1);
% CA_highImpactEvents_BaCR_periapses = cell(n_highTrajs,1);
% CA_highImpactEvents_BaCR_apoapses  = cell(n_highTrajs,1);
% CA_highImpactEvents_BaCR_zEquals0  = cell(n_highTrajs,1);
% % % 
% % % %%% Loop through trajectories and store the data in cell arrays
% % % for trajID = 1:n_lowTrajs
% % %     CA_lowImpactTrajData_BaCR{trajID}         = lowImpactTrajData_BaCR(lowImpactTrajData_BaCR(:,c_trajID) == trajID,:);
% % %     CA_lowImpactEvents_BaCR_periapses{trajID} = lowImpactEvents_BaCR_periapses(lowImpactEvents_BaCR_periapses(:,c_trajID) == trajID,:);
% % %     CA_lowImpactEvents_BaCR_apoapses{trajID}  = lowImpactEvents_BaCR_apoapses(lowImpactEvents_BaCR_apoapses(:,c_trajID) == trajID,:);
% % %     CA_lowImpactEvents_BaCR_zEquals0{trajID}  = lowImpactEvents_BaCR_zEquals0(lowImpactEvents_BaCR_zEquals0(:,c_trajID) == trajID,:);
% % % end
% % % 
% for trajID = 1:n_highTrajs
%     CA_highImpactTrajData_BaCR{trajID}         = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,c_trajID) == trajID,:);
%     CA_highImpactEvents_BaCR_periapses{trajID} = highImpactEvents_BaCR_periapses(highImpactEvents_BaCR_periapses(:,c_trajID) == trajID,:);
%     CA_highImpactEvents_BaCR_apoapses{trajID}  = highImpactEvents_BaCR_apoapses(highImpactEvents_BaCR_apoapses(:,c_trajID) == trajID,:);
%     CA_highImpactEvents_BaCR_zEquals0{trajID}  = highImpactEvents_BaCR_zEquals0(highImpactEvents_BaCR_zEquals0(:,c_trajID) == trajID,:);
% end

% % -------------------------------------------------
% %%% Combining data in large matrices for mass-plotting
% % -------------------------------------------------
% all_lowImpactTrajData_BaCR         = cell2mat(CA_lowImpactTrajData_BaCR);
% all_lowImpactEvents_BaCR_periapses = cell2mat(CA_lowImpactEvents_BaCR_periapses);
% all_lowImpactEvents_BaCR_apoapses  = cell2mat(CA_lowImpactEvents_BaCR_apoapses);
% all_lowImpactEvents_BaCR_zEquals0  = cell2mat(CA_lowImpactEvents_BaCR_zEquals0);
% 
% all_highImpactTrajData_BaCR         = cell2mat(CA_highImpactTrajData_BaCR);
% all_highImpactEvents_BaCR_periapses = cell2mat(CA_highImpactEvents_BaCR_periapses);
% all_highImpactEvents_BaCR_apoapses  = cell2mat(CA_highImpactEvents_BaCR_apoapses);
% all_highImpactEvents_BaCR_zEquals0  = cell2mat(CA_highImpactEvents_BaCR_zEquals0);

% ========================================================================
%%% Plotting isolated trajectories
% ========================================================================
if loadPlot_1yrData
    989
end
return
% -------------------------------------------------
%%% Discrete low trajs
% -------------------------------------------------
% figure; hold all
% title('Low Trajs')
% PlotBoi3_CR3Bn(20)
% axis equal
% plotSecondary(secondary)
% for kk = 1:1
%     plot3(CA_lowImpactTrajData_BaCR{kk}(:,c_x),CA_lowImpactTrajData_BaCR{kk}(:,c_y),CA_lowImpactTrajData_BaCR{kk}(:,c_z),'color',colors.blue,'linewidth',2)
%     plot3(CA_lowImpactEvents_BaCR_periapses{kk}(:,c_x),CA_lowImpactEvents_BaCR_periapses{kk}(:,c_y),CA_lowImpactEvents_BaCR_periapses{kk}(:,c_z),'k.','markersize',18)
%     plot3(CA_lowImpactEvents_BaCR_apoapses{kk}(:,c_x),CA_lowImpactEvents_BaCR_apoapses{kk}(:,c_y),CA_lowImpactEvents_BaCR_apoapses{kk}(:,c_z),'k.','markersize',18)
% %     plot3(CA_lowImpactEvents_zEquals0{kk}(:,c_x),CA_lowImpactEvents_zEquals0{kk}(:,c_y),CA_lowImpactEvents_zEquals0{kk}(:,c_z),'k.','markersize',18)
% end
% -------------------------------------------------
%%% Discrete low trajs
% -------------------------------------------------
if plot_allLowTrajs
    figure; hold all
    title('Low Trajs')
    PlotBoi3_CR3Bn(20)
    axis equal
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
    view(0,0)

    plotRange_data = 1:size(lowImpactTrajData_BaCR,1);
%     plotRange_per = 1:size(highImpactEvents_BaCR_periapses,1);
%     plotRange_apo = 1:size(highImpactEvents_BaCR_apoapses,1);

    plot3(lowImpactTrajData_BaCR(plotRange_data,c_x),lowImpactTrajData_BaCR(plotRange_data,c_y),lowImpactTrajData_BaCR(plotRange_data,c_z),'color',colors.blue,'linewidth',1)
end % plot_allHighTrajs

% -------------------------------------------------
%%% Discrete high trajs
% -------------------------------------------------
if plot_allHighTrajs
    figure; hold all
    title('High Trajs')
    PlotBoi3_CR3Bn(20)
    axis equal
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
    view(0,0)
    % % for kk = 1:length(CA_highImpactTrajData_BaCR)
    % %     plot3(CA_highImpactTrajData_BaCR{kk}(:,c_x),CA_highImpactTrajData_BaCR{kk}(:,c_y),CA_highImpactTrajData_BaCR{kk}(:,c_z),'color',colors.red,'linewidth',1)
    % %     plot3(CA_highImpactEvents_BaCR_periapses{kk}(:,c_x),CA_highImpactEvents_BaCR_periapses{kk}(:,c_y),CA_highImpactEvents_BaCR_periapses{kk}(:,c_z),'b.','markersize',18)
    % %     plot3(CA_highImpactEvents_BaCR_apoapses{kk}(:,c_x),CA_highImpactEvents_BaCR_apoapses{kk}(:,c_y),CA_highImpactEvents_BaCR_apoapses{kk}(:,c_z),'g.','markersize',18)
    % % %     plot3(CA_highImpactEvents_zEquals0{kk}(:,c_x),CA_highImpactEvents_zEquals0{kk}(:,c_y),CA_highImpactEvents_zEquals0{kk}(:,c_z),'k.','markersize',18)
    % % end

    plotRange_data = 1:size(highImpactTrajData_BaCR,1);
    plotRange_per = 1:size(highImpactEvents_BaCR_periapses,1);
    plotRange_apo = 1:size(highImpactEvents_BaCR_apoapses,1);

    plot3(highImpactTrajData_BaCR(plotRange_data,c_x),highImpactTrajData_BaCR(plotRange_data,c_y),highImpactTrajData_BaCR(plotRange_data,c_z),'color',colors.red,'linewidth',1)
    plot3(highImpactEvents_BaCR_periapses(plotRange_per,c_x),highImpactEvents_BaCR_periapses(plotRange_per,c_y),highImpactEvents_BaCR_periapses(plotRange_per,c_z),'b.','markersize',18)
    plot3(highImpactEvents_BaCR_apoapses(plotRange_apo,c_x),highImpactEvents_BaCR_apoapses(plotRange_apo,c_y),highImpactEvents_BaCR_apoapses(plotRange_apo,c_z),'g.','markersize',18)
    % plot3(highImpactEvents_BaCR_zEquals0(plotRange,c_x),highImpactEvents_BaCR_zEquals0(plotRange,c_y),highImpactEvents_BaCR_zEquals0(plotRange,c_z),'k.','markersize',18)
end % plot_allHighTrajs

% -------------------------------------------------
%%% All trajs
% -------------------------------------------------
% if plot_allTrajs
%     figure; hold all
%     title('All Trajs')
%     PlotBoi3_CR3Bn(20)
%     axis equal
%     plotSecondary(secondary)
%     for kk = 1:n_lowTrajs
%         plot3(CA_lowImpactTrajData_BaCR{kk}(:,c_x),CA_lowImpactTrajData_BaCR{kk}(:,c_y),CA_lowImpactTrajData_BaCR{kk}(:,c_z),'color',colors.blue,'linewidth',2)
%     end
%     for kk = 1:n_highTrajs
%         plot3(CA_highImpactTrajData_BaCR{kk}(:,c_x),CA_highImpactTrajData_BaCR{kk}(:,c_y),CA_highImpactTrajData_BaCR{kk}(:,c_z),'color',colors.red,'linewidth',2)
%     end
% 
%     xlim([0.97, 1.03])
%     ylim([-1, 1].*0.02)
% end % plot_allTrajs

% ========================================================================
%%% Plotting full sets
% ========================================================================
% -------------------------------------------------
%%% Periapse Mapping
% -------------------------------------------------
if plot_periapseMap
    figure; hold all
    title('Periapse Mapping')
    PlotBoi3_CR3Bn(20)
    axis equal
    p_low  = plot3(lowImpactEvents_BaCR_periapses(:,c_x),lowImpactEvents_BaCR_periapses(:,c_y),lowImpactEvents_BaCR_periapses(:,c_z),'b.','markersize',4);
    p_high = plot3(highImpactEvents_BaCR_periapses(:,c_x),highImpactEvents_BaCR_periapses(:,c_y),highImpactEvents_BaCR_periapses(:,c_z),'r.','markersize',4);
    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
end % plot_periapseMap
% -------------------------------------------------
%%% Apoapse Mapping
% -------------------------------------------------
if plot_apoapseMap
    figure; hold all
    title('Apoapse Mapping')
    PlotBoi3_CR3Bn(20)
    axis equal
    p_low  = plot3(lowImpactEvents_BaCR_apoapses(:,c_x),lowImpactEvents_BaCR_apoapses(:,c_y),lowImpactEvents_BaCR_apoapses(:,c_z),'b.','markersize',4);
    p_high = plot3(highImpactEvents_BaCR_apoapses(:,c_x),highImpactEvents_BaCR_apoapses(:,c_y),highImpactEvents_BaCR_apoapses(:,c_z),'r.','markersize',4);
    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
end % plot_apoapseMap
% -------------------------------------------------
%%% |R| at Periapsis Mapping
% -------------------------------------------------
if plot_RPeriapse
    RS_lowImpactEvents_periapsis = rowNorm(lowImpactEvents_BaCR_periapses(:,c_x:c_z) - [1-secondary.MR,0,0]);
    RS_highImpactEvents_periapsis = rowNorm(highImpactEvents_BaCR_periapses(:,c_x:c_z) - [1-secondary.MR,0,0]);

%     figure; hold all
%     title('|R| of Periapses')
%     PlotBoi2('RS$_{per}$','RS$_{per}$',20,'LaTex')
%     axis equal
%     p_low  = plot(RS_lowImpactEvents_periapsis,RS_lowImpactEvents_periapsis,'b.','markersize',8);
%     p_high = plot(RS_highImpactEvents_periapsis,RS_highImpactEvents_periapsis,'r.','markersize',8);
% 
%     p_rad  = plot([0, secondary.R_n],[secondary.R_n, secondary.R_n],'k','linewidth',1);
%     plot([secondary.R_n, secondary.R_n],[0, secondary.R_n],'k','linewidth',1)
% 
%     xlim([0, max([max(RS_highImpactEvents_periapsis), max(RS_lowImpactEvents_periapsis)])])
%     ylim([0, max([max(RS_highImpactEvents_periapsis), max(RS_lowImpactEvents_periapsis)])])
% 
%     [legh,objh] = legend([p_low, p_high, p_rad],'Low','High', 'Europa Radius');
%     lineh = findobj(objh,'type','line');
%     set(lineh,'linestyle','-','linewidth',3);
    figure; hold all
    title('|R| of Periapse vs latitude')
    PlotBoi2('RS$_{per}$','Latitude, $^\circ$',20,'LaTex')

    p_low  = plot(RS_lowImpactEvents_periapsis,lowImpactEvents_BaCR_periapses(:,c_lat),'b.','markersize',8);
    p_high = plot(RS_highImpactEvents_periapsis,highImpactEvents_BaCR_periapses(:,c_lat),'r.','markersize',8);

    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);

end
% -------------------------------------------------
%%% |R| at Apoapsis Mapping
% -------------------------------------------------
if plot_RApoapse
    RS_lowImpactEvents_apoapsis  = rowNorm(lowImpactEvents_BaCR_apoapses(:,c_x:c_z) - [1-secondary.MR,0,0]);
    RS_highImpactEvents_apoapsis  = rowNorm(highImpactEvents_BaCR_apoapses(:,c_x:c_z) - [1-secondary.MR,0,0]);

%     figure; hold all
%     title('|R| of Apoapses')
%     PlotBoi2('RS$_{apo}$','RS$_{apo}$',20,'LaTex')
%     axis equal
%     p_low  = plot(RS_lowImpactEvents_apoapsis,RS_lowImpactEvents_apoapsis,'b.','markersize',8);
%     p_high = plot(RS_highImpactEvents_apoapsis,RS_highImpactEvents_apoapsis,'r.','markersize',8);
% 
%     p_rad  = plot([0, secondary.R_n],[secondary.R_n, secondary.R_n],'k','linewidth',1);
%     plot([secondary.R_n, secondary.R_n],[0, secondary.R_n],'k','linewidth',1)
% 
%     xlim([0, max([max(RS_highImpactEvents_apoapsis), max(RS_lowImpactEvents_apoapsis)])])
%     ylim([0, max([max(RS_highImpactEvents_apoapsis), max(RS_lowImpactEvents_apoapsis)])])
% 
% 
%     [legh,objh] = legend([p_low, p_high, p_rad],'Low','High', 'Europa Radius');
%     lineh = findobj(objh,'type','line');
%     set(lineh,'linestyle','-','linewidth',3);
    figure; hold all
    title('|R| of Apoapses vs latitude')
    PlotBoi2('RS$_{apo}$','Latitude, $^\circ$',20,'LaTex')

    p_low  = plot(RS_lowImpactEvents_apoapsis,lowImpactEvents_BaCR_apoapses(:,c_lat),'b.','markersize',8);
    p_high = plot(RS_highImpactEvents_apoapsis,highImpactEvents_BaCR_apoapses(:,c_lat),'r.','markersize',8);

    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);

end % plot_RApoapse

% -------------------------------------------------
%%% X-Y Plane Crossing
% -------------------------------------------------
% if plot_xyCrossing
%     figure; hold all
%     title('X-Y Plane Crossing')
%     PlotBoi3_CR3Bn(20)
%     axis equal
%     p_low  = plot(lowImpactEvents_BaCR_zEquals0(:,c_x),lowImpactEvents_BaCR_zEquals0(:,c_y),'b.','markersize',4);
%     p_high = plot(highImpactEvents_BaCR_zEquals0(:,c_x),highImpactEvents_BaCR_zEquals0(:,c_y),'r.','markersize',4);
%     [legh,objh] = legend([p_low, p_high],'Low','High');
%     lineh = findobj(objh,'type','line');
%     set(lineh,'linestyle','-','linewidth',3);
% end

% -------------------------------------------------
%%% X-Y Plane Crossing - 3D
% -------------------------------------------------
if plot_xyCrossing_3D    
    figure; hold all
    title('X-Y Plane Crossing w/ a')
    axis normal
    PlotBoi3('$x_n$','$y_n$','$a$',20,'LaTeX')

    p_low  = plot3(lowImpactEvents_BaCR_zEquals0(:,c_x),lowImpactEvents_BaCR_zEquals0(:,c_y), lowImpactEvents_BaCR_zEquals0(:,c_a),'b.','markersize',4);
    p_high = plot3(highImpactEvents_BaCR_zEquals0(:,c_x),highImpactEvents_BaCR_zEquals0(:,c_y),highImpactEvents_BaCR_zEquals0(:,c_a),'r.','markersize',4);

    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    figure; hold all
    title('X-Y Plane Crossing w/ i')
    axis normal
    PlotBoi3('$x_n$','$y_n$','$i$',20,'LaTeX')

    p_low  = plot3(lowImpactEvents_BaCR_zEquals0(:,c_x),lowImpactEvents_BaCR_zEquals0(:,c_y), lowImpactEvents_BaCR_zEquals0(:,c_i),'b.','markersize',4);
    p_high = plot3(highImpactEvents_BaCR_zEquals0(:,c_x),highImpactEvents_BaCR_zEquals0(:,c_y),highImpactEvents_BaCR_zEquals0(:,c_i),'r.','markersize',4);

    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% Testing
% -------------------------------------------------
% figure; hold all
% title('|R| of Apoapses vs latitude')
% PlotBoi2('RS$_{apo}$','Latitude, $^\circ$',20,'LaTex')
% 
% p_low  = plot(RS_lowImpactEvents_apoapsis,lowImpactEvents_BaCR_apoapses(:,c_lat),'b.','markersize',8);
% p_high = plot(RS_highImpactEvents_apoapsis,highImpactEvents_BaCR_apoapses(:,c_lat),'r.','markersize',8);
% 
% [legh,objh] = legend([p_low, p_high],'Low','High');
% lineh = findobj(objh,'type','line');
% set(lineh,'linestyle','-','linewidth',3);


if ~loadPlot_1yrData
    %%% i vs a vs H
    figure; hold all
    PlotBoi3('$i$','$a$','$|h|$',20,'LaTex')
    ylim([6, 10].*1e-3)
    zlim([0 0.0006])
    p_low  = plot3(lowImpactTrajData_BaCR(:,c_i), lowImpactTrajData_BaCR(:,c_a), lowImpactTrajData_BaCR(:,c_H), 'b.', 'markersize',3);
    p_high = plot3(highImpactTrajData_BaCR(:,c_i), highImpactTrajData_BaCR(:,c_a), highImpactTrajData_BaCR(:,c_H), 'r.', 'markersize',3);
    [legh,objh] = legend([p_low, p_high],'Low','High');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end
% ========================================================================
%%% Checking end-times of high trajs to see if most were impacting or not
% ========================================================================

% % average end time of trajectories
% endtimes = NaN(n_highTrajs,1);
% for kk = 1:n_highTrajs
%     endtimes(kk) = CA_highImpactTrajData_BaCR{kk}(end,c_t);
% end
% max(endtimes)
% mean(endtimes)

% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)

% % % savepath = '/Users/lukebury/CU_Google_Drive/Documents/JulGit/Research/HighLatitudeLanding/GoodFullData';
% % % filename_high = fullfile(savepath, 'highImpactTrajsZEquals0_full.txt');
% % % filename_low  = fullfile(savepath, 'lowImpactTrajsZEquals0_full.txt');
% % %  
% % % f_low  = fopen(filename_low, 'wt');
% % % f_high = fopen(filename_high, 'wt');
% % %  
% % % fprintf(f_low, 'trajID,x,y,z,xd,yd,zd,t,JC,R,a,e,i,E,H,h_x,h_y,h_z,lat,lon\n');
% % % fprintf(f_high, 'trajID,x,y,z,xd,yd,zd,t,JC,R,a,e,i,E,H,h_x,h_y,h_z,lat,lon\n');
% % %  
% % % for kk = 1:size(all_lowImpactEvents_BaCR_zEquals0_allData,1)
% % %      fprintf(f_low, '%1.0d,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.3f,%1.3f\n',...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,1), all_lowImpactEvents_BaCR_zEquals0_allData(kk,2), all_lowImpactEvents_BaCR_zEquals0_allData(kk,3),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,4), all_lowImpactEvents_BaCR_zEquals0_allData(kk,5), all_lowImpactEvents_BaCR_zEquals0_allData(kk,6),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,7), all_lowImpactEvents_BaCR_zEquals0_allData(kk,8), all_lowImpactEvents_BaCR_zEquals0_allData(kk,9),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,10), all_lowImpactEvents_BaCR_zEquals0_allData(kk,11), all_lowImpactEvents_BaCR_zEquals0_allData(kk,12),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,13), all_lowImpactEvents_BaCR_zEquals0_allData(kk,14), all_lowImpactEvents_BaCR_zEquals0_allData(kk,15),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,16), all_lowImpactEvents_BaCR_zEquals0_allData(kk,17), all_lowImpactEvents_BaCR_zEquals0_allData(kk,18),...
% % %          all_lowImpactEvents_BaCR_zEquals0_allData(kk,19), all_lowImpactEvents_BaCR_zEquals0_allData(kk,20));
% % % end
% % %  
% % % for kk = 1:size(all_highImpactEvents_BaCR_zEquals0_allData,1)
% % %      fprintf(f_high, '%1.0d,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.3f,%1.3f\n',...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,1), all_highImpactEvents_BaCR_zEquals0_allData(kk,2), all_highImpactEvents_BaCR_zEquals0_allData(kk,3),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,4), all_highImpactEvents_BaCR_zEquals0_allData(kk,5), all_highImpactEvents_BaCR_zEquals0_allData(kk,6),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,7), all_highImpactEvents_BaCR_zEquals0_allData(kk,8), all_highImpactEvents_BaCR_zEquals0_allData(kk,9),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,10), all_highImpactEvents_BaCR_zEquals0_allData(kk,11), all_highImpactEvents_BaCR_zEquals0_allData(kk,12),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,13), all_highImpactEvents_BaCR_zEquals0_allData(kk,14), all_highImpactEvents_BaCR_zEquals0_allData(kk,15),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,16), all_highImpactEvents_BaCR_zEquals0_allData(kk,17), all_highImpactEvents_BaCR_zEquals0_allData(kk,18),...
% % %          all_highImpactEvents_BaCR_zEquals0_allData(kk,19), all_highImpactEvents_BaCR_zEquals0_allData(kk,20));
% % % end

% RS = rowNorm(highImpactTrajData_BaCR(:,c_x:c_z) - [1-secondary.MR,0,0]);
% indx_minR = find(RS==min(RS));
% susTrajID = highImpactTrajData_BaCR(indx_minR,1);
% 
% % indx_minR = find(RS_highIimpactEvents_periapsis==min(RS_highIimpactEvents_periapsis));
% 
% % susTrajID = highImpactEvents_BaCR_periapses(indx_minR,1);
% 
% susTraj = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,1)==susTrajID,:);
% 
% 
% figure; hold all
% PlotBoi3_CR3Bn(20)
% axis equal
% plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
% plot3(susTraj(:,2),susTraj(:,3),susTraj(:,4),'r')
% 
% 
% 
% 
% interestingSolIndx = find(highImpactTrajData_BaCR(:,c_i) == min(highImpactTrajData_BaCR(:,c_i)));
% interestingTrajID  = highImpactTrajData_BaCR(interestingSolIndx,1);
% interestingTraj    = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,1)==interestingTrajID,:);
% figure; hold all;
% title('min i')
% plot3(interestingTraj(:,2),interestingTraj(:,3),interestingTraj(:,4),'r')
% plotSecondary(secondary)
% axis equal
% PlotBoi3_CR3Bn(20)
% plot3(L123(1:2,1),zeros(2),zeros(2),'^','markeredgecolor',colors.grn,'markerfacecolor',colors.ltgrn,'markersize',10)
% 
% 
% interestingSolIndx = find(highImpactTrajData_BaCR(:,c_x) == max(highImpactTrajData_BaCR(:,c_x)));
% interestingTrajID  = highImpactTrajData_BaCR(interestingSolIndx,1);
% interestingTraj    = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,1)==interestingTrajID,:);
% figure; hold all;
% title('max x')
% plot3(interestingTraj(:,2),interestingTraj(:,3),interestingTraj(:,4),'r')
% plotSecondary(secondary)
% axis equal
% PlotBoi3_CR3Bn(20)
% plot3(L123(1:2,1),zeros(2),zeros(2),'^','markeredgecolor',colors.grn,'markerfacecolor',colors.ltgrn,'markersize',10)


% interestingSolIndx = find(highImpactTrajData_BaCR(:,c_x) == min(highImpactTrajData_BaCR(:,c_x)));
% interestingTrajID  = highImpactTrajData_BaCR(interestingSolIndx,1);
% interestingTraj    = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,1)==interestingTrajID,:);
% figure; hold all;
% title('min x')
% plot3(interestingTraj(:,2),interestingTraj(:,3),interestingTraj(:,4),'r')
% plotSecondary(secondary)
% axis equal
% PlotBoi3_CR3Bn(20)
% plot3(L123(1:2,1),zeros(2),zeros(2),'^','markeredgecolor',colors.grn,'markerfacecolor',colors.ltgrn,'markersize',10)


% interestingTrajID  = 36356;
% interestingTraj    = highImpactTrajData_BaCR(highImpactTrajData_BaCR(:,1)==interestingTrajID,:);
% figure; hold all;
% title('Weird, constant inclination')
% plot3(interestingTraj(:,2),interestingTraj(:,3),interestingTraj(:,4),'r')
% plotSecondary(secondary)
% axis equal
% PlotBoi3_CR3Bn(20)
% plot3(L123(1:2,1),zeros(2),zeros(2),'^','markeredgecolor',colors.grn,'markerfacecolor',colors.ltgrn,'markersize',10)
% 

% figure; hold all;
% plot3(lowImpactTrajData_BaCR(:,c_x),lowImpactTrajData_BaCR(:,c_y),lowImpactTrajData_BaCR(:,c_z),'b')
% plotSecondary(secondary)
% axis equal
% PlotBoi3_CR3Bn(20)
% plot3(L123(1:3,1),zeros(3),zeros(3),'^','markeredgecolor',colors.grn,'markerfacecolor',colors.ltgrn,'markersize',10)





