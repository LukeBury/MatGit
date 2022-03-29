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
dataPath = '/Volumes/LB_External_Drive/Research/High_Latitude_Landing_Study/Data/Full_Results/';
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
% ========================================================================
%%% Setup
% ========================================================================
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
L123 = EquilibriumPoints(secondary.MR,1,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

% -------------------------------------------------
%%% Load the data
% -------------------------------------------------
%%% Column specifiers
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
c_L_del     = 24;
c_G_del     = 25;
c_H_del     = 26;
c_isL1Cross = 27;



study_1   = false;
study_50  = false;
study_100 = false;
study_150 = true;

% tic
% data_L2_1yr = dlmread([dataPath, 'lowTrajs_50mps_neckFull_L2_125km_22v0s_1yr.txt'],',',1,0);
% toc
% save([dataPath,'lowTrajs_50mps_neckFull_L2_125km_22v0s_1yr','.mat'], 'data_L2_1yr','-v7.3')
% toc
% tic
% data_L2_1yr = dlmread([dataPath, 'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr.txt'],',',1,0);
% toc
% save([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr','.mat'], 'data_L2_1yr','-v7.3')
% toc
% return

if study_1
%     load([dataPath,'lowTrajs_1mps_neckFull_L1_250km_22v0s_4pi.mat']); % data_L1_4pi
    load([dataPath,'lowTrajs_1mps_neckFull_L1_250km_22v0s_1yr.mat']); % data_L1_1yr
%     load([dataPath,'lowTrajs_1mps_neckFull_L1_250km_22v0s_2yr.mat']); % data_L1_2yr
%     load([dataPath,'lowTrajs_1mps_neckFull_L2_2km_22v0s_4pi.mat']); % data_L2_4pi
%     load([dataPath,'lowTrajs_1mps_neckImpact_L2_2km_22v0s_4pi.mat']); % data_L2_4pi_impact
    load([dataPath,'lowTrajs_1mps_neckFull_L2_2km_22v0s_1yr.mat']); % data_L2_1yr
    load([dataPath,'highTrajs_1mps_70lat_1pi.mat']); % data_high_1pi
    data_title = '1 m/s';
    
elseif study_50
    load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_4pi.mat']); % data_L1_4pi
    load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr.mat']); % data_L1_1yr
%     load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_2yr.mat']); % data_L1_2yr
%     load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_3yr.mat']); % data_L1_3yr
    load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi.mat']); % data_L2_4pi
    load([dataPath,'lowTrajs_50mps_neckImpact_L2_250km_22v0s_4pi.mat']); % data_L2_4pi_impact
    load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr.mat']); % data_L2_1yr
%     load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_2yr.mat']); % data_L2_2yr
%     load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_3yr.mat']); % data_L2_3yr
%     load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_3yr_L1Crossings.mat']); % data_L2_3yr_L1Crossings
%     load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_4yr.mat']); % data_L2_4yr
    load([dataPath,'highTrajs_50mps_70lat_1pi.mat']); % data_high_1pi
    data_title = '50 m/s';
    
elseif study_100
    load([dataPath,'lowTrajs_100mps_neckFull_L1_250km_22v0s_4pi.mat']); % data_L1_4pi
    load([dataPath,'lowTrajs_100mps_neckFull_L1_250km_22v0s_1yr.mat']); % data_L1_1yr
    load([dataPath,'lowTrajs_100mps_neckFull_L2_250km_22v0s_4pi.mat']); % data_L2_4pi
    load([dataPath,'lowTrajs_100mps_neckImpact_L2_250km_22v0s_4pi.mat']); % data_L2_4pi_impact
    load([dataPath,'lowTrajs_100mps_neckFull_L2_250km_22v0s_1yr.mat']); % data_L2_1yr
    load([dataPath,'highTrajs_100mps_70lat_1pi.mat']); % data_high_1pi
    data_title = '100 m/s';
    
elseif study_150
    load([dataPath,'lowTrajs_150mps_neckFull_L1_250km_22v0s_4pi.mat']); % data_L1_4pi
    load([dataPath,'lowTrajs_150mps_neckFull_L1_250km_22v0s_1yr.mat']); % data_L1_1yr
    load([dataPath,'lowTrajs_150mps_neckFull_L2_250km_22v0s_4pi.mat']); % data_L2_4pi
    load([dataPath,'lowTrajs_150mps_neckImpact_L2_250km_22v0s_4pi.mat']); % data_L2_4pi_impact
    load([dataPath,'lowTrajs_150mps_neckFull_L2_250km_22v0s_1yr.mat']); % data_L2_1yr
    load([dataPath,'highTrajs_150mps_70lat_1pi.mat']); % data_high_1pi
    data_title = '150 m/s';
end

label_L1_4pi        = 'L1, 4pi';
label_L1_1yr        = 'L1, 1 year';
label_L2_4pi        = 'L2, 4pi';
label_L2_4pi_impact = 'L2, 4pi - impact';
label_L2_1yr        = 'L2, 1 year';
label_high_1pi      = 'Lat 70°+, 1pi';

color_L1_4pi        = colors.grn2;
color_L1_1yr        = colors.ltgrn;
color_L2_4pi        = colors.blue;
color_L2_4pi_impact = colors.brown;
color_L2_1yr        = colors.cyan;
color_high_1pi      = colors.red;







% %%% Paths to data files
% % 1 mps
% file_low_1mps_L1_1yr          = [dataPath, 'lowTrajs_1mps_neckFull_L1_250km_22v0s_1yr.txt'];
% file_low_1mps_L1_4pi          = [dataPath, 'lowTrajs_1mps_neckFull_L1_250km_22v0s_4pi.txt'];
% file_low_1mps_L2_1yr          = [dataPath, 'lowTrajs_1mps_neckFull_L2_2km_22v0s_1yr.txt'];
% file_low_1mps_L2_4pi          = [dataPath, 'lowTrajs_1mps_neckFull_L2_2km_22v0s_4pi.txt'];
% file_low_1mps_L2_4pi_impact   = [dataPath, 'lowTrajs_1mps_neckImpact_L2_2km_22v0s_4pi.txt'];
% 
% % 50 mps
% file_low_50mps_L1_4pi         = [dataPath, 'lowTrajs_50mps_neckFull_L1_250km_22v0s_4pi.txt'];   % '50 m/s, L1, 4?'
% file_low_50mps_L2_4pi         = [dataPath, 'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi.txt'];   % '50 m/s, L2, 4?'
% file_low_50mps_L2_4pi_impact  = [dataPath, 'lowTrajs_50mps_neckImpact_L2_250km_22v0s_4pi.txt']; % '50 m/s, L2, 4? - impact'
% file_low_50mps_L2_1yr         = [dataPath, 'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr.txt'];   % '50 m/s, L2, 1 year'
% file_high_50mps_70lat_1pi     = [dataPath, 'highTrajs_50mps_70lat_1pi.txt'];                    % '50 m/s, 70°+, 1?'
% 
% % 100 mps
% file_low_100mps_L1_4pi        = [dataPath, 'lowTrajs_100mps_neckFull_L1_250km_22v0s_4pi.txt'];
% file_low_100mps_L2_4pi        = [dataPath, 'lowTrajs_100mps_neckFull_L2_250km_22v0s_4pi.txt'];
% file_low_100mps_L2_4pi_impact = [dataPath, 'lowTrajs_100mps_neckImpact_L2_250km_22v0s_4pi.txt'];
% file_low_100mps_L2_1yr        = [dataPath, 'lowTrajs_100mps_neckFull_L2_250km_22v0s_1yr.txt']; % highest lat never left L1 and to return.. it reached the high lat naturally.
% file_high_100mps_70lat_1pi    = [dataPath, 'highTrajs_100mps_70lat_1pi.txt']; % roughly 40 min
% 
% % 150 mps
% file_low_150mps_L1_4pi        = [dataPath, 'lowTrajs_150mps_neckFull_L1_250km_22v0s_4pi.txt'];
% file_low_150mps_L2_4pi        = [dataPath, 'lowTrajs_150mps_neckFull_L2_250km_22v0s_4pi.txt'];
% file_low_150mps_L2_4pi_impact = [dataPath, 'lowTrajs_150mps_neckImpact_L2_250km_22v0s_4pi.txt'];
% file_low_150mps_L2_1yr        = [dataPath, 'lowTrajs_150mps_neckFull_L2_250km_22v0s_1yr.txt']; % highest lat in trajID:2796 ... X0 below
% file_high_150mps_70lat_1pi    = [dataPath, 'highTrajs_150mps_70lat_1pi.txt'];
% 
% %%% Choose data files
% file_1 = file_low_50mps_L1_4pi;
% file_2 = file_low_50mps_L2_4pi;
% 
% 
% if contains(file_1,'low') && contains(file_2,'low')
%     data1_label = 'Low set 1';
%     data1_color = colors.mag;
%     
%     data2_label = 'Low set 2';
%     data2_color = colors.blue;
% elseif contains(file_1,'low') && contains(file_2,'high')
%     data1_label = 'Low';
%     data1_color = colors.blue;
%     
%     data2_label = 'High';
%     data2_color = colors.red;
% elseif contains(file_1,'high') && contains(file_2,'low')
%     data1_label = 'High';
%     data1_color = colors.red;
%     
%     data2_label = 'Low';
%     data2_color = colors.blue;
% end
% 
% data1_name = split(file_1, '/');
% data1_name = data1_name{end};
% data2_name = split(file_2, '/');
% data2_name = data2_name{end};
% 
% %%% Load the data into matrices
% data1_BaCR = dlmread(file_1,',',1,0);
% data2_BaCR = dlmread(file_2,',',1,0);
% 
% n_data1 = size(data1_BaCR,1);
% n_data2 = size(data2_BaCR,2);
% 
% nTrajs_data1 = data1_BaCR(end-1,c_trajID);
% nTrajs_data2 = data2_BaCR(end-1,c_trajID);
% 
% %%% Add optional third data set
% warning('Third dataset is on')
% file_3       = file_low_50mps_L2_4pi_impact;
% data3_label  = 'Low set 3';
% data3_color  = colors.grn2;
% data3_name   = split(file_3, '/');
% data3_name   = data3_name{end};
% data3_BaCR   = dlmread(file_3,',',1,0);
% n_data3      = size(data3_BaCR,1);
% nTrajs_data3 = data3_BaCR(end-1,c_trajID);
% 
% warning('Fourth dataset is on')
% file_4       = file_low_50mps_L2_1yr;
% data4_label  = 'Low set 4';
% data4_color  = colors.brown;
% data4_name   = split(file_4, '/');
% data4_name   = data4_name{end};
% data4_BaCR   = dlmread(file_4,',',1,0);
% n_data4      = size(data4_BaCR,1);
% nTrajs_data4 = data4_BaCR(end-1,c_trajID);
% 
% warning('Fifth dataset is on')
% file_5       = file_high_50mps_70lat_1pi;
% data5_label  = 'High set 1';
% data5_color  = colors.red;
% data5_name   = split(file_5, '/');
% data5_name   = data5_name{end};
% data5_BaCR   = dlmread(file_5,',',1,0);
% n_data5      = size(data5_BaCR,1);
% nTrajs_data5 = data5_BaCR(end-1,c_trajID);


% ========================================================================
%%% If desired, create matrices of specific events from each data set
% ========================================================================

% data1_BaCR_apoapses  = data1_BaCR(data1_BaCR(:,c_isApo) == true, :);
% data1_BaCR_zEquals0  = data1_BaCR(data1_BaCR(:,c_isxyCross) == true, :);
% 
% data2_BaCR_apoapses  = data2_BaCR(data2_BaCR(:,c_isApo) == true, :);
% data2_BaCR_zEquals0  = data2_BaCR(data2_BaCR(:,c_isxyCross) == true, :);

% ========================================================================
%%% Separate into cell arrays if you want to easily see data
%%% trajectory-by-trajectory
% ========================================================================
if 1+1==1
    %%% Initialize cell arrays
    CA_data1_BaCR           = cell(n_data1,1);
    CA_data2_BaCR           = cell(n_data2,1);

    %%% Loop through trajectories and store the data in cell arrays
    for trajID = 1:n_data1
        CA_data1_BaCR{trajID}  = data1_BaCR(data1_BaCR(:,c_trajID) == trajID,:);
    end
    for trajID = 1:n_data2
        CA_data2_BaCR{trajID}  = data2_BaCR(data2_BaCR(:,c_trajID) == trajID,:);
    end
end
% ========================================================================
%%% For scanning columns of huge files
% ========================================================================
if 1+1==1
    fid = fopen(file_lowImpactTrajData);
    n = 0;
    maxLat = 0;
    while ischar(tline)
        n = n+1;
    %     disp(tline)

        % If past the header, split the line and check the latitude
        if n > 1
            x = strsplit(tline,',');
            lat = str2double(x{c_lat});
            if abs(lat) > maxLat
                maxLat = abs(lat);
            end
    %         disp(lat)
        end

        % Print status
        if mod(n, 100000) == 0
            fprintf('%1.0d\n',n)
        end

        tline = fgetl(fid);

    end
    fclose(fid);
end

% ========================================================================
%%% Compare the data and analyze!
% ========================================================================
% -------------------------------------------------
%%% General trajectory info
% -------------------------------------------------
if 1+1==1
%     fprintf('============================\n')
%     fprintf('Data1 - %s\n',data1_name)
%     fprintf('%1.1f deg: Maximum latitude\n', max(abs(data1_BaCR(:,c_lat))))
%     fprintf('%1.1f deg: Closest inlination to 90\n', data1_BaCR(logical(abs(abs(data1_BaCR(:,c_i)-90) == min(abs(abs(data1_BaCR(:,c_i))-90)))), c_i))
%     fprintf('============================\n')
%     fprintf('Data2 - %s\n',data2_name)
%     fprintf('%1.1f deg: Maximum latitude\n', max(abs(data2_BaCR(:,c_lat))))
%     fprintf('%1.1f deg: Closest inlination to 90\n', data2_BaCR(logical(abs(abs(data2_BaCR(:,c_i)-90) == min(abs(abs(data2_BaCR(:,c_i))-90)))), c_i))
%     fprintf('============================\n')
    
    fprintf('============================\n')
    fprintf('%s\n',label_L1_4pi)
    fprintf('%1.1f deg: Maximum latitude\n', max(abs(data_L1_4pi(:,c_lat))))
    fprintf('============================\n')
    fprintf('%s\n',label_L1_1yr)
    fprintf('%1.1f deg: Maximum latitude\n', max(abs(data_L1_1yr(:,c_lat))))
    fprintf('============================\n')
    fprintf('%s\n',label_L2_4pi)
    fprintf('%1.1f deg: Maximum latitude\n', max(abs(data_L2_4pi(:,c_lat))))
    fprintf('============================\n')
    fprintf('%s\n',label_L2_4pi_impact)
    fprintf('%1.1f deg: Maximum latitude\n', max(abs(data_L2_4pi_impact(:,c_lat))))
    fprintf('============================\n')
    fprintf('%s\n',label_L2_1yr)
    fprintf('%1.1f deg: Maximum latitude\n', max(abs(data_L2_1yr(:,c_lat))))
    fprintf('============================\n')

end

% -------------------------------------------------
%%% Plot the trajectories
% -------------------------------------------------
if 1+1==1
    figure; hold all
%     title(data_title)
    PlotBoi3_CR3Bn(23)
    axis equal
    xlim([L123(1,1)-secondary.R_n, L123(2,1)+secondary.R_n])
    
    p_high_1pi = plot3(data_high_1pi(:,c_x),data_high_1pi(:,c_y),data_high_1pi(:,c_z), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot3(data_L1_1yr(:,c_x),data_L1_1yr(:,c_y),data_L1_1yr(:,c_z), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot3(data_L2_1yr(:,c_x),data_L2_1yr(:,c_y),data_L2_1yr(:,c_z), '.','markersize',1, 'color', color_L2_1yr);
%     p_L1_1yr   = plot3(data_L1_3yr(indicies_L1_1styr,c_x),data_L1_3yr(indicies_L1_1styr,c_y),data_L1_3yr(indicies_L1_1styr,c_z), '.','markersize',1, 'color', color_L1_4pi);
%     p_L2_1yr   = plot3(data_L2_3yr(indicies_L2_1styr,c_x),data_L2_3yr(indicies_L2_1styr,c_y),data_L2_3yr(indicies_L2_1styr,c_z), '.','markersize',1, 'color', color_L2_4pi);
    
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
%     data_L1_4pi
%     data_L1_1yr
%     data_L2_4pi
%     data_L2_4pi_impact
%     data_L2_1yr
%     data_high_1pi
    prms.u = secondary.MR;
    prms.R2_n = secondary.R_n;
    JC_data = getJacobiConstant_ZH(data_L1_4pi(1,c_x:c_zd), prms);
    plotCR3BP_YZNeck(JC_data, prms.u , 1, 0, prms, colors.black, 1.5)
    plotCR3BP_YZNeck(JC_data, prms.u , 2, 0, prms, colors.black, 1.5)
    
    figure; hold all
    title(data_title)
    PlotBoi3_CR3Bn(23)
    axis equal
    p_high_1pi = plot3(data_high_1pi(:,c_x),data_high_1pi(:,c_y),data_high_1pi(:,c_z), '.','color', color_high_1pi);
    p_L1_1styr = plot3(data_L1_3yr(indicies_L1_1styr,c_x), data_L1_3yr(indicies_L1_1styr,c_y), data_L1_3yr(indicies_L1_1styr, c_z), '.', 'color', colors.grn2);
    p_L1_2ndyr = plot3(data_L1_3yr(indicies_L1_2ndyr,c_x), data_L1_3yr(indicies_L1_2ndyr,c_y), data_L1_3yr(indicies_L1_2ndyr, c_z), '.', 'color', colors.ltgrn);
    p_L1_3rdyr = plot3(data_L1_3yr(indicies_L1_3rdyr,c_x), data_L1_3yr(indicies_L1_3rdyr,c_y), data_L1_3yr(indicies_L1_3rdyr, c_z), '.', 'color', colors.ylwgrn);
    [legh,objh] = legend([p_high_1pi, p_L1_1styr, p_L1_2ndyr, p_L1_3rdyr],label_high_1pi,'L1, 1st year', 'L1, 2nd year', 'L1, 3rd year');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    figure; hold all
    title(data_title)
    PlotBoi3_CR3Bn(23)
    axis equal
    p_high_1pi = plot3(data_high_1pi(:,c_x),data_high_1pi(:,c_y),data_high_1pi(:,c_z), '.','color', color_high_1pi);
    p_L2_1styr = plot3(data_L2_4yr(indicies_L2_1styr,c_x), data_L2_4yr(indicies_L2_1styr,c_y), data_L2_4yr(indicies_L2_1styr, c_z), '.', 'color', colors.blue);
    p_L2_2ndyr = plot3(data_L2_4yr(indicies_L2_2ndyr,c_x), data_L2_4yr(indicies_L2_2ndyr,c_y), data_L2_4yr(indicies_L2_2ndyr, c_z), '.', 'color', colors.blue2);
    p_L2_3rdyr = plot3(data_L2_4yr(indicies_L2_3rdyr,c_x), data_L2_4yr(indicies_L2_3rdyr,c_y), data_L2_4yr(indicies_L2_3rdyr, c_z), '.', 'color', colors.ltblue);
    p_L2_4thyr = plot3(data_L2_4yr(indicies_L2_4thyr,c_x), data_L2_4yr(indicies_L2_4thyr,c_y), data_L2_4yr(indicies_L2_4thyr, c_z), '.', 'color', colors.purp);
    [legh,objh] = legend([p_high_1pi, p_L2_1styr, p_L2_2ndyr, p_L2_3rdyr, p_L2_4thyr],label_high_1pi,'L2, 1st year', 'L2, 2nd year', 'L2, 3rd year', 'L2, 4th year');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
%     p_high_1pi = plot(data_high_1pi(:,c_lon),data_high_1pi(:,c_lat), '.','markersize',1, 'color', color_high_1pi);
%     p_L1_1styr = plot(data_L1_3yr(indicies_L1_1styr,c_lon), data_L1_3yr(indicies_L1_1styr,c_lat),'.','markersize',1,'color',colors.grn2);
%     p_L1_2ndyr = plot(data_L1_3yr(indicies_L1_2ndyr,c_lon), data_L1_3yr(indicies_L1_2ndyr,c_lat),'.','markersize',1,'color',colors.ltgrn);
%     p_L1_3rdyr = plot(data_L1_3yr(indicies_L1_3rdyr,c_lon), data_L1_3yr(indicies_L1_3rdyr,c_lat),'.','markersize',1,'color',colors.ylwgrn);
%     p_L2_1styr = plot(data_L2_4yr(indicies_L2_1styr,c_lon), data_L2_4yr(indicies_L2_1styr,c_lat),'.','markersize',1,'color',colors.blue);
%     p_L2_2ndyr = plot(data_L2_4yr(indicies_L2_2ndyr,c_lon), data_L2_4yr(indicies_L2_2ndyr,c_lat),'.','markersize',1,'color',colors.blue2);
%     p_L2_3rdyr = plot(data_L2_4yr(indicies_L2_3rdyr,c_lon), data_L2_4yr(indicies_L2_3rdyr,c_lat),'.','markersize',1,'color',colors.ltblue);
%     p_L2_4thyr = plot(data_L2_4yr(indicies_L2_4thyr,c_lon), data_L2_4yr(indicies_L2_4thyr,c_lat),'.','markersize',1,'color',colors.purp);
%     [legh,objh] = legend([p_high_1pi, p_L1_1styr, p_L1_2ndyr, p_L1_3rdyr, p_L2_1styr, p_L2_2ndyr, p_L2_3rdyr, p_L2_4thyr],label_high_1pi,'L1, 1st year', 'L1, 2nd year', 'L1, 3rd year', 'L2, 1st year', 'L2, 2nd year', 'L2, 3rd year', 'L2, 4th year');
%     lineh = findobj(objh,'type','line');
%     set(lineh,'linestyle','-','linewidth',3);
    
    
    figure; hold all
    title(data_title)
    PlotBoi3_CR3Bn(23)
    axis equal
    xlim([L123(1,1)-secondary.R_n, L123(2,1)+secondary.R_n])
    
    p_high_1pi = plot3(data_high_1pi(:,c_x),data_high_1pi(:,c_y),data_high_1pi(:,c_z), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot3(data_L1_4pi(:,c_x),data_L1_4pi(:,c_y),data_L1_4pi(:,c_z), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot3(data_L2_4pi(:,c_x),data_L2_4pi(:,c_y),data_L2_4pi(:,c_z), '.','markersize',1, 'color', color_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_4pi, label_L2_4pi);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% Plot the lat/lon projections
% -------------------------------------------------
if 1+1==1
    figure; hold all
%     title(data_title)
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$',26,'LaTex')
    ylim([-90 90])
    xlim([-180 180])
    
% %     p_high_1pi = plot(data_high_1pi(:,c_lon),data_high_1pi(:,c_lat), '.','markersize',1, 'color', color_high_1pi);
%     p_high_1pi = plot(data_high_1pi(1:200:size(data_high_1pi,1),c_lon),data_high_1pi(1:200:size(data_high_1pi,1),c_lat), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot(data_L1_1yr(:,c_lon),data_L1_1yr(:,c_lat), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_lon),data_L2_1yr(:,c_lat), '.','markersize',1, 'color', color_L2_1yr);
%     
    [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
%     
%     
%     p_L2_2yr   = plot(data_L2_2yr(:,c_lon),data_L2_2yr(:,c_lat), '.','markersize',1, 'color', colors.turq);
%     p_L1_2yr   = plot(data_L1_2yr(:,c_lon),data_L1_2yr(:,c_lat), '.','markersize',1, 'color', colors.drkgrn);
%     
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr, p_L1_2yr, p_L2_2yr],label_L1_1yr, label_L2_1yr, 'L1, 2 years','L2, 2 years');
%     lineh = findobj(objh,'type','line');
%     set(lineh,'linestyle','-','linewidth',3);
    
    
    p_L1_3yr = plot(data_L1_3yr(:,c_lon),data_L1_3yr(:,c_lat), '.','markersize',1, 'color', colors.black);
    p_L1_2yr = plot(data_L1_2yr(:,c_lon),data_L1_2yr(:,c_lat), '.','markersize',1, 'color', colors.drkgrn);
    p_L1_1yr = plot(data_L1_1yr(:,c_lon),data_L1_1yr(:,c_lat), '.','markersize',1, 'color', color_L1_1yr);
    [legh,objh] = legend([p_L1_1yr, p_L1_2yr, p_L1_3yr],label_L1_1yr, 'L1, 2 years','L1, 3 years');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);


%     dataLatLon_L1_2ndYr = data_L1_2yr(data_L1_2yr(:,c_t) >= 31557600/tNorm, [c_lat, c_lon, c_t]);
%     dataLatLon_L2_2ndYr = data_L2_2yr(data_L2_2yr(:,c_t) >= 31557600/tNorm, [c_lat, c_lon, c_t]);
    p_L1_1styr = plot(data_L1_1yr(:,c_lon),data_L1_1yr(:,c_lat), '.','markersize',1, 'color', colors.drkgrn); 
    p_L1_2ndyr = plot(dataLatLon_L1_2ndYr(:,2), dataLatLon_L1_2ndYr(:,1), '.', 'markersize', 1, 'color', color_L1_1yr);
    [legh,objh] = legend([p_L1_1styr, p_L1_2ndyr],'L1, 1st year', 'L1, 2nd year');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    p_L2_1styr = plot(data_L2_1yr(:,c_lon),data_L2_1yr(:,c_lat), '.','markersize',1, 'color', colors.drkgrn); 
    p_L2_2ndyr = plot(dataLatLon_L2_2ndYr(:,2), dataLatLon_L2_2ndYr(:,1), '.', 'markersize', 1, 'color', color_L1_1yr);
    [legh,objh] = legend([p_L2_1styr, p_L2_2ndyr],'L2, 1st year', 'L2, 2nd year');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    
    figure; hold all
    title(data_title)
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$',23,'LaTex')
    ylim([-90 90])
    p_high_1pi = plot(data_high_1pi(1:200:size(data_high_1pi,1),c_lon),data_high_1pi(1:200:size(data_high_1pi,1),c_lat), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot(data_L1_4pi(:,c_lon),data_L1_4pi(:,c_lat), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_4pi(:,c_lon),data_L2_4pi(:,c_lat), '.','markersize',1, 'color', color_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_4pi, label_L2_4pi);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
% dataLatLon_L1_1stYr = data_L1_3yr( data_L1_3yr(:,c_t) < 31557600/tNorm,               [c_lat, c_lon, c_t]);
% dataLatLon_L1_2ndYr = data_L1_3yr( data_L1_3yr(:,c_t) >= 31557600/tNorm,               [c_lat, c_lon, c_t]);
% dataLatLon_L1_2ndYr = dataLatLon_L1_2ndYr(dataLatLon_L1_2ndYr(:,3) < 63115200/tNorm,   :);
% dataLatLon_L1_3rdYr = data_L1_3yr(data_L1_3yr(:,c_t) >= 63115200/tNorm,                [c_lat, c_lon, c_t]);
% p_L1_3rdyr = plot(dataLatLon_L1_1stYr(:,2),dataLatLon_L1_1stYr(:,1), '.','markersize',1, 'color', colors.black);
% p_L1_2ndyr = plot(dataLatLon_L1_2ndYr(:,2),dataLatLon_L1_2ndYr(:,1), '.','markersize',1, 'color', colors.drkgrn);
% % p_L1_1styr = plot(dataLatLon_L1_3rdYr(:,2),dataLatLon_L1_3rdYr(:,1), '.','markersize',1, 'color', color_L1_1yr);
% [legh,objh] = legend([p_L1_1styr, p_L1_2ndyr, p_L1_3rdyr],'L1, 1st year','L1, 2nd year', 'L1, 3rd year');
% lineh = findobj(objh,'type','line');
% set(lineh,'linestyle','-','linewidth',3);


indicies_L1_1styr = data_L1_3yr(:,c_t) < 31557600/tNorm;
indicies_L1_2ndyr = and(data_L1_3yr(:,c_t) >= 31557600/tNorm, data_L1_3yr(:,c_t) < 63115200/tNorm);
indicies_L1_3rdyr = data_L1_3yr(:,c_t) >= 63115200/tNorm;

indicies_L2_1styr = data_L2_3yr(:,c_t) < 31557600/tNorm;
indicies_L2_2ndyr = and(data_L2_3yr(:,c_t) >= 31557600/tNorm, data_L2_3yr(:,c_t) < 63115200/tNorm);
indicies_L2_3rdyr = data_L2_3yr(:,c_t) >= 63115200/tNorm;




figure; hold all
title(data_title)
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$',23,'LaTex')
ylim([-90 90])
xlim([-180 180])

color_high_1pi = [221, 208, 255]./255;
% p_high_1pi = plot(data_high_1pi(1:200:end,c_lon),data_high_1pi(1:200:end,c_lat), '.','markersize',0.01, 'color', color_high_1pi);
p_high_1pi = plot(data_high_1pi(:,c_lon),data_high_1pi(:,c_lat), '.','markersize',0.01, 'color', color_high_1pi);
p_L1_1styr = plot(data_L1_3yr(indicies_L1_1styr,c_lon), data_L1_3yr(indicies_L1_1styr,c_lat),'.','markersize',1,'color',colors.grn2);
p_L1_2ndyr = plot(data_L1_3yr(indicies_L1_2ndyr,c_lon), data_L1_3yr(indicies_L1_2ndyr,c_lat),'.','markersize',1,'color',colors.ltgrn);
p_L1_3rdyr = plot(data_L1_3yr(indicies_L1_3rdyr,c_lon), data_L1_3yr(indicies_L1_3rdyr,c_lat),'.','markersize',1,'color',colors.ylwgrn);
p_L2_1styr = plot(data_L2_3yr(indicies_L2_1styr,c_lon), data_L2_3yr(indicies_L2_1styr,c_lat),'.','markersize',1,'color',colors.blue);
p_L2_2ndyr = plot(data_L2_3yr(indicies_L2_2ndyr,c_lon), data_L2_3yr(indicies_L2_2ndyr,c_lat),'.','markersize',1,'color',colors.blue2);
p_L2_3rdyr = plot(data_L2_3yr(indicies_L2_3rdyr,c_lon), data_L2_3yr(indicies_L2_3rdyr,c_lat),'.','markersize',1,'color',colors.ltblue);
% p_L2_4thyr = plot(data_L2_4yr(indicies_L2_4thyr,c_lon), data_L2_4yr(indicies_L2_4thyr,c_lat),'.','markersize',1,'color',colors.purp);
[legh,objh] = legend([p_high_1pi, p_L1_1styr, p_L1_2ndyr, p_L1_3rdyr, p_L2_1styr, p_L2_2ndyr, p_L2_3rdyr],label_high_1pi,'L1, 1st year', 'L1, 2nd year', 'L1, 3rd year', 'L2, 1st year', 'L2, 2nd year', 'L2, 3rd year');
% [legh,objh] = legend([ p_L1_1styr, p_L1_2ndyr, p_L1_3rdyr, p_L2_1styr, p_L2_2ndyr, p_L2_3rdyr],'L1, 1st year', 'L1, 2nd year', 'L1, 3rd year', 'L2, 1st year', 'L2, 2nd year', 'L2, 3rd year');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);




end


% -------------------------------------------------
%%% Plot Periapses
% -------------------------------------------------
if 1+1==1
    figure; hold all
    axis equal
%     title(data_title)
    PlotBoi3_CR3Bn(23)
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
    xlim([L123(1,1)-secondary.R_n, L123(2,1)+secondary.R_n])
    ylim([-1 1].*0.015)
    view(0,90)
    
    
%     p_L1_1yr = plot3(data_L1_1yr(data_L1_1yr(:,c_isPer) == true,c_x),data_L1_1yr(data_L1_1yr(:,c_isPer) == true,c_y),data_L1_1yr(data_L1_1yr(:,c_isPer) == true,c_z),'.','color',color_L1_1yr,'markersize',1);
%     p_L2_1yr = plot3(data_L2_1yr(data_L2_1yr(:,c_isPer) == true,c_x),data_L2_1yr(data_L2_1yr(:,c_isPer) == true,c_y),data_L2_1yr(data_L2_1yr(:,c_isPer) == true,c_z),'.','color',color_L2_1yr,'markersize',1);
    p_high_1pi = plot3(data_high_1pi(data_high_1pi(:,c_isPer) == true,c_x),data_high_1pi(data_high_1pi(:,c_isPer) == true,c_y),data_high_1pi(data_high_1pi(:,c_isPer) == true,c_z),'.','color',color_high_1pi,'markersize',1);
    
    plot3(L123(1:2,1),zeros(2,1),zeros(2,1),'^','markeredgecolor',colors.black,'markerfacecolor',colors.black,'markersize',10);
    
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr, p_high_1pi], label_L1_1yr, label_L2_1yr, label_high_1pi);
    [legh,objh] = legend([ p_high_1pi], label_high_1pi); zlim([-1 1].*5.7e-3);
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr], label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% Plot Apoapses
% -------------------------------------------------
if 1+1==1
    figure; hold all
    axis equal
%     title(data_title)
    PlotBoi3_CR3Bn(23)
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
    xlim([L123(1,1)-secondary.R_n, L123(2,1)+secondary.R_n])
    ylim([-1 1].*0.015)
    view(0,90)
    
    
    p_L1_1yr = plot3(data_L1_1yr(data_L1_1yr(:,c_isApo) == true,c_x),data_L1_1yr(data_L1_1yr(:,c_isApo) == true,c_y),data_L1_1yr(data_L1_1yr(:,c_isApo) == true,c_z),'.','color',color_L1_1yr,'markersize',1);
    p_L2_1yr = plot3(data_L2_1yr(data_L2_1yr(:,c_isApo) == true,c_x),data_L2_1yr(data_L2_1yr(:,c_isApo) == true,c_y),data_L2_1yr(data_L2_1yr(:,c_isApo) == true,c_z),'.','color',color_L2_1yr,'markersize',1);
%     p_high_1pi = plot3(data_high_1pi(data_high_1pi(:,c_isApo) == true,c_x),data_high_1pi(data_high_1pi(:,c_isApo) == true,c_y),data_high_1pi(data_high_1pi(:,c_isApo) == true,c_z),'.','color',color_high_1pi,'markersize',1);
    
    plot3(L123(1:2,1),zeros(2,1),zeros(2,1),'^','markeredgecolor',colors.black,'markerfacecolor',colors.black,'markersize',10);
    
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr, p_high_1pi], label_L1_1yr, label_L2_1yr, label_high_1pi);
    [legh,objh] = legend([p_L1_1yr, p_L2_1yr], label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% Plot 3D X-Y Crossings
% -------------------------------------------------
if 1+1==1
    figure; hold all
    PlotBoi3('$x_n$','$y_n$','$i$, $^\circ$',23,'LaTex')
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], colors.blue, 0.0)
    
    p1 = plot3(data1_BaCR_zEquals0(:,c_x),data1_BaCR_zEquals0(:,c_y),data1_BaCR_zEquals0(:,c_i),'.','color',data1_color,'markersize',1);
    p2 = plot3(data2_BaCR_zEquals0(:,c_x),data2_BaCR_zEquals0(:,c_y),data2_BaCR_zEquals0(:,c_i),'.','color',data2_color,'markersize',1);
    
    [legh,objh] = legend([p1, p2], data1_label, data2_label);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end


% -------------------------------------------------
%%% |R| at Periapsis Mapping
% -------------------------------------------------
if 1+1==1
    RS_data1_periapsis = rowNorm(data1_BaCR_periapses(:,c_x:c_z) - [1-secondary.MR,0,0]);
    RS_data2_periapsis = rowNorm(data2_BaCR_periapses(:,c_x:c_z) - [1-secondary.MR,0,0]);

    figure; hold all
    title('|R| of Periapse vs latitude')
    PlotBoi2('RS$_{per}$','Latitude, $^\circ$',23,'LaTex')

    p1  = plot(RS_data1_periapsis,data1_BaCR_periapses(:,c_lat),'.','color',data1_color,'markersize',8);
    p2 = plot(RS_data2_periapsis,data2_BaCR_periapses(:,c_lat),'.','color',data2_color,'markersize',8);

    [legh,objh] = legend([p1, p2],data1_label,data2_label);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% |R2| vs Latitude ... and ... z vs Latitude
% -------------------------------------------------
if 1+1==1
%     RS_data1 = rowNorm(data1_BaCR(:,c_x:c_z) - [1-secondary.MR,0,0]);
%     RS_data2 = rowNorm(data2_BaCR(:,c_x:c_z) - [1-secondary.MR,0,0]);
%     RS_data4 = rowNorm(data4_BaCR(:,c_x:c_z) - [1-secondary.MR,0,0]);
%     RS_data5 = rowNorm(data5_BaCR(:,c_x:c_z) - [1-secondary.MR,0,0]);
% 
%     figure; hold all
%     title('|R| vs latitude')
%     PlotBoi2('RS$_{per}$','Latitude, $^\circ$',23,'LaTex')
%     xlim([0, (L123(2,1)-(1-secondary.MR))])
% 
%     p5 = plot(RS_data5,data5_BaCR(:,c_lat),'.','color',data5_color,'markersize',2);
%     p1 = plot(RS_data1,data1_BaCR(:,c_lat),'.','color',data1_color,'markersize',2);
%     p4 = plot(RS_data4,data4_BaCR(:,c_lat),'.','color',data4_color,'markersize',2);
%     p2 = plot(RS_data2,data2_BaCR(:,c_lat),'.','color',data2_color,'markersize',2);
%     pS = plot([secondary.R_n, secondary.R_n], [-90, 90], 'k', 'linewidth',1);
%     
%     [legh,objh] = legend([p1, p2, p4,p5, pS], data1_label ,data2_label, data4_label, data5_label, 'Surface');
%     lineh = findobj(objh,'type','line');
%     set(lineh,'linestyle','-','linewidth',3);

    Rs_L1_1yr = rowNorm(data_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L2_1yr = rowNorm(data_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    
    figure; hold all
%     title('r_2 vs latitude')
    PlotBoi2('$r_2$','Latitude, $^\circ$',23,'LaTex')
    xlim([0, (L123(2,1)-(1-secondary.MR))])
    
%     p_high_1pi = plot(Rs_high_1pi, data_high_1pi(:,c_lat), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr = plot(Rs_L1_1yr, data_L1_1yr(:,c_lat), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr = plot(Rs_L2_1yr, data_L2_1yr(:,c_lat), '.','markersize',1, 'color', color_L2_1yr);
    
    [legh,objh] = legend([p_L1_1yr, p_L2_1yr], label_L1_1yr ,label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    figure; hold all
    PlotBoi2('$|z_n|$','$|$Latitude$|$, $^\circ$',23,'LaTex')
    p_L1_1yr = plot(abs(data_L1_1yr(:,c_z)), abs(data_L1_1yr(:,c_lat)), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr = plot(abs(data_L2_1yr(:,c_z)), abs(data_L2_1yr(:,c_lat)), '.','markersize',1, 'color', color_L2_1yr);
    p_poles  = plot([secondary.R_n, secondary.R_n],[0 max(abs(data_L1_1yr(:,c_lat)))],'k','linewidth',1);
    [legh,objh] = legend([p_L1_1yr, p_L2_1yr, p_poles], label_L1_1yr ,label_L2_1yr, 'R_{Europa}');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    ylim([0 max(abs(data_L1_1yr(:,c_lat)))])
end

% -------------------------------------------------
%%% i vs a vs H
% -------------------------------------------------
if 1+1==1
    figure; hold all
    PlotBoi3('$i$','$a$','$|h|$',23,'LaTex')

    p1 = plot3(data1_BaCR_zEquals0(:,c_i),data1_BaCR_zEquals0(:,c_a),data1_BaCR_zEquals0(:,c_H),'.','color',data1_color,'markersize',1);
    p2 = plot3(data2_BaCR_zEquals0(:,c_i),data2_BaCR_zEquals0(:,c_a),data2_BaCR_zEquals0(:,c_H),'.','color',data2_color,'markersize',1);
    
    [legh,objh] = legend([p1, p2], data1_label ,data2_label);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% norm(Hx,Hy) vs H .... and vs i .... and vs R
% -------------------------------------------------
if 1+1==1
    %%% Calculate HxHx data
    HxHy_high_1pi = rowNorm(data_high_1pi(:,c_hx:c_hy));
    HxHy_L1_1yr   = rowNorm(data_L1_1yr(:,c_hx:c_hy));
    HxHy_L2_1yr   = rowNorm(data_L2_1yr(:,c_hx:c_hy));

    figure; hold all
    title(data_title)
    PlotBoi2('$|H|$','$|H_xH_y|$',23,'LaTex')
    
    p_L1_1yr   = plot(data_L1_1yr(:,c_H),HxHy_L1_1yr, '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_H),HxHy_L2_1yr, '.','markersize',1, 'color', color_L2_1yr);
    p_high_1pi = plot(data_high_1pi(:,c_H),HxHy_high_1pi, '.','markersize',1, 'color', color_high_1pi);
    
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    

    % ----------------
    
    figure; hold all
    title(data_title)
    PlotBoi2('$i$','$|H_xH_y|$',23,'LaTex')
    
    p_L1_1yr   = plot(data_L1_1yr(:,c_i),HxHy_L1_1yr, '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_i),HxHy_L2_1yr, '.','markersize',1, 'color', color_L2_1yr);
    p_high_1pi = plot(data_high_1pi(:,c_i),HxHy_high_1pi, '.','markersize',1, 'color', color_high_1pi);

    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    % ----------------
    if ~exist('Rs_L1_1yr','var')
        Rs_L1_1yr = rowNorm(data_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    end
    if ~exist('Rs_L2_1yr','var')
        Rs_L2_1yr = rowNorm(data_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    end
    if ~exist('Rs_high_1pi','var')
        Rs_high_1pi = rowNorm(data_high_1pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    end
    
    figure; hold all
    title(data_title)
    PlotBoi2('$R$','$|H_xH_y|$',23,'LaTex')
    
    p_L1_1yr   = plot(Rs_L1_1yr,HxHy_L1_1yr, '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(Rs_L2_1yr,HxHy_L2_1yr, '.','markersize',1, 'color', color_L2_1yr);
    p_high_1pi = plot(Rs_high_1pi,HxHy_high_1pi, '.','markersize',1, 'color', color_high_1pi);

    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
end



% -------------------------------------------------
%%% Delaunay variables
% -------------------------------------------------
if 1+1==1
    % L & G
    figure; hold all
    title(data_title)
    PlotBoi2('$L$','$G$',23,'LaTex')
    p_high_1pi = plot(data_high_1pi(:,c_L_del), data_high_1pi(:,c_G_del), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot(data_L1_1yr(:,c_L_del), data_L1_1yr(:,c_G_del), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_L_del), data_L2_1yr(:,c_G_del), '.','markersize',1, 'color', color_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    % L & H
    figure; hold all
    title(data_title)
    PlotBoi2('$L$','$H$',23,'LaTex')
    p_high_1pi = plot(data_high_1pi(:,c_L_del), data_high_1pi(:,c_H_del), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot(data_L1_1yr(:,c_L_del), data_L1_1yr(:,c_H_del), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_L_del), data_L2_1yr(:,c_H_del), '.','markersize',1, 'color', color_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    % G & H
    figure; hold all
    title(data_title)
    PlotBoi2('$H$','$G$',23,'LaTex')
    p_high_1pi = plot(data_high_1pi(:,c_H_del), data_high_1pi(:,c_G_del), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot(data_L1_1yr(:,c_H_del), data_L1_1yr(:,c_G_del), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_H_del), data_L2_1yr(:,c_G_del), '.','markersize',1, 'color', color_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
end

% -------------------------------------------------
%%% Latitude vs inclination
% -------------------------------------------------
if 1+1==1
    figure; hold all
    title(data_title)
    PlotBoi2('$i_{SCI}$, $^\circ$','Latitude, $^\circ$',23,'LaTex')
    
    p_high_1pi = plot(data_high_1pi(:,c_i),data_high_1pi(:,c_lat), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot(data_L1_1yr(:,c_i),data_L1_1yr(:,c_lat), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_i),data_L2_1yr(:,c_lat), '.','markersize',1, 'color', color_L2_1yr);
    
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
%     data_L1_4pi
%     data_L1_1yr
%     data_L2_4pi
%     data_L2_4pi_impact
%     data_L2_1yr
%     data_high_1pi
end


% -------------------------------------------------
%%% i_SCR vs R (takes ~12 minutes to create all these vectors)
% -------------------------------------------------
if 1+1==1
    if ~exist('Rs_L1_1yr','var')
        Rs_L1_1yr = rowNorm(data_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    end
    if ~exist('Rs_L2_1yr','var')
        Rs_L2_1yr = rowNorm(data_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    end
    if ~exist('Rs_high_1pi','var')
        Rs_high_1pi = rowNorm(data_high_1pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    end
    
    if ~exist('is_SCR_L1_1yr','var')
        is_SCR_L1_yr = NaN(size(data_L1_1yr,1),1);
        for kk = 1:length(is_SCR_L1_yr)
            [a,e,i,raan,w,ta] = ECI2OE(data_L1_1yr(kk,c_x:c_z)-[1-secondary.MR,0,0],data_L1_1yr(kk,c_xd:c_zd),secondary.MR);
            is_SCR_L1_yr(kk) = i*180/pi;
        end
    end
    if ~exist('is_SCR_L2_1yr','var')
        is_SCR_L2_yr = NaN(size(data_L2_1yr,1),1);
        for kk = 1:length(is_SCR_L2_yr)
            [a,e,i,raan,w,ta] = ECI2OE(data_L2_1yr(kk,c_x:c_z)-[1-secondary.MR,0,0],data_L2_1yr(kk,c_xd:c_zd),secondary.MR);
            is_SCR_L2_yr(kk) = i*180/pi;
        end
    end
    if ~exist('is_SCR_high_1pi','var')
        is_SCR_high_1pi = NaN(size(data_high_1pi,1),1);
        for kk = 1:length(is_SCR_high_1pi)
            [a,e,i,raan,w,ta] = ECI2OE(data_high_1pi(kk,c_x:c_z)-[1-secondary.MR,0,0],data_high_1pi(kk,c_xd:c_zd),secondary.MR);
            is_SCR_high_1pi(kk) = i*180/pi;
        end
    end
    
    
    figure; hold all
%     title(data_title)
    
    xlim([0, (L123(2,1)-(1-secondary.MR))*2])
    
    PlotBoi2('R$_n$','$i_{SCR}$, $^\circ$',23,'LaTex')
    p_L1_1yr   = plot(Rs_L1_1yr, is_SCR_L1_yr, '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(Rs_L2_1yr, is_SCR_L2_yr, '.','markersize',1, 'color', color_L2_1yr);
    p_high_1pi = plot(Rs_high_1pi, is_SCR_high_1pi, '.','markersize',1, 'color', color_high_1pi);
    pS = plot([secondary.R_n, secondary.R_n], [0, 180], 'k', 'linewidth',1);
    
%     PlotBoi3('R$_n$','$i_{SCR}$, $^\circ$', 'Longitude, $^\circ$' ,23,'LaTex')
%     p_L1_1yr   = plot3(Rs_L1_1yr, is_SCR_L1_yr, data_L1_1yr(:,c_lon), '.','markersize',1, 'color', color_L1_1yr);
%     p_L2_1yr   = plot3(Rs_L2_1yr, is_SCR_L2_yr, data_L2_1yr(:,c_lon), '.','markersize',1, 'color', color_L2_1yr);
%     p_high_1pi = plot3(Rs_high_1pi, is_SCR_high_1pi, data_high_1pi(:,c_lon), '.','markersize',1, 'color', color_high_1pi);
%     pS = plot([secondary.R_n, secondary.R_n], [0, 180], 'k', 'linewidth',1);
%     
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, pS], label_high_1pi  ,label_L1_1yr, label_L2_1yr,'Surface');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    % Longitude vs R
    figure; hold all
%     title(data_title)
    xlim([0, (L123(2,1)-(1-secondary.MR))])
    ylim([-1 1].*180)
    PlotBoi2('R$_n$','Longitude, $^\circ$',23,'LaTex')
    p_L1_1yr   = plot(Rs_L1_1yr(:), data_L1_1yr(:,c_lon), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(Rs_L2_1yr(1:20:end), data_L2_1yr(1:20:end,c_lon), '.','markersize',1, 'color', color_L2_1yr);
    p_high_1pi = plot(Rs_high_1pi(1:200:end), data_high_1pi(1:200:end,c_lon), '.','markersize',1, 'color', color_high_1pi);
    pS = plot([secondary.R_n, secondary.R_n], [-180, 180], 'k', 'linewidth',1);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, pS], label_high_1pi  ,label_L1_1yr, label_L2_1yr,'Surface');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    % Latitude vs R
    figure; hold all
%     title('r_2 vs latitude')
    PlotBoi2('$R_n$','Latitude, $^\circ$',23,'LaTex')
    xlim([0, (L123(2,1)-(1-secondary.MR))])
    ylim([-1 1].*90)
    p_high_1pi = plot(Rs_high_1pi(1:200:end), data_high_1pi(1:200:end,c_lat), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr = plot(Rs_L1_1yr, data_L1_1yr(:,c_lat), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr = plot(Rs_L2_1yr, data_L2_1yr(:,c_lat), '.','markersize',1, 'color', color_L2_1yr);
    pS = plot([secondary.R_n, secondary.R_n], [-90, 90], 'k', 'linewidth',1);
    [legh,objh] = legend([p_L1_1yr, p_L2_1yr, pS], label_L1_1yr ,label_L2_1yr,'Surface');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    
    
    figure; hold all
%     title(data_title)
    
    xlim([0, (L123(2,1)-(1-secondary.MR))*2])
    
    PlotBoi2('R$_n$','$i_{SCI}$, $^\circ$',23,'LaTex')
    p_L1_1yr   = plot(Rs_L1_1yr, data_L1_1yr(:,c_i), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(Rs_L2_1yr, data_L2_1yr(:,c_i), '.','markersize',1, 'color', color_L2_1yr);
    p_high_1pi = plot(Rs_high_1pi, data_high_1pi(:,c_i), '.','markersize',1, 'color', color_high_1pi);
    pS = plot([secondary.R_n, secondary.R_n], [0, 180], 'k', 'linewidth',1);
    
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, pS], label_high_1pi  ,label_L1_1yr, label_L2_1yr,'Surface');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
end


% -------------------------------------------------
%%% Latitude over time
% -------------------------------------------------
if 1+1 == 1
    figure; hold all
    title(data_title)
    PlotBoi2('$t$, years','Latitude, $^\circ$',23,'LaTex')

    p_L1_3yr = plot(data_L1_3yr(:,c_t).*tNorm/(86400*365.25), data_L1_3yr(:,c_lat),'.','markersize',1,'color',color_L1_1yr);
    p_L2_4yr = plot(data_L2_4yr(:,c_t).*tNorm/(86400*365.25), data_L2_4yr(:,c_lat),'.','markersize',1,'color',color_L2_1yr);

    [legh,objh] = legend([p_L1_3yr, p_L2_4yr], 'L1, 3 years'  ,'L2, 4 years');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end


% -------------------------------------------------
%%% z-energy analysis
% -------------------------------------------------
if 1+1 == 1
    if ~exist('zKEnergy_high_1pi','var')
        zKEnergy_high_1pi = NaN(size(data_high_1pi,1),1);
        for kk = 1:length(zKEnergy_high_1pi)
            zKEnergy_high_1pi(kk) = (data_high_1pi(kk,c_zd)^2)/2;
        end
    end
    if ~exist('zKEnergy_L1_1yr','var')
        zKEnergy_L1_1yr = NaN(size(data_L1_1yr,1),1);
        for kk = 1:length(zKEnergy_L1_1yr)
            zKEnergy_L1_1yr(kk) = (data_L1_1yr(kk,c_zd)^2)/2;
        end
    end
    if ~exist('zKEnergy_L2_1yr','var')
        zKEnergy_L2_1yr = NaN(size(data_L2_1yr,1),1);
        for kk = 1:length(zKEnergy_L2_1yr)
            zKEnergy_L2_1yr(kk) = (data_L2_1yr(kk,c_zd)^2)/2;
        end
    end
    
    if ~exist('is_SCR_L1_1yr','var')
        is_SCR_L1_yr = NaN(size(data_L1_1yr,1),1);
        for kk = 1:length(is_SCR_L1_yr)
            [a,e,i,raan,w,ta] = ECI2OE(data_L1_1yr(kk,c_x:c_z)-[1-secondary.MR,0,0],data_L1_1yr(kk,c_xd:c_zd),secondary.MR);
            is_SCR_L1_yr(kk) = i*180/pi;
        end
    end
    if ~exist('is_SCR_L2_1yr','var')
        is_SCR_L2_yr = NaN(size(data_L2_1yr,1),1);
        for kk = 1:length(is_SCR_L2_yr)
            [a,e,i,raan,w,ta] = ECI2OE(data_L2_1yr(kk,c_x:c_z)-[1-secondary.MR,0,0],data_L2_1yr(kk,c_xd:c_zd),secondary.MR);
            is_SCR_L2_yr(kk) = i*180/pi;
        end
    end
    if ~exist('is_SCR_high_1pi','var')
        is_SCR_high_1pi = NaN(size(data_high_1pi,1),1);
        for kk = 1:length(is_SCR_high_1pi)
            [a,e,i,raan,w,ta] = ECI2OE(data_high_1pi(kk,c_x:c_z)-[1-secondary.MR,0,0],data_high_1pi(kk,c_xd:c_zd),secondary.MR);
            is_SCR_high_1pi(kk) = i*180/pi;
        end
    end
    
figure; hold all 
PlotBoi2('$\dot{z}^2/2$','Latitude, $^\circ$',23,'LaTex')
% title(data_title)
p_high =  plot(zKEnergy_high_1pi(1:200:end),data_high_1pi(1:200:end,c_lat),'.','color',color_high_1pi);
p_L1 = plot(zKEnergy_L1_1yr,data_L1_1yr(:,c_lat),'.','color',color_L1_1yr);
p_L2 = plot(zKEnergy_L2_1yr,data_L2_1yr(:,c_lat),'.','color',color_L2_1yr);
% [legh,objh] = legend([p_high, p_L1, p_L2], label_high_1pi, label_L1_1yr, label_L2_1yr);
[legh,objh] = legend([p_high, p_L1, p_L2], label_high_1pi, label_L1_1yr, label_L2_1yr);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);
xlim([0 0.035])
ylim([-1 1].*90)


figure; hold all 
PlotBoi2('$\dot{z}^2/2$','$i_{SCI}$, $^\circ$',23,'LaTex')
title(data_title)
p_high = plot(zKEnergy_high_1pi,data_high_1pi(:,c_i),'.','color',color_high_1pi);
p_L1 = plot(zKEnergy_L1_1yr,data_L1_1yr(:,c_i),'.','color',color_L1_4pi);
p_L2 = plot(zKEnergy_L2_1yr,data_L2_1yr(:,c_i),'.','color',color_L2_4pi);
[legh,objh] = legend([p_high, p_L1, p_L2], label_high_1pi, label_L1_1yr, label_L2_1yr);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);
xlim([0 0.05])

figure; hold all 
PlotBoi2('$\dot{z}^2/2$','$i_{SCR}$, $^\circ$',23,'LaTex')
p_high = plot(zKEnergy_high_1pi, is_SCR_high_1pi,'.','color',color_high_1pi);
p_L1 = plot(zKEnergy_L1_1yr, is_SCR_L1_yr,'.','color',color_L1_4pi);
p_L2 = plot(zKEnergy_L2_1yr, is_SCR_L2_yr,'.','color',color_L2_4pi);
[legh,objh] = legend([p_high, p_L1, p_L2], label_high_1pi, label_L1_1yr, label_L2_1yr);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);

end


% -------------------------------------------------
%%% Double checking things
% -------------------------------------------------
if 2 + 2 == 2
% % 
% % 
high_range = 1:1000:size(data_high_1pi,1);
L1_range = 1:1000:size(data_L1_1yr,1);
L2_range = 1:100:size(data_L2_1yr,1);
% % 
% c_prm = c_a;
% % 
figure; hold all 
% xlim([0 2*(L123(2,1)-(1-secondary.MR))])
PlotBoi2('Data1','Data2',23,'LaTex')
% % % p_high = plot(Rs_high_1pi(high_range), data_high_1pi(high_range,c_prm),'.','color',color_high_1pi);
% % % p_L1   = plot(Rs_L1_1yr(L1_range), data_L1_1yr(L1_range,c_prm),'.','color',color_L1_1yr);
% % % p_L2   = plot(Rs_L2_1yr(L2_range), data_L2_1yr(L2_range,c_prm),'.','color',color_L2_1yr);
% % 
% % % p_high = plot(data_high_1pi(high_range, c_lat), zdOverV_high(high_range),'.','color',color_high_1pi);
% % % p_L1   = plot(data_L1_1yr(L1_range, c_lat), zdOverV_L1(L1_range),'.','color',color_L1_1yr);
% % % p_L2   = plot(data_L2_1yr(L2_range, c_lat), zdOverV_L2(L2_range),'.','color',color_L2_1yr);
% % 
% % p_high = plot(data_high_1pi(:, c_lat), HxyOverH_high(:),'.','color',color_high_1pi);
% % p_L1   = plot(data_L1_1yr(:, c_lat), HxyOverH_L1(:),'.','color',color_L1_1yr);
% % p_L2   = plot(data_L2_1yr(:, c_lat), HxyOverH_L2(:),'.','color',color_L2_1yr);
% % 
% % p_high = plot(Rs_high_1pi, HxyOverH_high(:),'.','color',color_high_1pi);
% % p_L1   = plot(Rs_L1_1yr, HxyOverH_L1(:),'.','color',color_L1_1yr);
% % p_L2   = plot(Rs_L2_1yr, HxyOverH_L2(:),'.','color',color_L2_1yr);
% 
% p_high = plot(Rs_high_1pi, HxHy_high_1pi(:),'.','color',color_high_1pi);
% p_L1   = plot(Rs_L1_1yr, HxHy_L1_1yr(:),'.','color',color_L1_1yr);
% p_L2   = plot(Rs_L2_1yr, HxHy_L2_1yr(:),'.','color',color_L2_1yr);
% % 
% % % p_high = plot(Rs_high_1pi(:), data_high_1pi(:,c_prm),'.','color',color_high_1pi);
% % % p_L1   = plot(Rs_L1_1yr(:), data_L1_1yr(:,c_prm),'.','color',color_L1_1yr);
% % % p_L2   = plot(Rs_L2_1yr(:), data_L2_1yr(:,c_prm),'.','color',color_L2_1yr);
% % 

p_high = plot(data_high_1pi(high_range, c_z),     data_high_1pi(high_range, c_zd),'.','color',color_high_1pi);
p_L1   = plot(data_L1_1yr(L1_range,c_z),          data_L1_1yr(L1_range,c_zd),'.','color',color_L1_1yr);
p_L2   = plot(data_L2_1yr(L2_range,c_z),          data_L2_1yr(L2_range,c_zd),'.','color',color_L2_1yr);

[legh,objh] = legend([p_high, p_L1, p_L2], label_high_1pi, label_L1_1yr, label_L2_1yr);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);


% % % 
% % % % zd / V
% % % zdOverV_high = zeros(size(data_high_1pi,1),1);
% % % zdOverV_L1   = zeros(size(data_L1_1yr,1),1);
% % % zdOverV_L2   = zeros(size(data_L2_1yr,1),1);
% % % 
% % % for kk = 1:size(data_high_1pi,1)
% % %     zdOverV_high(kk) = data_high_1pi(kk, c_zd) / norm(data_high_1pi(kk,c_xd:c_zd));
% % % end
% % % 
% % % for kk = 1:size(data_L1_1yr,1)
% % %     zdOverV_L1(kk) = data_L1_1yr(kk, c_zd) / norm(data_L1_1yr(kk,c_xd:c_zd));
% % % end
% % % 
% % % for kk = 1:size(data_L2_1yr,1)
% % %     zdOverV_L2(kk) = data_L2_1yr(kk, c_zd) / norm(data_L2_1yr(kk,c_xd:c_zd));
% % % end
% % % 
% % % 
% % Hxy / H
% HxyOverH_high = zeros(size(data_high_1pi,1),1);
% HxyOverH_L1   = zeros(size(data_L1_1yr,1),1);
% HxyOverH_L2   = zeros(size(data_L2_1yr,1),1);
% 
% for kk = 1:size(data_high_1pi,1)
%     HxyOverH_high(kk) = norm(data_high_1pi(kk,c_hx:c_hy)) / norm(data_high_1pi(kk,c_hx:c_hz));
% end
% 
% for kk = 1:size(data_L1_1yr,1)
%     HxyOverH_L1(kk) = norm(data_L1_1yr(kk,c_hx:c_hy)) / norm(data_L1_1yr(kk,c_hx:c_hz));
% end
% 
% for kk = 1:size(data_L2_1yr,1)
%     HxyOverH_L2(kk) = norm(data_L2_1yr(kk,c_hx:c_hy)) / norm(data_L2_1yr(kk,c_hx:c_hz));
% end





zEnergy_high = zeros(size(data_high_1pi,1),1);
zEnergy_L1   = zeros(size(data_L1_1yr,1),1);
zEnergy_L2   = zeros(size(data_L2_1yr,1),1);

for kk = 1:size(data_high_1pi,1)
    zEnergy_high(kk) = (data_high_1pi(kk,c_zd)^2)/2 - (secondary.MR)/data_high_1pi(kk,c_z);
end

for kk = 1:size(data_L1_1yr,1)
    zEnergy_L1(kk) = (data_L1_1yr(kk,c_zd)^2)/2 - (secondary.MR)/data_L1_1yr(kk,c_z);
end

for kk = 1:size(data_L2_1yr,1)
    zEnergy_L2(kk) = (data_L2_1yr(kk,c_zd)^2)/2 - (secondary.MR)/data_L2_1yr(kk,c_z);
end


% zEnergyRatio_high = zeros(size(data_high_1pi,1),1);
% zEnergyRatio_L1   = zeros(size(data_L1_1yr,1),1);
% zEnergyRatio_L2   = zeros(size(data_L2_1yr,1),1);
% 
% for kk = 1:size(data_high_1pi,1)
% %     zEnergyRatio_high(kk) = ((data_high_1pi(kk,c_zd)^2)/2 - (secondary.MR)/data_high_1pi(kk,c_z)) / ...
% %                             ((norm(data_high_1pi(kk,c_xd:c_zd))^2)/2 - (secondary.MR)/Rs_high_1pi(kk));
%     zEnergyRatio_high(kk) = ((data_high_1pi(kk,c_zd)^2)/2 - (secondary.MR)/Rs_high_1pi(kk)) / ...
%                             ((norm(data_high_1pi(kk,c_xd:c_zd))^2)/2 - (secondary.MR)/Rs_high_1pi(kk));
% end
% 
% for kk = 1:size(data_L1_1yr,1)
% %     zEnergyRatio_L1(kk) = ((data_L1_1yr(kk,c_zd)^2)/2 - (secondary.MR)/data_L1_1yr(kk,c_z)) / ....
% %                           ((norm(data_L1_1yr(kk,c_xd:c_zd))^2)/2 - (secondary.MR)/Rs_L1_1yr(kk));
%     zEnergyRatio_L1(kk) = ((data_L1_1yr(kk,c_zd)^2)/2 - (secondary.MR)/Rs_L1_1yr(kk)) / ....
%                           ((norm(data_L1_1yr(kk,c_xd:c_zd))^2)/2 - (secondary.MR)/Rs_L1_1yr(kk));
% end
% 
% for kk = 1:size(data_L2_1yr,1)
% %     zEnergyRatio_L2(kk) = ((data_L2_1yr(kk,c_zd)^2)/2 - (secondary.MR)/data_L2_1yr(kk,c_z)) / ...
% %                           ((norm(data_L2_1yr(kk,c_xd:c_zd))^2)/2 - (secondary.MR)/Rs_L2_1yr(kk));
%     zEnergyRatio_L2(kk) = ((data_L2_1yr(kk,c_zd)^2)/2 - (secondary.MR)/Rs_L2_1yr(kk)) / ...
%                           ((norm(data_L2_1yr(kk,c_xd:c_zd))^2)/2 - (secondary.MR)/Rs_L2_1yr(kk));
% end


% R vs z
% figure; hold all
% xlim([0 2*(L123(2,1)-(1-secondary.MR))])
% PlotBoi2('$R_n$','$z_n$',23,'LaTex')
% p_L1   = plot(Rs_L1_1yr,abs(data_L1_1yr(:,c_z)),'.','color',color_L1_1yr);
% p_L2   = plot(Rs_L2_1yr, abs(data_L2_1yr(:,c_z)),'.','color',color_L2_1yr);
% p_L2_line = plot([0.020486982866603 0.020486982866603],[0, 6e-3],'k');
% p_L1_line = plot([0.020210923466789 0.020210923466789],[0, 6e-3],'--k');
% [legh,objh] = legend([p_L1, p_L2, p_L1_line, p_L2_line], label_L1_1yr, label_L2_1yr,'R_{L1}','R_{L2}');
% lineh = findobj(objh,'type','line');
% set(lineh(1:4),'linestyle','-','linewidth',3);

% R vs zEnergyRatio
% figure; hold all
% xlim([0 2*(L123(2,1)-(1-secondary.MR))])
% PlotBoi2('$R_n$','$z$ Energy Ratio',23,'LaTex')
% p_high   = plot(Rs_high_1pi,zEnergyRatio_high,'.','color',color_high_1pi);
% p_L1   = plot(Rs_L1_1yr,zEnergyRatio_L1,'.','color',color_L1_1yr);
% p_L2   = plot(Rs_L2_1yr, zEnergyRatio_L2,'.','color',color_L2_1yr);
% p_L2_line = plot([0.020486982866603 0.020486982866603],[1, -1],'k');
% p_L1_line = plot([0.020210923466789 0.020210923466789],[1, -1],'--k');
% [legh,objh] = legend([p_high, p_L1, p_L2, p_L1_line, p_L2_line], label_high_1pi, label_L1_1yr, label_L2_1yr,'R_{L1}','R_{L2}');
% lineh = findobj(objh,'type','line');
% set(lineh(1:6),'linestyle','-','linewidth',3);


% % R vs z vs lat
% figure; hold all
% % xlim([0 2*(L123(2,1)-(1-secondary.MR))])
% PlotBoi3('Longitude, $^\circ$', 'Latitude, $^\circ$','$H_z$',23,'LaTex')
% p_high = plot3(data_high_1pi(high_range,c_lon), data_high_1pi(high_range,c_lat), data_high_1pi(high_range, c_hz),'.','color',color_high_1pi);
% p_L1   = plot3(data_L1_1yr(L1_range,c_lon), data_L1_1yr(L1_range,c_lat), data_L1_1yr(L1_range, c_hz),'.','color',color_L1_1yr);
% p_L2   = plot3(data_L2_1yr(L2_range,c_lon), data_L2_1yr(L2_range,c_lat), data_L2_1yr(L2_range, c_hz),'.','color',color_L2_1yr);
% % p_high = plot3(data_high_1pi(high_range,c_lon), data_high_1pi(high_range,c_lat), HxyOverH_high(high_range),'.','color',color_high_1pi);
% % p_L1   = plot3(data_L1_1yr(L1_range,c_lon), data_L1_1yr(L1_range,c_lat), HxyOverH_L1(L1_range),'.','color',color_L1_1yr);
% % p_L2   = plot3(data_L2_1yr(L2_range,c_lon), data_L2_1yr(L2_range,c_lat), HxyOverH_L2(L2_range),'.','color',color_L2_1yr);
% 
% [legh,objh] = legend([p_high, p_L1, p_L2], label_high_1pi, label_L1_1yr, label_L2_1yr);
% lineh = findobj(objh,'type','line');
% set(lineh,'linestyle','-','linewidth',3);




if 1+1 == 1
%     Rmax = (L123(2,1) - (1-secondary.MR))*1.1;
    Rmax = max(Rs_high_1pi);
    
%     nearIndices_high = Rs_high_1pi<Rmax;
    nearIndices_L1   = Rs_L1_1yr<Rmax;
    nearIndices_L2   = Rs_L2_1yr<Rmax;
    
%     nearData_high = data_high_1pi(nearIndices_high,:);
    nearData_L1 = data_L1_1yr(nearIndices_L1,:);
    nearData_L2 = data_L2_1yr(nearIndices_L2,:);
    
    
    nearRange_L1 = 1:1000:size(nearData_L1,1);
    nearRange_L2 = 1:1000:size(nearData_L2,1);
    
    HxyOverH_near_L1   = zeros(size(nearData_L1,1),1);
    HxyOverH_near_L2   = zeros(size(nearData_L2,1),1);
    for kk = 1:size(nearData_L1,1)
        HxyOverH_near_L1(kk) = norm(nearData_L1(kk,c_hx:c_hy)) / norm(nearData_L1(kk,c_hx:c_hz));
    end
    for kk = 1:size(nearData_L2,1)
        HxyOverH_near_L2(kk) = norm(nearData_L2(kk,c_hx:c_hy)) / norm(nearData_L2(kk,c_hx:c_hz));
    end
    
    
    
    
    
    
    zKEnergyRatio_high = zeros(size(data_high_1pi,1),1);
    zKEnergyRatio_L1   = zeros(size(data_L1_1yr,1),1);
    zKEnergyRatio_L2   = zeros(size(nearData_L2,1),1);
    
    for kk = 1:size(data_high_1pi,1)
        zKEnergyRatio_high(kk) = ((data_high_1pi(kk,c_zd)^2)/2) / ...
                                ((norm(data_high_1pi(kk,c_xd:c_zd))^2)/2 );
    end
    
    for kk = 1:size(nearData_L1,1)
        zKEnergyRatio_L1(kk) = ((nearData_L1(kk,c_zd)^2)/2 ) / ....
                              ((norm(nearData_L1(kk,c_xd:c_zd))^2)/2 );
    end
    
    for kk = 1:size(nearData_L2,1)
        zKEnergyRatio_L2(kk) = ((nearData_L2(kk,c_zd)^2)/2 ) / ...
                              ((norm(nearData_L2(kk,c_xd:c_zd))^2)/2);
    end
    
    
    
    
    
    
    
    % R vs z vs lat
    figure; hold all
    % xlim([0 2*(L123(2,1)-(1-secondary.MR))])
    PlotBoi3('Longitude, $^\circ$', 'Latitude, $^\circ$','$Data3$',23,'LaTex')
%     p_high = plot3(data_high_1pi(high_range,c_lon), data_high_1pi(high_range,c_lat), data_high_1pi(high_range, c_i),'.','color',color_high_1pi);
%     p_L1   = plot3(nearData_L1(nearRange_L1,c_lon), nearData_L1(nearRange_L1,c_lat), nearData_L1(nearRange_L1, c_i),'.','color',color_L1_1yr);
%     p_L2   = plot3(nearData_L2(nearRange_L2,c_lon), nearData_L2(nearRange_L2,c_lat), nearData_L2(nearRange_L2, c_i),'.','color',color_L2_1yr);
    p_high = plot3(data_high_1pi(high_range,c_lon), data_high_1pi(high_range,c_lat), zKEnergyRatio_high(high_range),'.','color',color_high_1pi);
    p_L1   = plot3(nearData_L1(nearRange_L1,c_lon), nearData_L1(nearRange_L1,c_lat), zKEnergyRatio_L1(nearRange_L1),'.','color',color_L1_1yr);
    p_L2   = plot3(nearData_L2(nearRange_L2,c_lon), nearData_L2(nearRange_L2,c_lat), zKEnergyRatio_L2(nearRange_L2),'.','color',color_L2_1yr);

    [legh,objh] = legend([p_high, p_L1, p_L2], label_high_1pi, label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    % R vs z
    
        
    
%     PlotBoi3('R$_n$','$i_{SCR}$, $^\circ$', 'Longitude, $^\circ$' ,23,'LaTex')
%     p_L1_1yr   = plot3(Rs_L1_1yr, is_SCR_L1_yr, data_L1_1yr(:,c_lon), '.','markersize',1, 'color', color_L1_1yr);
%     p_L2_1yr   = plot3(Rs_L2_1yr, is_SCR_L2_yr, data_L2_1yr(:,c_lon), '.','markersize',1, 'color', color_L2_1yr);
%     p_high_1pi = plot3(Rs_high_1pi, is_SCR_high_1pi, data_high_1pi(:,c_lon), '.','markersize',1, 'color', color_high_1pi);
%     pS = plot([secondary.R_n, secondary.R_n], [0, 180], 'k', 'linewidth',1);
%     
    
    
    figure; hold all
    xlim([0 1*(L123(2,1)-(1-secondary.MR))])
    PlotBoi2('$R_n$','$z_n$',23,'LaTex')
    p_L1_1yr   = plot(Rs_L1_1yr, data_L1_1yr(:,c_z),'.','color',color_L1_1yr);
    p_L2_1yr   = plot(Rs_L2_1yr, data_L2_1yr(:,c_z),'.','color',color_L2_1yr);
    pS = plot([secondary.R_n, secondary.R_n], [-1, 1].*6e-3, 'k', 'linewidth',1);
    [legh,objh] = legend([p_L1_1yr, p_L2_1yr, pS]  ,label_L1_1yr, label_L2_1yr,'Surface');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    traj = data_L2_1yr(data_L2_1yr(:,c_trajID)==462, :);
    figure; hold all
    PlotBoi3_CR3Bn(24)
    kk = 842;
    plot3(traj(1:kk,c_x),traj(1:kk,c_y),traj(1:kk,c_z),'color',colors.blue2,'linewidth',1.5);
    plotBody3(secondary.R_n,[1-secondary.MR,0,0],colors.blue2,0)
    axis equal
    zlim([-1 1]*1.2e-2)
    
    
    figure; hold all
    JC = L2FlyoverVelocity_2_JC(50, secondary.MR, L123(2,:), vNorm, 1);
    prms.u = secondary.MR;
    prms.n = 1;
    prms.R2_n = secondary.R_n;
    plotCR3BP_YZNeck( JC, secondary.MR , 1, 0, prms, color_L1_1yr, 2);
    plotCR3BP_YZNeck( JC, secondary.MR , 2, 0, prms, color_L2_1yr, 2);
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    axis equal
    PlotBoi3_CR3Bn(23)
    legend('L1 Neck','L2 Neck')
end

end
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









% L2 and High might be separate here
% figure; hold all
%     % xlim([0 2*(L123(2,1)-(1-secondary.MR))])
%     PlotBoi3('Longitude, $^\circ$', 'Latitude, $^\circ$','$Data3$',23,'LaTex')
%     p_high = plot3(data_high_1pi(:,c_lon), data_high_1pi(:,c_lat), data_high_1pi(:, c_i),'.','color',color_high_1pi);
% %     p_L1   = plot3(nearData_L1(nearRange_L1,c_lon), nearData_L1(nearRange_L1,c_lat), nearData_L1(nearRange_L1, c_i),'.','color',color_L1_1yr);
%     p_L2   = plot3(nearData_L2(:,c_lon), nearData_L2(:,c_lat), nearData_L2(:, c_i),'.','color',color_L2_1yr);
