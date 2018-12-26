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
%%% Setup
% ========================================================================

% -----------------------------
% Setup
% -----------------------------
%%% Setting 3B systems to loop through
systems = {bodies.earth,  bodies.moon,       colors.std.blue;...
          bodies.mars,    bodies.phobos,     colors.std.blue;...
          bodies.mars,    bodies.deimos,     colors.std.blue;...
          bodies.jupiter, bodies.io,         colors.std.blue;... 
          bodies.jupiter, bodies.europa,     colors.std.blue;...
          bodies.jupiter, bodies.ganymede,   colors.std.blue;...
          bodies.jupiter, bodies.callisto,   colors.std.blue;...
          bodies.saturn,  bodies.enceladus,  colors.std.blue;...
          bodies.saturn,  bodies.mimas,      colors.std.blue;...
          bodies.saturn,  bodies.rhea,       colors.std.blue;...
          bodies.saturn,  bodies.titan,      colors.std.blue;...
          bodies.uranus,  bodies.titania,    colors.std.blue;...
          bodies.uranus,  bodies.oberon,     colors.std.blue;...
          bodies.neptune, bodies.triton,     colors.std.blue};


%%% Setting column accessers
c_p    = 1;
c_s    = 2;
c_col  = 3;
c_Rn   = 4;
c_SR   = 5;
c_Rn_n = 6;
c_SR_n = 7;
c_MR   = 8;
% -----------------------------
% Running Analysis
% -----------------------------
%%% Looping through systems and building table of L123 differences
for kk = 1:size(systems,1)
    %%% Reassigning data for clarity
    primary   = systems{kk,c_p};
    secondary = systems{kk,c_s};
    
    %%% Assigning rNorm
    rNorm = secondary.a; % n <-> km
    
    %%% Acquire Lagrange points of both models
    L123 = EquilibriumPoints(secondary.MR,1:3);

    %%% Finding ratio between R_n and interior space
    x_int = L123(2,1) - L123(1,1);
    sizeRatio = secondary.R_n/x_int;
    
    %%% Putting info into cell array for ease
    systems{kk,c_Rn} = secondary.R_n;
    systems{kk,c_SR} = sizeRatio;
    systems{kk,c_MR} = secondary.MR;
    
    figure(1)
    subplot(2,1,1); hold all
    plot(secondary.MR, sizeRatio, '.', 'markersize', 18, 'color', systems{kk,c_col});
    PlotBoi2('Mass Ratio','Size Ratio',16,'LaTex')
    set(gca,'xscale','log')
    
    subplot(2,1,2); hold all
    plot(secondary.MR, secondary.R_n, '.', 'markersize', 18, 'color', systems{kk,c_col});
    PlotBoi2('Mass Ratio','R$_n$',16,'LaTex')
    set(gca,'xscale','log')
    
    figure(2); hold all
    plot(secondary.R_n, sizeRatio, '.', 'markersize', 18, 'color', systems{kk,c_col});
    PlotBoi2('R$_n$','Size Ratio',16,'LaTex')
    
end

Rns = [systems{:,c_Rn}];
SRs = [systems{:,c_SR}];

Rn_normalizer = max(Rns);
SR_normalizer = max(SRs);

names = {};
%%% Looping through systems and looking at normalized Rn and SR magnitudes
for kk = 1:size(systems,1)
    systems{kk,c_Rn_n} = systems{kk,c_Rn}/Rn_normalizer;
    systems{kk,c_SR_n} = systems{kk,c_SR}/SR_normalizer;
    
    names{kk} = systems{kk,c_s}.title;
    
    figure(3); hold all
    
    p1 = plot(kk, systems{kk,c_Rn_n}, 'o', 'markersize',10, 'markeredgecolor', colors.std.black, 'markerfacecolor', colors.std.blue);
    p2 = plot(kk, systems{kk,c_SR_n}, 'o', 'markersize',10, 'markeredgecolor', colors.std.black, 'markerfacecolor', colors.std.red);
    
end


Rn_ns = [systems{:,c_Rn_n}];
SR_ns = [systems{:,c_SR_n}];
MRs   = [systems{:,c_MR}];

x = linspace(1,size(systems,1),size(systems,1));

legend([p1 p2], 'R_n','SR')
set(gca,'xticklabel',names)
PlotBoi2('Moon','Normalized R$_n$ and Size Ratio',16,'LaTex')

R = corrcoef(Rns,SRs)
fprintf('The correlation coefficient between these two parameters\n')
fprintf('is apprximately -0.4, which indicates slight negative\n')
fprintf('correlation\n')

mat = [Rn_ns;SR_ns; x];
mat_sort_Rn = sortrows(mat',1)';
mat_sort_MR = sortrows([MRs; x]',1)';



% -----------------------------
%%% Plotting sorted by size ratio
% -----------------------------
figure; hold all
p1 = plot(x,mat_sort_Rn(1,:), 'o', 'markersize',10, 'markeredgecolor', colors.std.black, 'markerfacecolor', colors.std.blue);
p2 = plot(x,mat_sort_Rn(2,:), 'o', 'markersize',10, 'markeredgecolor', colors.std.black, 'markerfacecolor', colors.std.red);
set(gca,'xticklabel',names(mat_sort_Rn(3,:)),'xtick',x)
xtickangle(270);
PlotBoi2('Moon','Normalized R$_n$ and Size Ratio',16,'LaTex')

fit1 = polyfit(x,mat_sort_Rn(1,:),1);

fit2 = polyfit(x,mat_sort_Rn(2,:),1);

plot(x,fit1(2) + fit1(1)*x,'color',colors.std.blue);
plot(x,fit2(2) + fit2(1)*x,'color',colors.std.red);

legend([p1 p2], 'R_n','SR')


% -----------------------------
%%% Plotting sorted by mass ratio
% -----------------------------
figure; hold all
p1 = plot(x,mat_sort_MR(1,:), 'o', 'markersize',10, 'markeredgecolor', colors.std.black, 'markerfacecolor', colors.std.blue);
set(gca,'xticklabel',names(mat_sort_MR(2,:)),'xtick',x)
xtickangle(270);
PlotBoi2('Moon','Mass Ratio',16,'LaTex')
setLogPlot()
legend([p1], 'Mass Ratio')


% %%% Setting 3B systems to loop through
% % systems = {bodies.earth,  bodies.moon,       colors.std.orange;...
% %           bodies.mars,    bodies.phobos,     colors.std.grn;...
% %           bodies.mars,    bodies.deimos,     colors.std.blue;...
% %           bodies.jupiter, bodies.io,         colors.std.brown;... 
% %           bodies.jupiter, bodies.europa,     colors.std.red;...
% %           bodies.jupiter, bodies.ganymede,   colors.std.black;...
% %           bodies.saturn,  bodies.enceladus,  colors.std.mag;...
% %           bodies.saturn,  bodies.titan,      colors.std.pink;...
% %           bodies.neptune, bodies.triton,     colors.std.purp};
% systems = {bodies.mars,    bodies.phobos,     colors.std.grn;...
%           bodies.mars,    bodies.deimos,     colors.std.blue;...
%           bodies.jupiter, bodies.io,         colors.std.brown;... 
%           bodies.jupiter, bodies.europa,     colors.std.red;...
%           bodies.jupiter, bodies.ganymede,   colors.std.black;...
%           bodies.saturn,  bodies.enceladus,  colors.std.mag;...
%           bodies.saturn,  bodies.titan,      colors.std.pink;...
%           bodies.neptune, bodies.triton,     colors.std.purp};
% 
% %%% Setting column accessers
% c_p = 1;
% c_s = 2;
% c_col = 3;
% c_Rn  = 4;
% c_SR  = 5;
% c_Rn_n = 6;
% c_SR_n = 7;
% % -----------------------------
% % Running Analysis
% % -----------------------------
% %%% Looping through systems and building table of L123 differences
% for kk = 1:size(systems,1)
%     %%% Reassigning data for clarity
%     primary   = systems{kk,c_p};
%     secondary = systems{kk,c_s};
%     
%     %%% Assigning rNorm
%     rNorm = secondary.a; % n <-> km
%     
%     %%% Acquire Lagrange points of both models
%     L123 = EquilibriumPoints(secondary.MR,1:3);
% 
%     %%% Finding ratio between R_n and interior space
%     x_int = L123(2,1) - L123(1,1);
%     sizeRatio = secondary.R_n/x_int;
%     
%     %%% Putting info into cell array for ease
%     systems{kk,c_Rn} = secondary.R_n;
%     systems{kk,c_SR} = sizeRatio;
%     
%     figure(1)
%     subplot(2,1,1); hold all
%     plot(secondary.MR, sizeRatio, '.', 'markersize', 18, 'color', systems{kk,c_col});
%     PlotBoi2('Mass Ratio','Size Ratio',16,'LaTex')
%     set(gca,'xscale','log')
%     
%     subplot(2,1,2); hold all
%     plot(secondary.MR, secondary.R_n, '.', 'markersize', 18, 'color', systems{kk,c_col});
%     PlotBoi2('Mass Ratio','R$_n$',16,'LaTex')
%     set(gca,'xscale','log')
%     
%     figure(2); hold all
%     plot(secondary.R_n, sizeRatio, '.', 'markersize', 18, 'color', systems{kk,c_col});
%     PlotBoi2('R$_n$','Size Ratio',16,'LaTex')
%     
% end
% 
% Rns = [systems{:,c_Rn}];
% SRs = [systems{:,c_SR}];
% 
% Rn_normalizer = max(Rns);
% SR_normalizer = max(SRs);
% 
% names = {};
% %%% Looping through systems and looking at normalized Rn and SR magnitudes
% for kk = 1:size(systems,1)
%     systems{kk,c_Rn_n} = systems{kk,c_Rn}/Rn_normalizer;
%     systems{kk,c_SR_n} = systems{kk,c_SR}/SR_normalizer;
%     
%     names{kk} = systems{kk,c_s}.name;
%     
%     figure(3); hold all
%     
%     p1 = plot(kk, systems{kk,c_Rn_n}, 'o', 'markersize',10, 'markeredgecolor', colors.std.black, 'markerfacecolor', colors.std.blue);
%     p2 = plot(kk, systems{kk,c_SR_n}, 'o', 'markersize',10, 'markeredgecolor', colors.std.black, 'markerfacecolor', colors.std.red);
%     
% end
% 
% legend([p1 p2], 'R_n','SR')
% set(gca,'xticklabel',names)
% PlotBoi2('Moon','Normalized R$_n$ and Size Ratio',16,'LaTex')
% 
% R = corrcoef(Rns,SRs)
% fprintf('The correlation coefficient between these two parameters is -0.4971\n')
% fprintf('Which indicates slight negative correlation')
% 
% mat = [[systems{:,c_Rn_n}];[systems{:,c_SR_n}]; linspace(1,size(systems,1),size(systems,1))];
% mat = sortrows(mat',1)';
% 
% x = linspace(1,size(systems,1),size(systems,1));
% 
% 
% figure; hold all
% p1 = plot(x,mat(1,:), 'o', 'markersize',10, 'markeredgecolor', colors.std.black, 'markerfacecolor', colors.std.blue);
% p2 = plot(x,mat(2,:), 'o', 'markersize',10, 'markeredgecolor', colors.std.black, 'markerfacecolor', colors.std.red);
% set(gca,'xticklabel',names(mat(3,:)))
% PlotBoi2('Moon','Normalized R$_n$ and Size Ratio',16,'LaTex')
% 
% fit1 = polyfit(x,mat(1,:),1)
% 
% fit2 = polyfit(x,mat(2,:),1)
% 
% plot(x,fit1(2) + fit1(1)*x,'color',colors.std.blue);
% plot(x,fit2(2) + fit2(1)*x,'color',colors.std.red);
% 
% legend([p1 p2], 'R_n','SR')














