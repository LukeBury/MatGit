clear
clc
close all
addpath(genpath('/Users/CU_Google_Drive/lukebury/Documents/MatGit/mbin'))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Running
% ========================================================================


% % these forms are equal for hohman
Jt_H = @(r) sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1);
Jt_Bi= @(r,l) sqrt(2*l/(1+l)) - 1 + sqrt(2*r/(l*(l+r))) - sqrt(2/(l*(1+l))) + sqrt(2*l/(r*(l+r))) - 1/sqrt(r);


r = 1:.01:10;

figure; hold all
% --------- l = 2
l = 2;
f_Jt_H  = sqrt(2./(r.*(1+r))).*(r-1) - (1./sqrt(r)).*(sqrt(r)-1);
f_Jt_Bi = sqrt(2.*l./(1+l)) - 1 + sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l))) + sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r);

plot(r,f_Jt_H,'-','linewidth',1.5,'color', colors.std.black);
p2 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.blue);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.blue)

% --------- l = 3
l = 3;
f_Jt_Bi = sqrt(2.*l./(1+l)) - 1 + sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l))) + sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r);

p3 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.red);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.red)

% --------- l = 4
l = 4;
f_Jt_Bi = sqrt(2.*l./(1+l)) - 1 + sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l))) + sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r);

p4 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.mag);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.mag)

% --------- l = 5
l = 5;
f_Jt_Bi = sqrt(2.*l./(1+l)) - 1 + sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l))) + sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r);

p5 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.purp);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.purp)

% --------- l = 6
l = 6;
f_Jt_Bi = sqrt(2.*l./(1+l)) - 1 + sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l))) + sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r);

p6 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.brown);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.brown)


PlotBoi2('r Value ($r_2/r_1$)','Cost of Transfer',14,'LaTex')
legend([p2 p3 p4 p5 p6],'l = 2','l = 3','l = 4','l = 5','l = 6')


fprintf('It appears that the cost of the Bi-Elliptic transfer is always\n less when 1 < l < r in this range and when l=r, the costs are equal\n')

