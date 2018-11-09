clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================
run_symbolicWork = 1;
run_plotWork     = 1;

save_Figure = 1;
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Paths
% ========================================================================
path_OTfigures = '/Users/lukebury/CU_Google_Drive/Documents/School/CU/Courses/5-6020_Optimal_Trajectories/HW/HW_Submission/Figures_OT/';

% ========================================================================
%%% Symbolic work
% ========================================================================
if run_symbolicWork == 1
    
    
syms r l eps real
% JH_sym = sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1); % for increasing orbit size
JB_sym = sqrt(2*l/(1+l)) - 1 + sqrt(2*r/(l*(l+r))) - sqrt(2/(l*(1+l))) + 1/sqrt(r) - sqrt(2*l/(r*(l+r)));

l = r - eps;
JB_meps = subs(JB_sym);
dfde_meps = diff(JB_meps,eps);
eps = 0;
eval_dfde_meps = subs(dfde_meps);

pretty(simplify(eval_dfde_meps))
989
return
rr = 1:.01:50;
J_meps = (1./rr).^(3/2) - 2^(1/2)./(2.*(rr./(rr + 1)).^(1/2).*(rr + 1).^2) - (2^(1/2).*(2.*rr + 1))./(2.*rr.^2.*(1./(rr.*(rr + 1))).^(1/2).*(rr + 1).^2);
figure(1); hold all
p1 = plot(rr,J_meps,'b');

% -------------------------

clear r l eps 
syms r l eps real
JB_sym = sqrt(2*l/(1+l)) - 1 + sqrt(2*r/(l*(l+r))) - sqrt(2/(l*(1+l))) + 1/sqrt(r) - sqrt(2*l/(r*(l+r)));
l = 1 + eps;
JB_peps = subs(JB_sym);
dfde_peps = diff(JB_peps,eps);
eps = 0;
eval_dfde_peps = subs(-dfde_peps);

pretty(simplify(eval_dfde_peps))

rr = 1:.01:50;
% J_peps = 1 - (2^(1/2).*rr.*(rr + 2))./(2.*(rr./(rr + 1)).^(1/2).*(rr + 1).^2) - 2^(1/2)./(2*(1./(rr.*(rr + 1))).^(1/2).*(rr + 1).^2);
J_peps = 2^(1/2)./(2*(1./(rr.*(rr + 1))).^(1/2).*(rr + 1).^2) + (2^(1/2).*rr.*(rr + 2))./(2.*(rr./(rr + 1)).^(1/2).*(rr + 1).^2) - 1;
% figure(2)
p2 = plot(rr,J_peps,'r');

legend([p1 p2], 'l = r - \epsilon','l = 1 + \epsilon')
PlotBoi2('r','Cost',16,'LaTex')








end % run_symbolicWork
% ========================================================================
%%% Plotting Section
% ========================================================================
if run_plotWork == 1
% % these forms are equal for hohman
Jt_H = @(r) sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1);
Jt_Bi= @(r,l) abs(sqrt(2*l/(1+l)) - 1) + abs(sqrt(2*r/(l*(l+r))) - sqrt(2/(l*(1+l)))) + abs(sqrt(2*l/(r*(l+r))) - 1/sqrt(r));


r = 1:.001:18;

figure; hold all
% --------- Hohmann
f_Jt_H  = sqrt(2./(r.*(1+r))).*(r-1) - (1./sqrt(r)).*(sqrt(r)-1);
hohmanCost = plot(r,f_Jt_H,'-','linewidth',1.5,'color', colors.std.black);

% --------- l = 1
l = 1;
f_Jt_Bi = abs(sqrt(2.*l./(1+l)) - 1) + abs(sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l)))) + abs(sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r));

p1 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.grn);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.grn)

% --------- l = 2
l = 2;
f_Jt_Bi = abs(sqrt(2.*l./(1+l)) - 1) + abs(sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l)))) + abs(sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r));

p2 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.blue);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.blue)

% --------- l = 3
l = 3;
f_Jt_Bi = abs(sqrt(2.*l./(1+l)) - 1) + abs(sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l)))) + abs(sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r));

p3 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.red);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.red)

% --------- l = 4
l = 4;
f_Jt_Bi = abs(sqrt(2.*l./(1+l)) - 1) + abs(sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l)))) + abs(sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r));

p4 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.mag);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.mag)

% --------- l = 5
l = 5;
f_Jt_Bi = abs(sqrt(2.*l./(1+l)) - 1) + abs(sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l)))) + abs(sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r));

p5 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.purp);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.purp)

% --------- l = 6
l = 6;
f_Jt_Bi = abs(sqrt(2.*l./(1+l)) - 1) + abs(sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l)))) + abs(sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r));

p6 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.brown);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.brown)

% --------- l = 16
l = 16;
f_Jt_Bi = abs(sqrt(2.*l./(1+l)) - 1) + abs(sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l)))) + abs(sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r));

p16 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.maglt);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.maglt)


PlotBoi2('r Value ($r_2/r_1$)','Cost of Transfer',14,'LaTex')
legend([hohmanCost p1 p2 p3 p4 p5 p6 p16],'Hohmann','Bi-E,  l = 1','Bi-E,  l = 2','Bi-E,  l = 3','Bi-E,  l = 4','Bi-E,  l = 5','Bi-E,  l = 6','Bi-E,  l = 16')


title('Left of the hohmann, l>r     Right of the hohmann, r > l')


if save_Figure == 1
    figName = [path_OTfigures,'OT_HW_P2.png'];
    saveas(gcf,figName)
end


end % run_plotWork










