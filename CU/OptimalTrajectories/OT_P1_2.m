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

save_Figures = 1;
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
    
    
% % % syms r l eps real
% % % %%% cost of hohmann
% % % % JH_sym = sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1); % for increasing orbit size
% % % %%% cost of Bi-El
% % % JB_sym = sqrt(2*l/(1+l)) - 1 + sqrt(2*r/(l*(l+r))) - sqrt(2/(l*(1+l))) + 1/sqrt(r) - sqrt(2*l/(r*(l+r)));
% % % 
% % % %%% sub value for l ... looking at l values just slightly less than r
% % % l = r - eps;
% % % 
% % % %%% subbing in l
% % % JB_meps = subs(JB_sym); % (meps ~ r-eps)
% % % 
% % % %%% differentiating wrt epsilon
% % % dfde_meps = diff(JB_meps,eps); % (meps ~ r-eps)
% % % 
% % % %%% evaluating at epsilon = 0 (meps ~ r-eps)
% % % eps = 0;
% % % eval_dfde_meps = simplify(subs(dfde_meps));
% % % 
% % % %%% printing result of evaluating at epsilon = 0
% % % pretty(eval_dfde_meps)
% % % 
% % % %%% Plotting that^ result as a function of r
% % % rr = 1:.01:100;
% % % J_meps = @(r) (1./r).^(3./2) - 2^(1./2)./(2.*(r./(r + 1)).^(1./2).*(r + 1).^2) - (2^(1./2).*(2.*r + 1))./(2.*r.^2.*(1./(r.*(r + 1))).^(1./2).*(r + 1).^2);
% % % figure; hold all
% % % p1 = plot(rr,J_meps(rr),'b','linewidth',2);
% % % 
% % % % J_meps = vpa(subs(eval_dfde_meps,r,rr));
% % % % figure; hold all
% % % % p1 = plot(rr,J_meps,'b','linewidth',2);
% % % 
% % % % -------------------------
% % % 
% % % clear r l eps 
% % % syms r l eps real
% % % %%% cost of Bi-El
% % % JB_sym = sqrt(2*l/(1+l)) - 1 + sqrt(2*r/(l*(l+r))) - sqrt(2/(l*(1+l))) + 1/sqrt(r) - sqrt(2*l/(r*(l+r)));
% % % 
% % % %%% sub value for l ... looking at l values just larger than 1
% % % l = 1 + eps;
% % % 
% % % %%% subbing in l
% % % JB_peps = subs(JB_sym); % (peps ~ 1+eps)
% % % 
% % % %%% differentiating wrt epsilon
% % % dfde_peps = diff(JB_peps,eps); % (peps ~ 1+eps)
% % % 
% % % %%% evaluating at epsilon = 0 (peps ~ 1+eps)
% % % eps = 0;
% % % eval_dfde_peps = simplify(subs(dfde_peps)); % (peps ~ 1+eps)
% % % 
% % % %%% printing result of evaluating at epsilon = 0
% % % pretty(eval_dfde_peps) % (peps ~ 1+eps)
% % % 
% % % %%% Plotting that^ result as a function of r
% % % J_peps1 = @(r) -(1 - (2^(1/2).*r.*(r + 2))./(2.*(r./(r + 1)).^(1./2).*(r + 1).^2) - 2.^(1./2)./(2.*(1./(r.*(r + 1))).^(1./2).*(r + 1).^2));
% % % 
% % % p2 = plot(rr,J_peps1(rr),'r','linewidth',2);
% % % 
% % % legend([p1 p2], 'l = r - \epsilon','l = 1 + \epsilon')
% % % PlotBoi2('$r$','Cost',16,'LaTex')


syms r l eps real
%%% cost of hohmann
% JH_sym = sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1); % for increasing orbit size
%%% cost of Bi-El
JB_sym = sqrt(2*l/(1+l)) - 1 + sqrt(2*r/(l*(l+r))) - sqrt(2/(l*(1+l))) + 1/sqrt(r) - sqrt(2*l/(r*(l+r)));

dfdJ = diff(JB_sym,l);

eps = 0;
l = r-eps;
dfdJ_eval = simplify(subs(dfdJ));
dfdJ_eval
pretty(dfdJ_eval)

J_meps = @(r) 2.^(1/2)./(2.*(r./(r + 1)).^(1/2).*(r + 1).^2) - (1./r).^(3/2) + (2.^(1/2).*(2.*r + 1))./(2.*r.^2.*(1./(r.*(r + 1))).^(1/2).*(r + 1).^2);

rr = 1:.01:100;
figure; hold all
p1 = plot(rr,J_meps(rr),'b','linewidth',2);

% -------------------------
% -------------------------

clear r l eps 
syms r l eps real
JB_sym = sqrt(2*l/(1+l)) - 1 + sqrt(2*r/(l*(l+r))) - sqrt(2/(l*(1+l))) + 1/sqrt(r) - sqrt(2*l/(r*(l+r)));

dfdJ = diff(JB_sym,l);

eps = 0;
l = 1+eps;
dfdJ_eval = simplify(subs(dfdJ));
dfdJ_eval
pretty(-dfdJ_eval)

J_peps = @(r) -(1 - (2.^(1/2).*r.*(r + 2))./(2.*(r./(r + 1)).^(1/2).*(r + 1).^2) - 2.^(1/2)./(2.*(1./(r.*(r + 1))).^(1/2).*(r + 1).^2));

rr = 1:.01:100;
p2 = plot(rr,J_peps(rr),'r','linewidth',2);

legend([p1 p2], 'l = r - \epsilon','l = 1 + \epsilon')
PlotBoi2('$r$','Cost',16,'LaTex')


if save_Figures == 1
    figName = [path_OTfigures,'OT_HW_P1-2_1.png'];
    saveas(gcf,figName)
end


return


end % run_symbolicWork
% ========================================================================
%%% Plotting Section
% ========================================================================
if run_plotWork == 1
% % these forms are equal for hohman
Jt_H = @(r) sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1);
Jt_Bi= @(r,l) abs(sqrt(2*l/(1+l)) - 1) + abs(sqrt(2*r/(l*(l+r))) - sqrt(2/(l*(1+l)))) + abs(sqrt(2*l/(r*(l+r))) - 1/sqrt(r));


r = 1:.001:18;

figure('position',[440 463 774 332])
% --------- Hohmann
f_Jt_H  = @(r) sqrt(2./(r.*(1+r))).*(r-1) - (1./sqrt(r)).*(sqrt(r)-1);
subplot(1,2,1); hold all
hohmanCost = plot(r,f_Jt_H(r),'-','linewidth',2,'color', colors.std.black);

subplot(1,2,2); hold all
hohmanCost = plot(r,f_Jt_H(r),'-','linewidth',2,'color', colors.std.black);
% % --------- l = 1
% l = 1;
% f_Jt_Bi = abs(sqrt(2.*l./(1+l)) - 1) + abs(sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l)))) + abs(sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r));
% 
% p1 = plot(r,f_Jt_Bi,'linewidth',1.5,'color',colors.std.grn);
% plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.grn)

% --------- l = 2
l = 2;
f_Jt_Bi = @(r,l) abs(sqrt(2.*l./(1+l)) - 1) + abs(sqrt(2.*r./(l.*(l+r))) - sqrt(2./(l.*(1+l)))) + abs(sqrt(2.*l./(r.*(l+r))) - 1./sqrt(r));

subplot(1,2,1); hold all
p2 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.blue);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.blue)

subplot(1,2,2); hold all
p2 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.blue);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.blue)
% % --------- l = 3
% l = 3;
% 
% p3 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.red);
% plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.red)

% --------- l = 4
l = 4;

subplot(1,2,1)
p4 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.mag);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.mag)

subplot(1,2,2)
p4 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.mag);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.mag)
% % --------- l = 5
% l = 5;
% 
% p5 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.purp);
% plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.purp)

% --------- l = 6
l = 6;

subplot(1,2,1)
p6 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.brown);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.brown)

subplot(1,2,2)
p6 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.brown);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.brown)
% --------- l = 15.5817
l = 15;

subplot(1,2,1)
p7 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.grn);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.grn)

subplot(1,2,2)
p7 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.grn);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.grn)
% % --------- l = 15.6
% l = 15.6;
% 
% subplot(1,2,1)
% p8 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.red);
% plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.red)
% 
% subplot(1,2,2)
% p8 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.red);
% plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.red)

% --------- l = 17.5
l = 17;

subplot(1,2,1)
p9 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.purp);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.purp)

subplot(1,2,2)
p9 = plot(r,f_Jt_Bi(r,l),'linewidth',1.5,'color',colors.std.purp);
plot([l l],[0 .7],'--','linewidth',1.5,'color',colors.std.purp)

subplot(1,2,1)
PlotBoi2('$r$ Value ($r_2/r_1$)','Cost of Transfer',14,'LaTex')

subplot(1,2,2)
PlotBoi2('$r$ Value ($r_2/r_1$)','',14,'LaTex')
legend([hohmanCost p2 p4 p6 p7 p9],'Hohmann','Bi-E,  l = 2','Bi-E,  l = 4','Bi-E,  l = 6','Bi-E,  l = 15.5817', 'Bi-E,  l = 17.5')
xlim([13 18])
ylim([0.5335 0.5385])
title('Zoomed View')



if save_Figures == 1
    figName = [path_OTfigures,'OT_HW_P1-2_2.png'];
    saveas(gcf,figName)
end


end % run_plotWork










