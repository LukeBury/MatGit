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

% syms u r1 r2 r3 real
% assume(r2>r3)
% assume(r3>r2)
% 
% J_H = sqrt(2*u*r2/((r1+r2)*r1)) - sqrt(u/r1) + sqrt(u/r2) - sqrt(2*u*r1/((r1+r2)*r2));
% 
% J_Bi = sqrt(2*u*r3/((r1+r3)*r1)) - sqrt(u/r1) + sqrt(2*u*r2/((r2+r3)*r3)) - sqrt(2*u*r1/((r1+r3)*r3)) + sqrt(2*u*r3/((r2+r3)*r2)) - sqrt(u/r2);
% 
% normalizer = sqrt(u/r1);
% 
% Jt_H = J_H / normalizer;
% 
% Jt_Bi = J_Bi / normalizer;
% 
% r = r2/r1;
% l = r3/r1;


% --------------------------------
% notes = @(r) sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1) - (sqrt(2)-1)*(1 + 1/sqrt(r));
% 
% 
% me = @(r) sqrt(2/(1+r)) + (1/r)*(2-sqrt(2)) - sqrt(2)/(r*sqrt(1+r)) - sqrt(2/r);
% 
% 
% % Target solution
% sol = @(r) r^3 - (7 + 4*sqrt(2))*r^2 + (3 + 4*sqrt(2))*r - 1;
% 
% 
% % these forms are equal for hohman
% Jt_H_HW    = @(r) sqrt(2*r/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)));
% Jt_H_Notes = @(r) sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1);
% 
% Jt_P = @(r) (sqrt(2)-1)*(1 + 1/sqrt(r))

% --------------------------------
% % % % 
% % % % syms r real
% % % % assume(r>1)
% % % % 
% % % % inequality = (sqrt(2*r/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)))) - ((sqrt(2)-1)*(1 + 1/sqrt(r)));
% % % % inequality_f = @(r) (sqrt(2*r/(1+r)) - 1 + 1/sqrt(r) - sqrt(2/(r*(1+r)))) - ((sqrt(2)-1)*(1 + 1/sqrt(r)));
% % % % 
r = 0:.001:12;
inequalitySolution = r.^3 - (7 + 4*sqrt(2)).*r.^2 + (3 + 4*sqrt(2)).*r - 1;




inequalityFromEquations = (sqrt(2.*r./(1+r)) - 1 + 1./sqrt(r) - sqrt(2./(r.*(1+r)))) - ((sqrt(2)-1).*(1 + 1./sqrt(r)));
% % % % 
% % % % % figure
% % % % 
% % % % 
% % % % 
% % % % 
% % % % figure
% % % % subplot(1,2,1); hold all
% % % % inEq = plot(r,inequalityFromEquations,'b','linewidth',1.5);
% % % % % inSo = plot(x,inequalitySolution,'r','linewidth',1.5);
% % % % PlotBoi2('r','Cost',14,'LaTex')
% % % % 
% % % % subplot(1,2,2); hold all
% % % % inEq = plot(r,inequalityFromEquations,'b','linewidth',1.5);
% % % % inSo = plot(r,inequalitySolution,'r','linewidth',1.5);
% % % % 
% % % % plot([min(r) max(r)],[0 0],'k')
% % % % 
% % % % legend([inSo inEq], 'Given Solution','From Equations')
% % % % PlotBoi2('r','Cost',14,'LaTex')
% % % % ylim([-10 10])
% % % % set(gcf,'Position',[-989 344 560 420]);

%%%% Comparing current equation
y = ((sqrt(2)*((r-1)./(sqrt(r+1)) - sqrt(r) - 1) + 2)).*(r-0.5715).*(r-0.1466).*11.1374.*3;
% y = ((sqrt(2)*((r-1)./(sqrt(r+1)) - sqrt(r) - 1) + 2))

figure; hold all
inEq = plot(r,inequalityFromEquations,'b','linewidth',1.5);
inSo = plot(r,inequalitySolution,'r','linewidth',1.5);
myTest = plot(r,y,'g','linewidth',1.5);

plot([min(r) max(r)],[0 0],'k')

legend([inSo inEq myTest], 'Given Solution','From Equations','My Test')
PlotBoi2('r','Cost',14,'LaTex')
ylim([-5 5])
set(gcf,'Position',[-989 344 560 420]);



% 
% x = 0:.01:10;
% y = 1 - 1./x;
% figure
% plot(x,y)





