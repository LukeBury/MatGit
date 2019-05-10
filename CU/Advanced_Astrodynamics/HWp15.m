clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================
run_15a = 0;
run_15b = 0;
run_15c_1 = 0;
run_15c_2 = 0;
run_15f = 1;
%% =======================================================================
%%% 15-a
% ========================================================================
if run_15a == 1
syms a as i w e E real

expression  = (1-e*cos(E))*(sin(w)*(cos(E)-e) + cos(w)*sin(E)*sqrt(1-e^2));

expressionIndefiniteIntegral = int(expression,E);
pretty(simplify(expressionIndefiniteIntegral))

expressionDefiniteIntegral = int(expression,E,[0 2*pi])

fullFinal = expressionDefiniteIntegral * (a*as*sin(i))/(2*pi);
pretty(simplify(fullFinal))


end % run_15a == 1


%% =======================================================================
%%% 15-b
% ========================================================================
if run_15b == 1
clear
syms a as i w e E f n real

%%% aDot
fprintf('aDot:\n')
R = a*(1-e^2)*sin(i)*sin(w+f)*as / (1 + e*cos(f));

dRdf = diff(R,f);
dRdf = simplify(dRdf);
% pretty(dRdf)

dfdE = (1 + e*cos(f))/(sqrt(1-e^2));
dEds = 1 / (1 - e*((e + cos(f)) / (1 + e*cos(f))));
% cos(E) = (e + cos(f)) / (1 + e*cos(f))

dRds = dRdf * dfdE * dEds;

aDot = (2/(n*a)) * dRds;
aDot = simplify(aDot);
pretty(aDot)
%%% eDot
fprintf('eDot:\n')
dRdw = diff(R,w);
dRdw = simplify(dRdw);

eDot = (1 / (n*a*a*e)) * ((1-e^2)*dRds - sqrt(1-e^2)*dRdw);
eDot = simplify(eDot);
pretty(eDot)

end % run_15a == 1

%% =======================================================================
%%% 15-c1
% ========================================================================
if run_15c_1 == 1
clear
syms a as i w e E f n real

eqn = (cos(w)*(cos(E)-e)/(1-e*cos(E)) - sin(w)*sqrt(1-e^2)*sin(E)/(1-e*cos(E)) + e*cos(w))*(1 - e*cos(E));

eqn_int = int(eqn,E);

delta_a = eqn_int * 2*sin(i)*as/(n*n*sqrt(1-e^2));
delta_a = simplify(delta_a);
pretty(delta_a)


n = 20; e = 0.3; i = pi/4; w = 3.4; f = 2.2; as = 1.4; E = 3.15;
delta_a_me = @(n,e,i,w,f,as,E) (4*as*cos(E/2)*sin(i)*(sin(E/2)*cos(w) - e^2*sin(E/2)*cos(w) + cos(E/2)*sin(w)*(1 - e^2)^(1/2)))/(n^2*(1 - e^2)^(1/2));
delta_a_real = @(n,e,i,w,f,as,E) (2*as*sin(i)/(n*n))*(sqrt(1-e^2)*cos(w)*sin(E) - sin(w)*(1-cos(E)));


delta_a_me2 = @(n,e,i,w,f,as,E) (2*as*sin(i)/(n*n*sqrt(1-e^2))) * (cos(w)*sin(E) + sin(w)*sqrt(1-e^2)*(cos(E)-1) - e*e*cos(w)*sin(E));
delta_a_me(n,e,i,w,f,as,E)
delta_a_real(n,e,i,w,f,as,E)
delta_a_me2(n,e,i,w,f,as,E)
end % run_15c_1 == 1

%% =======================================================================
%%% 15-c2
% ========================================================================
if run_15c_2 == 1
clear
syms a as i w e E f n real

eqn = cos(w) - e*cos(w)*cos(E) + cos(w)*(cos(E)-e)*cos(E) - sin(w)*sqrt(1-e^2)*sin(E)*cos(E);

eqn_int = int(eqn,E);

pretty(simplify(eqn_int))

delta_e = eqn_int * as*sqrt(1-e^2)*sin(i)/(n*n*a);

n = 20; e = 0.3; i = pi/4; w = 3.4; a = 2.2; as = 1.4; E = 3.15;
delta_e_me = @(n,e,i,w,a,as,E) (as*sin(i)*(1 - e^2)^(1/2)*((3*E*cos(w))/2 - cos(E/2)*sin(E/2)*cos(w) - 2*cos(E/2)^2*sin(w)*(1 - e^2)^(1/2) + 2*cos(E/2)^4*sin(w)*(1 - e^2)^(1/2) + 2*cos(E/2)^3*sin(E/2)*cos(w) - 4*e*cos(E/2)*sin(E/2)*cos(w)))/(a*n^2);
delta_e_real = @(n,e,i,w,a,as,E) (as*sqrt(1-e^2)*sin(i)/(n*n*a))*((3/2)*cos(w)*E - 2*e*cos(w)*sin(E) + (1/4)*cos(w)*sin(2*E) - (sqrt(1-e^2)/4)*sin(w)*(1-cos(2*E)));

delta_e_me(n,e,i,w,a,as,E)
delta_e_real(n,e,i,w,a,as,E)
% test_eqn = @(E) (3*E*cos(w))/2 - cos(E/2)*sin(E/2)*cos(w) - 2*cos(E/2)^2*sin(w)*(1 - e^2)^(1/2) + 2*cos(E/2)^4*sin(w)*(1 - e^2)^(1/2) + 2*cos(E/2)^3*sin(E/2)*cos(w) - 4*e*cos(E/2)*sin(E/2)*cos(w);
% 
% test2_eqn = test_eqn(E) - test_eqn(0);
% pretty(simplify(test2_eqn))


end % run_15c_2 == 1



%% =======================================================================
%%% 15-f
% ========================================================================
if run_15f == 1
clear
syms a as i w e E f n real

eqn = (sqrt(1-e^2)*cos(w)*sin(E) - sin(w) + sin(w)*cos(E))*(1-e*cos(E));

eqn_int = int(eqn,E,[0, 2*pi]);

delta_a_bar = simplify(eqn_int * (as*sin(i)/(pi*n*n)));
pretty(delta_a_bar)

end % run_15f == 1





