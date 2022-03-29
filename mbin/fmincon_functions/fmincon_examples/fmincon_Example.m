close all
clear all
clc
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))


global MU AU Sun_mu
MU = 1.32712440018e11;
AU = 1.49597870700e8;
Sun_mu = 1.32712440018e11;

%This program will run fmincon to optimize a trajectory from Earth to Mars
%based on several inputs. 

%X is the state vector. X(1) = Earth Departure Date. X(2) = Mars Arrival
%Date. 
X = [2457400, 2457600]; 

%These are the lower bounds (lb) and upper bounds (ub) on the searching
%region. Fmincon will only find solutions within the dates bounded by lb
%and ub.
lb = [X(1) - 100; X(2) - 100];
ub = lb + [500; 800];

%%%%%%%%%%%%%% Options for fmincon %%%%%%%%%%%%%%%%%%%%%%
%See Help fminopts for setting up the options
fminopts = optimoptions(@fmincon,'TolCon', 1e-12, 'TolX', 1e-8, 'Algorith', 'sqp', 'Display', 'iter');

%These are linear inequalities. Not necessary for this simple problem. 
A = []; b = []; Aeq = []; beq = [];

fnom = fmincon_example_function(X);  %The nominal value of the function evalulation at the initial guess. 
f2   = @(X)fmincon_example_function(X);

%Running fmincon
[Xoutc, fvalc] = fmincon(f2, X, A, b, Aeq, beq, lb, ub, [], fminopts);


%%%%%%%%%%%%%%%% Compute the performance of the trajectory %%%%%%%%%%%
[C3Launch, Vinf_in, TOFd] = fmincon_example_function_output(Xoutc);
fprintf('Trajectory Performance\n')
fprintf('**************************************\n')
fprintf('C3 = %3.3f km^2/s^2\n', C3Launch)
fprintf('V-infinity arrival at Mars = %2.3f km/s\n', Vinf_in);
fprintf('TOF = %2.3f days\n', TOFd)
