% ========================================================================
%%% Description
% ========================================================================
% For the CR3BP, symbolically calculate equations of motion, gravitational 
% potential, modified Jacobi constant, dF/dX ("A") matrix, and each method 
% of normalized mean motion according to [1] when various zonal harmonic
% terms (up to J20) for each body are considered
%
% Tip: For a more human-readable print-out of any symbolic equation, use
% Matlab's 'pretty' function. For example, type "pretty(JC)" into the
% command line for a more readable version of the resulting Jacobi constant
% equation
%
% Created: 09/01/19
% Author : Luke Bury, luke.bury@colorado.edu
%
% ========================================================================
%%% References
% ========================================================================
% 1) Luke Bury and Jay McMahon. The Effect of Zonal Harmonics on Dynamical
% Structures in the Circular Restricted Three-Body Problem Near the
% Secondary Body. Celestrial Mechanics and Dynamical Astronomy, 132(45),
% 2020.
%
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
ticWhole = tic;

% ========================================================================
%%% Options/Setup
% ========================================================================
% -------------------------------------------------
%%% ZH options
% -------------------------------------------------
%%% Make vectors of the ZH terms desired for each body. For example, to
%%% include J2, J4, and J6 of the primary, and J2 and J4 of the secondary, 
%%% the vectors would look like this:
% primaryZHRange   = [2, 4, 6];
% secondaryZHRange = [2, 4];

primaryZHRange   = [2, 4];
secondaryZHRange = [2];


% ========================================================================
%%% Symoblically generate equations
% ========================================================================
% -------------------------------------------------
%%% Create symbolic variables
% -------------------------------------------------
%%% General variables. 'real' tells Matlab these are all real numbers
syms x y z xd yd zd mu n r1 r2 R1n R2n real

%%% Zonal harmonic terms for each body. 'p' denotes primary and 's' denotes
%%% secondary. 'real' tells Matlab these are all real numbers
syms J2p J3p J4p J5p J6p J7p J8p J9p J10p J11p J12p J13p J14p J15p J16p J17p J18p J19p J20p real
syms J2s J3s J4s J5s J6s J7s J8s J9s J10s J11s J12s J13s J14s J15s J16s J17s J18s J19s J20s real

% -------------------------------------------------
%%% Build potential functions
% -------------------------------------------------
%%% Initialize gravitational potential functions with Newtonian relations
Vp = (1-mu)/r1;
Vs = mu/r2;

%%% Build Vp (potential due to primary) by looping through primary-related
%%% zonal harmonic terms
for ZHp_i = primaryZHRange
    %%% Grab the current zonal harmonic term  
    if ZHp_i == 2
        J_i = J2p;
    elseif ZHp_i == 3
        J_i = J3p;
    elseif ZHp_i == 4
        J_i = J4p;
    elseif ZHp_i == 5
        J_i = J5p;
    elseif ZHp_i == 6
        J_i = J6p;
    elseif ZHp_i == 7
        J_i = J7p;
    elseif ZHp_i == 8
        J_i = J8p;
    elseif ZHp_i == 9
        J_i = J9p;
    elseif ZHp_i == 10
        J_i = J10p;
    elseif ZHp_i == 11
        J_i = J11p;
    elseif ZHp_i == 12
        J_i = J12p;
    elseif ZHp_i == 13
        J_i = J13p;
    elseif ZHp_i == 14
        J_i = J14p;
    elseif ZHp_i == 15
        J_i = J15p;
    elseif ZHp_i == 16
        J_i = J16p;
    elseif ZHp_i == 17
        J_i = J17p;
    elseif ZHp_i == 18
        J_i = J18p;
    elseif ZHp_i == 19
        J_i = J19p;
    elseif ZHp_i == 20
        J_i = J20p;
    end
    
    %%% Create the potential corresponding to the current zonal harmonic
    %%% term
    newTerm_p = -((1-mu)/r1) * ((R1n^ZHp_i)/(r1^ZHp_i))* J_i * legendreP(ZHp_i, z/r1);
    
    %%% Append new term to primary potential
    Vp = Vp + simplify(newTerm_p);
end


%%% Build Vs (potential due to secondary) by looping through primary-related
%%% zonal harmonic terms
for ZHs_i = secondaryZHRange
    %%% Grab the current zonal harmonic term  
    if ZHs_i == 2
        J_i = J2s;
    elseif ZHs_i == 3
        J_i = J3s;
    elseif ZHs_i == 4
        J_i = J4s;
    elseif ZHs_i == 5
        J_i = J5s;
    elseif ZHs_i == 6
        J_i = J6s;
    elseif ZHs_i == 7
        J_i = J7s;
    elseif ZHs_i == 8
        J_i = J8s;
    elseif ZHs_i == 9
        J_i = J9s;
    elseif ZHs_i == 10
        J_i = J10s;
    elseif ZHs_i == 11
        J_i = J11s;
    elseif ZHs_i == 12
        J_i = J12s;
    elseif ZHs_i == 13
        J_i = J13s;
    elseif ZHs_i == 14
        J_i = J14s;
    elseif ZHs_i == 15
        J_i = J15s;
    elseif ZHs_i == 16
        J_i = J16s;
    elseif ZHs_i == 17
        J_i = J17s;
    elseif ZHs_i == 18
        J_i = J18s;
    elseif ZHs_i == 19
        J_i = J19s;
    elseif ZHs_i == 20
        J_i = J20s;
    end
    
    %%% Create the potential corresponding to the current zonal harmonic
    %%% term
    newTerm_s = -((mu)/r2) * ((R2n^ZHs_i)/(r2^ZHs_i))* J_i * legendreP(ZHs_i, z/r2);
    
    %%% Append new term to the secondary potential
    Vs = Vs + simplify(newTerm_s);
end

%%% Combine potentials from each body
V = Vp + Vs;

%%% Build position vectors and substitute into V
r1 = ((x+mu)^2 + y^2 + z^2)^(1/2);
r2 = ((x-1+mu)^2 + y^2 + z^2)^(1/2);

V = subs(V);

%%% Create system gravitational potential
U_CR3BP = (1/2)*(n^2)*(x^2 + y^2) + V;
% U_2BI = V;

%%% Create Jacobi Constant
JC = -(xd^2 + yd^2 +zd^2) + 2*U_CR3BP;

%%% Differentiate system gravitational potential
dUdx = diff(U_CR3BP,x);
dUdy = diff(U_CR3BP,y);
dUdz = diff(U_CR3BP,z);

dUdx = simplify(dUdx);
dUdy = simplify(dUdy);
dUdz = simplify(dUdz);

%%% Create equations of motion
xdd = 2*n*yd + dUdx;
ydd = -2*n*xd + dUdy;
zdd = dUdz;

EOM_CR3BP = [xdd; ydd; zdd];

%%% Calculate the Jacobian of the dynamics and state
Amat_dFdX = simplify(jacobian([xd;yd;zd; EOM_CR3BP], [x;y;z;xd;yd;zd]));


% ========================================================================
%%% Starting process over so it's all correct for mean motion process
% ========================================================================
% -------------------------------------------------
%%% Redefine symbollic variables and build potentials of each body in
%%% relative terms
% -------------------------------------------------
syms m1 m2 R1n R2n J2p J2s r12 x y z real

%%% Build gravitational potential functions
V1 = m2/r12;
V2 = -m1/r12;

%%% Build V1 (how the primary moves because of the secondary) by looping 
%%% through secondary-related zonal harmonic terms
for ZHs_i = secondaryZHRange
    %%% Grab the current zonal harmonic term  
    if ZHs_i == 2
        J_i = J2s;
    elseif ZHs_i == 3
        J_i = J3s;
    elseif ZHs_i == 4
        J_i = J4s;
    elseif ZHs_i == 5
        J_i = J5s;
    elseif ZHs_i == 6
        J_i = J6s;
    elseif ZHs_i == 7
        J_i = J7s;
    elseif ZHs_i == 8
        J_i = J8s;
    elseif ZHs_i == 9
        J_i = J9s;
    elseif ZHs_i == 10
        J_i = J10s;
    elseif ZHs_i == 11
        J_i = J11s;
    elseif ZHs_i == 12
        J_i = J12s;
    elseif ZHs_i == 13
        J_i = J13s;
    elseif ZHs_i == 14
        J_i = J14s;
    elseif ZHs_i == 15
        J_i = J15s;
    elseif ZHs_i == 16
        J_i = J16s;
    elseif ZHs_i == 17
        J_i = J17s;
    elseif ZHs_i == 18
        J_i = J18s;
    elseif ZHs_i == 19
        J_i = J19s;
    elseif ZHs_i == 20
        J_i = J20s;
    end
    
    %%% Create the potential corresponding to the current zonal harmonic
    %%% term
    newTerm_s = -(m2/r12) * ((R2n^ZHs_i)/(r12^ZHs_i))* J_i * legendreP(ZHs_i, z/r12);
    
    %%% Append new term to V1 potential
    V1 = V1 + simplify(newTerm_s);
end

%%% Build V2 (how the secondary moves because of the primary) by looping 
%%% through primary-related zonal harmonic terms
for ZHp_i = primaryZHRange
    %%% Grab the current zonal harmonic term  
    if ZHp_i == 2
        J_i = J2p;
    elseif ZHp_i == 3
        J_i = J3p;
    elseif ZHp_i == 4
        J_i = J4p;
    elseif ZHp_i == 5
        J_i = J5p;
    elseif ZHp_i == 6
        J_i = J6p;
    elseif ZHp_i == 7
        J_i = J7p;
    elseif ZHp_i == 8
        J_i = J8p;
    elseif ZHp_i == 9
        J_i = J9p;
    elseif ZHp_i == 10
        J_i = J10p;
    elseif ZHp_i == 11
        J_i = J11p;
    elseif ZHp_i == 12
        J_i = J12p;
    elseif ZHp_i == 13
        J_i = J13p;
    elseif ZHp_i == 14
        J_i = J14p;
    elseif ZHp_i == 15
        J_i = J15p;
    elseif ZHp_i == 16
        J_i = J16p;
    elseif ZHp_i == 17
        J_i = J17p;
    elseif ZHp_i == 18
        J_i = J18p;
    elseif ZHp_i == 19
        J_i = J19p;
    elseif ZHp_i == 20
        J_i = J20p;
    end
    
    %%% Create the potential corresponding to the current zonal harmonic
    %%% term
    newTerm_p = (m1/r12) * ((R1n^ZHp_i)/(r12^ZHp_i))* J_i * legendreP(ZHp_i, z/r12);
    
    %%% Append new term to V2 potential
    V2 = V2 + simplify(newTerm_p);
end

%%% Set z=0 since mean motion is a planar concern
z = 0;

%%% Substitute z=0 into the potentials
V1 = subs(V1);
V2 = subs(V2);

%%% Define relative position vector
r12 = sqrt(x^2 + y^2);

%%% Substitude in this position vector so that we can take derivate d/dx
V1 = subs(V1);
V2 = subs(V2);

%%% Take d/dx
r1dd = diff(V1, x);
r2dd = diff(V2, x);

%%% Acquire relative acceleration
r12dd = r2dd - r1dd;

%%% Divide out the extra x that pops out in derivative. In this analytical
%%% problem, we only need to compare the 'x' component of the dynamics and
%%% kinematics, and certain parts cancel out, including a full 'r12', of
%%% which this x is the first component
r12dd_divx = r12dd / x;

%%% Setting y=0 and x=1 allows us to substitute in '1' for r12
y=0;
x=1;

%%% Turn mass values into normalized parameters with mass ratio
m1 = 1-mu;
m2 = mu;

%%% Solve for mean motion [1]
n_squared = subs(r12dd_divx);
n = simplify(sqrt(n_squared));


% ========================================================================
%%% Print results
% ========================================================================
fprintf('------------------------------------------------------------------------------------------------------------\n')
fprintf('Primary ZH terms:   ')
fprintf('%d ',primaryZHRange)
fprintf('\n')
fprintf('Secondary ZH terms: ')
fprintf('%d ',secondaryZHRange)
fprintf('\n')
fprintf('------------------------------------------------------------------------------------------------------------\n')
fprintf('Key:\n')
fprintf('mu               = CR3BP mass ratio\n')
fprintf('[x,y,z,xd,yd,zd] = Classical barycentric CR3BP state\n')
fprintf('n                = Normalized (dimensionless) mean motion\n')
fprintf('Jnb              = Zonal harmonic term where n = zonal harmonic term number, b = ''p'' for primary body or ''s'' for secondary body\n')
fprintf('a                = semimajor axis (km)\n')
fprintf('m1, m2           = masses of primary bodies in kg\n')
fprintf('G                = Gravitational parameter in km^3 * kg^-1 * s^-2\n')
fprintf('R1n, R2n         = Normalized radius of primary and secondary body\n')
fprintf('n_real           = Measured mean motion of system in s^-1\n')
fprintf('------------------------------------------------------------------------------------------------------------\n')
fprintf('CR3BP Equations of Motion:\nxdd = %s\nydd = %s\nzdd = %s\n',EOM_CR3BP(1),EOM_CR3BP(2),EOM_CR3BP(3))
fprintf('------------------------------------------------------------------------------------------------------------\n')
fprintf('CR3BP Gravitational Potential:\nU = %s\n',U_CR3BP)
fprintf('------------------------------------------------------------------------------------------------------------\n')
fprintf('Jacobi Constant:\nJC = %s\n',JC)
fprintf('------------------------------------------------------------------------------------------------------------\n')
fprintf('''A'' matrix (dF/dX) complete (6x6 results stored in Amat_dFdX)\n')
fprintf('------------------------------------------------------------------------------------------------------------\n')
fprintf('Normalized mean motion (two methods):\n\n')
fprintf('Ephemeris-based method [1]:\n')
fprintf('n = n_real*sqrt((a^3)/(G*(m1+m2)))\n\n\n')
fprintf('Theory-based method [1]:\n')
fprintf('n = %s\n', n)
fprintf('------------------------------------------------------------------------------------------------------------\n')

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\n(Elapsed time: %1.4f seconds)\n',tocWhole)


