% ========================================================================
%%% Description
% ========================================================================
% This script will generate and print the CR3BP equations of motion,
% gravitational potential, and modified Jacobi constant for a system with
% zonal harmonics of the primary and secondary body considered

% Created: 09/01/10
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
ticWhole = tic;

% ========================================================================
%%% Input
% ========================================================================
% -------------------------------------------------
%%% Zonal Harmonic options
% -------------------------------------------------
%%% Choose zonal harmonic terms for each body (as a vector)
primaryZHRange   = [2, 4, 6];
secondaryZHRange = [2];

% ========================================================================
%%% Generate Equations
% ========================================================================
% -------------------------------------------------
%%% Create symbolic variables
% -------------------------------------------------
syms x y z xd yd zd u R1 R2 r1 r2 real
syms J2p J3p J4p J5p J6p J7p J8p J9p J10p J11p J12p J13p J14p J15p J16p J17p J18p J19p J20p real
syms J2s J3s J4s J5s J6s J7s J8s J9s J10s J11s J12s J13s J14s J15s J16s J17s J18s J19s J20s real

% -------------------------------------------------
%%% Build potential functions
% -------------------------------------------------
%%% Build gravitational potential functions
Vp = (1-u)/r1;
Vs = u/r2;

%%% Build Vp
for ZHp_i = primaryZHRange
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

    newTerm_p = -((1-u)/r1) * ((R1^ZHp_i)/(r1^ZHp_i))* J_i * legendreP(ZHp_i, z/r1);
    
    Vp = Vp + simplify(newTerm_p);
end


%%% Build Vs
for ZHs_i = secondaryZHRange
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

    newTerm_s = -((u)/r2) * ((R2^ZHs_i)/(r2^ZHs_i))* J_i * legendreP(ZHs_i, z/r2);
    
    Vs = Vs + simplify(newTerm_s);
end

%%% Combine
V = Vp + Vs;

%%% Build position vectors 
r1 = ((x+u)^2 + y^2 + z^2)^(1/2);
r2 = ((x-1+u)^2 + y^2 + z^2)^(1/2);

V = subs(V);

%%% Create system gravitational potential
U_CR3BP = (1/2)*(x^2 + y^2) + V;
U_2BI = V;


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
xdd = 2*yd + dUdx;
ydd = -2*xd + dUdy;
zdd = dUdz;

EOM_CR3BP = [xdd; ydd; zdd];

fprintf('CR3BP Equations of Motion:\n%s\n%s\n%s\n\n',EOM_CR3BP(1),EOM_CR3BP(2),EOM_CR3BP(3))
fprintf('CR3BP Gravitational Potential:\n%s\n\n',U_CR3BP)
fprintf('Jacobi Constant:\n%s\n\n',JC)

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\n(Elapsed time: %1.4f seconds)\n',tocWhole)

