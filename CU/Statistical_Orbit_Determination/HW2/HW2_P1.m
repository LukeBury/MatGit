clear
clc
close all

% ------------------------------------------------------------------------
%%% Creating Variabls and Equations
% ------------------------------------------------------------------------
syms R_x R_y R_z Rs_x Rs_y Rs_z V_x V_y V_z Vs_x Vs_y Vs_z

R = [R_x, R_y, R_z];
Rs = [Rs_x, Rs_y, Rs_z];
V = [V_x, V_y, V_z];
Vs = [Vs_x, Vs_y, Vs_z];

X = [R V];

p = ((R_x - Rs_x)^2 + (R_y - Rs_y)^2 + (R_z - Rs_z)^2)^(1/2)

pdot = dot((R - Rs),(V - Vs)) / p


jacobian(p,X)

jacobian(pdot,X)'






