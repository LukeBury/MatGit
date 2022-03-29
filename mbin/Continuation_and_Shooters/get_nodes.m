function [t_nodes, X_nodes] = get_nodes(X_in, t_in, N, integrator, options, prms)
%%% Description
%       Returns N nodes along a designated trajectory
%       
% ------------------------------------------------------------------------
%%% Inputs
%       X_in       - [6x1] X0 state vector
%       t_in       - [1x2] initial and final time 
%       N          - [1x1] Number of nodes (for POs, needs to be N_des + 1)
%       integrator - [func] Integrator for nodes (ex: '@Int_CR3Bn_ZH')
%       options    - [struct] Abs and Rel accuracy for integrator
%       prms       - [struct] Struct of necessary parameters
% ------------------------------------------------------------------------
%%% Outputs
%       X_out - [Nx6] State vectors at evenly spaced times
%       t_out - [Nx1] evenly space node times 
% ------------------------------------------------------------------------
% Created: 6/4/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Create nodes
% -------------------------------------------------
[t_nodes, X_nodes] = ode113(integrator, linspace(t_in(1), t_in(2), N), X_in, options, prms);

end % function