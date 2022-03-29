function JCs_vec = getJacobiConstant_ZH(states, prms)
%%% Description
%       For determining the Jacobi constants of states in the CR3BP, with
%       the option of accounting for J2-J6 of either body. If the input
%       structure "prms" includes any of these terms, they will be
%       accounted for in the modified Jacobi constant.
%       
% ------------------------------------------------------------------------
%%% Inputs
%       states - [nx6] Matrix of states to evaluate Jacobi Constants of
%       prms - [struct] fields: u, n, R1, R2, J2p, J2s, J3p, J3s, J4p, J4s
% ------------------------------------------------------------------------
%%% Outputs
%       JCs_vec - [nx1] Vector of jacobi constants corresponding to the
%                       input states
% ------------------------------------------------------------------------
% Created: 09/25/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================

if isequal(size(states), [6,1])
    states = states';
end
% -------------------------------------------------
%%% Preallocate output
% -------------------------------------------------
n_states = size(states, 1);
JCs_vec = NaN(n_states, 1);

% -------------------------------------------------
%%% Scan fields of prms and compute necessary components
% -------------------------------------------------
%%% Shortcut/clarity variables
u = prms.u;

if isfield(prms, 'R1')
    R1 = prms.R1;
end
if isfield(prms, 'R2')
    R2 = prms.R2;
end

%%% Loop through states
for kk = 1:n_states
    % -----------------------------
    % Set shortcut/clarity variables
    % -----------------------------
    x  = states(kk,1);
    y  = states(kk,2);
    z  = states(kk,3);
    xd = states(kk,4);
    yd = states(kk,5);
    zd = states(kk,6);
    
    % -----------------------------
    % Compute zonal harmonic terms
    % -----------------------------
    if isfield(prms, 'J2p')
        JC_J2p = -(prms.J2p*R1^2*(u - 1)*((u + x)^2 + y^2 - 2*z^2))/((u + x)^2 + y^2 + z^2)^(5/2);
    else
        JC_J2p = 0;
    end

    if isfield(prms, 'J2s')
        JC_J2s = (prms.J2s*R2^2*u*((u + x - 1)^2 + y^2 - 2*z^2))/((u + x - 1)^2 + y^2 + z^2)^(5/2);
    else
        JC_J2s = 0;
    end
    
    if isfield(prms, 'J3p')
        JC_J3p = -(prms.J3p*R1^3*z*(u - 1)*(3*(u + x)^2 + 3*y^2 - 2*z^2))/((u + x)^2 + y^2 + z^2)^(7/2);
    else
        JC_J3p = 0;
    end

    if isfield(prms, 'J3s')
        JC_J3s = (prms.J3s*R2^3*u*z*(3*(u + x - 1)^2 + 3*y^2 - 2*z^2))/((u + x - 1)^2 + y^2 + z^2)^(7/2);
    else
        JC_J3s = 0;
    end

    if isfield(prms, 'J4p')
        JC_J4p = (prms.J4p*R1^4*(u - 1)*(3*((u + x)^2 + y^2 + z^2)^2 + 35*z^4 - 30*z^2*((u + x)^2 + y^2 + z^2)))/(4*((u + x)^2 + y^2 + z^2)^(9/2));
    else
        JC_J4p = 0;
    end

    if isfield(prms, 'J4s')
        JC_J4s = -(prms.J4s*R2^4*u*(3*((u + x - 1)^2 + y^2 + z^2)^2 + 35*z^4 - 30*z^2*((u + x - 1)^2 + y^2 + z^2)))/(4*((u + x - 1)^2 + y^2 + z^2)^(9/2));
    else
        JC_J4s = 0;
    end

    if isfield(prms, 'J5p')
        JC_J5p = (prms.J5p*R1^5*z*(u - 1)*(15*((u + x)^2 + y^2 + z^2)^2 + 63*z^4 - 70*z^2*((u + x)^2 + y^2 + z^2)))/(4*((u + x)^2 + y^2 + z^2)^(11/2));
    else
        JC_J5p = 0;
    end

    if isfield(prms, 'J5s')
        JC_J5s = -(prms.J5s*R2^5*u*z*(15*((u + x - 1)^2 + y^2 + z^2)^2 + 63*z^4 - 70*z^2*((u + x - 1)^2 + y^2 + z^2)))/(4*((u + x - 1)^2 + y^2 + z^2)^(11/2));
    else
        JC_J5s = 0;
    end
    
    if isfield(prms, 'J6p')
        JC_J6p = -(prms.J6p*R1^6*(u - 1)*(5*((u + x)^2 + y^2 + z^2)^3 - 105*z^2*((u + x)^2 + y^2 + z^2)^2 - 231*z^6 + 315*z^4*((u + x)^2 + y^2 + z^2)))/(8*((u + x)^2 + y^2 + z^2)^(13/2));
    else
        JC_J6p = 0;
    end

    if isfield(prms, 'J6s')
        JC_J6s = (prms.J6s*R2^6*u*(5*((u + x - 1)^2 + y^2 + z^2)^3 - 105*z^2*((u + x - 1)^2 + y^2 + z^2)^2 - 231*z^6 + 315*z^4*((u + x - 1)^2 + y^2 + z^2)))/(8*((u + x - 1)^2 + y^2 + z^2)^(13/2));
    else
        JC_J6s = 0;
    end
    
    %%% Add terms together
    JC_vanilla = (prms.n^2)*(x^2 + y^2) - xd^2 - yd^2 - zd^2 - (2*(u - 1))/((u + x)^2 + y^2 + z^2)^(1/2) + (2*u)/((u + x - 1)^2 + y^2 + z^2)^(1/2);
    JC_ZH      = JC_J2p + JC_J2s + JC_J3p + JC_J3s + JC_J4p + JC_J4s + JC_J5p + JC_J5s + JC_J6p + JC_J6s;
    
    JCs_vec(kk) = JC_vanilla + JC_ZH;
end

end

