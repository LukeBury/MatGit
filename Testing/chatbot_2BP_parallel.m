
states = [1, 1, 1, 1e-7, 0, 0].*100000;

% Define gravitational parameter of Earth
mu = 3.986e14;

% Define time step for propagation
dt = 60;



% Loop over each state and propagate in parallel
parfor i = 1:size(states, 1)
    % Extract initial position and velocity
    r0 = states(i, 1:3);
    v0 = states(i, 4:6);
    
    % Set initial time and state for ode113
    t0 = 0;
    y0 = [r0, v0];
    
    % Propagate trajectory using ode113
    [t, y] = ode113(keplerianDynamics, [t0, dt], y0);
    
    % Extract position and velocity from solution
    r = y(end, 1:3);
    v = y(end, 4:6);
    
    % Store results in structure
    results(i).r = r;
    results(i).v = v;
end

% Define function for Keplerian dynamics
function dy = keplerianDynamics(t, y)
    % Extract position and velocity
    r = y(1:3);
    v = y(4:6);
    
    % Compute acceleration using Kepler's laws
    a = -mu * r / norm(r)^3;
    
    % Set derivative of state vector
    dy = [v, a];
end