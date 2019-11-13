function [tout, yout] = ode78e(F, t0, tfinal, y0, tol, trace, evfcn, etol)
% [tout, yout] = ode78e(F, t0, tfinal, y0, tol, trace, evfcn, etol)
% USE: Runge-Kutta-Fehlberg 7-8 integrator algorithm. Backward and forward
%      integration capability. Event function can be specified that terminates
%      integration.
% IN : F         :function handle for equations of motion
% t0        :initial time  [scalar]
% tfinal    :final time [scalar]
% y0        :initial state [nx1]
% tol       :integrator absolute tolerance [scalar] {1e-6}
% trace     :debugging switch, 1/0 [boolean] {0}
% evfcn     :event function handle {N/A}
%            for event functions with parameter inputs, make sure to use the format:
%            evfcn = @(t,y) fname(t, y, par1,..., parN);
% OUT:
% tout      :output time array [mx1]
% yout      :output state array [mxn]
% 
% LOG
% 08/03/2015 % Brian D. Anderson  Added general event function, backwards integration 
%   Removed superfluous input.
% 07/20/2015 % Travis Swenson %   Added event function
%   Improved accuracy and speed of event trigger.
% 05/19/2001 % Marc Compere   %   Backward integration.
% 10/06/1999 % Jonathan Essen %   Original Code.
% 

% initialize continuation variable
integrating = 1;

% The Fehlberg coefficients:
alpha   = [ 2./27. 1/9 1/6 5/12 .5 5/6 1/6 2/3 1/3 1 0 1 ]';
beta    = [ ...
[2/27       0       0      0        0         0       0         0     0      0     0 0 0 ]
[1/36       1/12    0      0        0         0       0         0     0      0     0 0 0 ]
[1/24       0       1/8    0        0         0       0         0     0      0     0 0 0 ]
[5/12       0       -25/16 25/16    0         0       0         0     0      0     0 0 0 ]
[.05        0       0      .25      .2        0       0         0     0      0     0 0 0 ]
[-25/108    0       0      125/108  -65/27    125/54  0         0     0      0     0 0 0 ]
[31/300     0       0      0        61/225    -2/9    13/900    0     0      0     0 0 0 ]
[2          0       0      -53/6    704/45    -107/9  67/90     3     0      0     0 0 0 ]
[-91/108    0       0      23/108   -976/135  311/54  -19/60    17/6  -1/12  0     0 0 0 ]
[2383/4100  0       0      -341/164 4496/1025 -301/82 2133/4100 45/82 45/164 18/41 0 0 0 ]
[3/205      0       0      0        0         -6/41   -3/205    -3/41 3/41   6/41  0 0 0 ]
[-1777/4100 0       0      -341/164 4496/1025 -289/82 2193/4100 51/82 33/164 12/41 0 1 0 ]...
]';
chi     = [0 0 0 0 0 34/105 9/35 9/35 9/280 9/280 0 41/840 41/840 ]';
psi     = [1 0 0 0 0 0      0    0    0     0     1 -1     -1     ]';
pow     = 1 / 8;
if nargin < 8,      etol    = 1e-13;                    end
if nargin < 7,      event   = 0;        else event = 1; end
if nargin < 6,      trace   = 0;                        end
if nargin < 5,      tol     = 1.e-6;                    end

% Initializations
t       = t0;
hmax    = (tfinal - t) / 2.5;
% hmin    = (tfinal - t) / 800000000;  % Original Tweak
% hmin    = (tfinal - t) / 8000000000000;  % Travis' Tweak
hmin    = (tfinal - t) / 1e22;  % Brian's Tweak 
h       = (tfinal - t) / 100;  % Tweak this?
y       = y0(:);
f       = y * zeros(1,13);
tout    = t;
yout    = y.';
tau     = tol * max(norm(y, 'inf'), 1);

if trace
    %  clc, t, h, y
    clc, t, y
end

% evaluate event function at initial state
if event
    evf = feval(evfcn,t,y);
end

if tfinal > t0 %FORWARD INTEGRATION
    restep = 0;
    % The main loop
    while (t < tfinal) & (h >= hmin) & integrating
        if t + h > tfinal, h = tfinal - t; end
        
        % Compute the slopes
        f(:,1)  = feval(F,t,y);
        for j = 1: 12
            f(:, j+1)   = feval(F, t + alpha(j) * h, y + h * f * beta(:, j));
        end
        
        % Truncation error term
        gamma1  = h * 41 / 840 * f * psi;
        
        % Estimate the error and the acceptable error
        delta   = norm(gamma1, 'inf');
        tau     = tol*max(norm(y, 'inf'), 1.0);
        
        % Update the solution only if the error is acceptable
        if delta <= tau
            t       = t + h;
            y       = y + h * f * chi;
            tout    = [tout; t];
            yout    = [yout; y.'];
        end
        if trace
            %        home, t, h, y
            home, t, y
        end
        
        % determine event trigger
        restep = 0;
        if (event == 1)
            % compute event function, save previous value
            evf_prev    = evf;
            evf         = feval(evfcn, t, y);
            
            % determine sign change in event function
            if sign(evf) ~= sign(evf_prev) && sign(evf_prev) ~= 0
                % If we need to update algorithm
                if abs(evf) > etol
                    % Delete last step to zoom in and refine
                    tout(end)   = [];
                    yout(end,:) = [];
                    t           = t - h;
                    y           = y - h * f * chi;
                    evf         = evf_prev;
                    restep      = 1;
                else
                    integrating = 0;
                end
            end
        end
        
        % Update the step size
        if delta ~= 0.0 && restep == 0 %avoid adaptive step size when refining event
            h   = min([hmax, 0.8 * h * (tau / delta) ^ pow]);
        end
        
        % reduce step size when backtracking to find exact event trigger
        if restep == 1;
            hmax    = h / 5;
            h       = hmax;
        end
    end
    
    if (t < tfinal)&(event==0)
        disp('SINGULARITY LIKELY.')
        t
    end
else %BACKWARD INTEGRATION
    restep = 0;
    % The main loop
    while (t > tfinal) & (h <= hmin) & integrating
        if t + h < tfinal, h = tfinal - t; end
        
        % Compute the slopes
        f(:,1)  = feval(F,t,y);
        for j = 1: 12
            f(:, j+1)   = feval(F, t + alpha(j) * h, y + h * f * beta(:, j));
        end
        
        % Truncation error term
        gamma1  = h * 41 / 840 * f * psi;
        
        % Estimate the error and the acceptable error
        delta   = norm(gamma1, 'inf');
        tau     = tol*max(norm(y, 'inf'), 1.0);
        
        % Update the solution only if the error is acceptable
        if delta <= tau
            t       = t + h;
            y       = y + h * f * chi;
            tout    = [tout; t];
            yout    = [yout; y.'];
        end
        if trace
            %        home, t, h, y
            home, t, y
        end
        
        % determine event trigger
        restep = 0;
        if (event == 1)
            % compute event function, save previous value
            evf_prev    = evf;
            evf         = feval(evfcn, t, y);
            
            % determine sign change in event function
            if sign(evf) ~= sign(evf_prev) && sign(evf_prev) ~= 0
                % If we need to update algorithm
                if abs(evf) > etol
                    % Delete last step to zoom in and refine
                    tout(end)   = [];
                    yout(end,:) = [];
                    t           = t - h;
                    y           = y - h * f * chi;
                    evf         = evf_prev;
                    restep      = 1;
                else
                    integrating = 0;
                end
            end
        end
        
        % Update the step size
        if delta ~= 0.0 && restep == 0 %avoid adaptive step size when refining event
            h   = max([hmax, 0.8 * h * (tau / delta) ^ pow]);
        end
        
        % reduce step size when backtracking to find exact event trigger
        if restep == 1;
            hmax    = h / 5;
            h       = hmax;
        end
    end
    
    if (t > tfinal)&(event==0)
        disp('SINGULARITY LIKELY.')
        t
    end
end


if event && (abs(evf) > etol)
    abs(evf)
    warning('event tolerance not achieved with current minimum stepsize')
end

end
