%%% Demonstrate numerical integration %%%

clear all; %close all

% ============= %
% --- TRUTH --- %
% ============= %

% An initial condition:
X0      = [0; 0];       % [m]
X0_dot  = [5; 10];      % [m/s]
X0_ddot = [0; -9.81];   % [m/s^2]

% Time bounds
t0      = 0;            % [s]
tf      = 2*10/9.81;    % [s]

% Fine time vector for truth
t = [t0 : 0.001 : tf];

% The analytic truth
X(1,:) = X0(1) + X0_dot(1)*t+ (1/2)*X0_ddot(1)*t.^2;
X(2,:) = X0(2) + X0_dot(2)*t+ (1/2)*X0_ddot(2)*t.^2;

% Plot the truth
plot(X(1,:),X(2,:),'-b')
hold on

% =========================== %
% --- EULER APPROXIMATION --- %
% =========================== %
clear X
X_ddot = X0_ddot;

% Define a time step
dt = 0.1;                   % }--- THE CHANGE IN INDEPENDENT VARIABLE           % }..................
                                                                                % }
% Integrate                                                                     % }
t = t0;                     % }--- THE INITIAL TIME                             % }
X(:,1) = X0;                % }--- THE INITIAL STATE (X and X_DOT)              % }
X_dot(:,1) = X0_dot;        % }-/                                               % }
i = 2;                                                                          % }
last = 0;                                                                       % }
while true                                                                      % }--- These are common to EVERY numerical integrator
    X(:,i) = X(:,i-1) + X_dot(:,i-1)*dt;        % }--- THE SYSTEM DYNAMICS      % }
    X_dot(:,i) = X_dot(:,i-1) + X_ddot*dt;      % }-/                           % }
                                                                                % }
    if last==1                                                                  % }
        break                                                                   % }
    end                                                                         % }
                                                                                % }
    t = t + dt;                                                                 % }
    if t > tf                                   % }--- THE STOPPING CONDITION   % }..................
        t = tf;
        last = 1;
    end
    i = i+1;
end
    
plot(X(1,:),X(2,:),'.-r')   