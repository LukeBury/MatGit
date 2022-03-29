%This function takes in a state vector composed of the launch date and the
%arrival date. It uses the dates in a Lambert function and returns the cost function f. fmincon will try to minimize
%the value of f. Users wil need to isert their own algorithms for
%determining the planetary positions and Lambert solutions. 
function f = fmincon_example_function(X);

global MU AU Sun_mu

launch  = X(1);
arrival = X(2);
TOFd    = arrival - launch;

%Insert your own Meeus algorithm to compute the positions of Earth at
%launch and Mars at Arrival


% % % [Vdepart, Varrive, Rdepart, Rarrive, TOF, Vplanet1, Vplanet2] = lambert_extend('earth', 'mars', launch, arrival);
[Rdepart, Vplanet1] = JulianDate_2_HeliocentricState(launch, 'Earth');
[Rarrive, Vplanet2] = JulianDate_2_HeliocentricState(launch, 'Mars');
[Vdepart, Varrive, ~] = lambertSolver(Rdepart, Rarrive, TOFd*86400, 0, 0, Sun_mu);

VLaunch_out = Vdepart - Vplanet1;
VLaunch  = norm(VLaunch_out);
C3Launch = VLaunch^2;

Vinf_in_vec = Varrive - Vplanet2;
Vinf_in     = norm(Vinf_in_vec);


%I've included several Cost functions. Choose your optimality criteria, set
%your weights and see what happens.

% %Minimize C3 only
% f = C3Launch;
% 
% 
% %Minimize C3 and V-infinity
% w1 = 1;
% w2 = 2;
% f = w1+C3Launch + w2*Vinf_in;
% 
% %Minimize C3/Arrival Vinfinity and TOF
% w1 = 10;
% w2 = 5;
% w3 = 1;
% f = w1+C3Launch + w2*Vinf_in + w3*TOFd;

%Try to normalize the parameters based on the maximum allowable
C3max   = 25;
Vinfmax = 7;
TOFdmax = 250;
f = C3Launch/C3max + Vinf_in/Vinfmax + TOFd/TOFdmax;
