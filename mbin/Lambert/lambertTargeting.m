function [v1,v2] = lambertTargeting(r1, r2, tof, transferType, direction, mu, printAssumptions)
%%% Description
%       Solves lamberts problem
% --------------------------------------------------------------
%%% Inputs
%       r1               - Initial position vector (km) [3x1]
%       r2               - Target position vector (km) [3x1]
%       tof              - Time of flight between positions (seconds)
%       transferType     - Type 1,2,3,4,5... (integer ... 1 and 2 are for
%                           0-rev solutions; >2 corresponds to one or more
%                           complete revolutions
%       direction        - +1 (prograde) or -1 (retrograde)
%       mu               - gravitational parameter of primary body (km^3/s^2)
%       printAssumptions - print out assumptions (0-no, 1-yes)
% --------------------------------------------------------------
%%% Outputs
%       v1 - Required velocity at initial position
%       v2 - Resulting velocity at target position
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
%%% Setting clock
ticLambert = tic;

%%% Position magnitudes
R1 = norm(r1);
R2 = norm(r2);

%%% Change in true anomaly
% delta_ta = acos(dot(r1,r2)/(norm(r1)*norm(r2))); % rad
if printAssumptions == 1
	warning('Making asumption that planets are in ecliptic')
end
ta1 = atan2(r1(2),r1(1));
ta2 = atan2(r2(2),r2(1));
delta_ta = ta2 - ta1;

%%% Making sure 0 < delta_ta < 2pi
if delta_ta > 2*pi
    while delta_ta > 2*pi
        delta_ta = delta_ta - 2*pi;
    end
elseif delta_ta < 0
    while delta_ta < 0
        delta_ta = delta_ta + 2*pi;
    end
end

cos_ta = dot(r1,r2) / (R1*R2);
if printAssumptions == 1
    warning('Does it work to handle the cos-ta like this for retrograde?')
end
A = direction * sqrt(R1*R2*(1 + cos_ta));

%%% Whether or not trajectory is possible
if delta_ta == 0 || A == 0
    warning('Trajectory cannot be computed')
    return
end

%%% Determining number of revolutions
if transferType <= 2
    nRevs = 0;
else
    nRevs = floor((transferType-1)/2);
end

%%% Initializing
if nRevs == 0 % If no complete revolutions
    psi = 0;
    psi_hi  = 4*pi*pi;
    psi_lo = -4*pi;
    
else % If 1 or more complete revolutions
    %%% Set initial bounds for psi
    psi_hi  = 4*pi*pi*((nRevs+1)^2);
    psi_lo = 4*pi*pi*(nRevs^2);
    
    %%% Create test-vector of psi values
    psi_test = linspace(psi_lo + 1e-3, psi_hi - 1e-3, 1000);
    tof_test = [];
    
    %%% Loop through psi_test values to find the minimum TOF solution
    for kk = 1:length(psi_test)
        psi_kk = psi_test(kk);
        
        if psi_kk > (1e-6)
            c2_test = (1.0 - cos(sqrt(psi_kk)))/psi_kk;
            c3_test = (sqrt(psi_kk) - sin(sqrt(psi_kk)))/sqrt(psi_kk^3);
        elseif psi_kk < (-1e-6)
            c2_test = (1.0 - cosh(sqrt(-psi_kk)))/psi_kk;
            c3_test = (sinh(sqrt(-psi_kk)) - sqrt(-psi_kk))/sqrt((-psi_kk)^3);
        else
            c2_test = 1/2;
            c3_test = 1/6;
        end
        
        y_test = R1 + R2 + A*(psi_kk*c3_test - 1)/sqrt(c2_test);
        chi_test = sqrt(y_test/c2_test);
        
        %%% Append TOF
        tof_test = [tof_test; ((chi_test^3)*c3_test + A*sqrt(y_test)) / sqrt(mu)];
    end
    
    %%% Find minimum TOF solution
    tof_min = min(tof_test);
    tof_min_index = find(tof_test == tof_min);
    psi_bound = psi_test(tof_min_index);
    
    %%% Use minimum TOF solution to form bounds on psi based on transfer type
    if rem(transferType,2) == 0 % If 'even' type ... positive slope
        psi_hi = 4*pi*pi*((nRevs+1)^2);
        psi_lo = psi_bound;
    elseif rem(transferType,2) == 1 % if 'odd' type ... negative slope
        psi_hi = psi_bound;
        psi_lo = 4*pi*pi*(nRevs^2);
    end
    
    %%% Set new initialized value of psi
    psi = (psi_hi + psi_lo)/2;
        
end

%%% Initialize other parameters
c2 = 1/2;
c3 = 1/6;
tof_est = 1e10;
timeFlag = 0;
N = 0.8;

while abs(tof_est - tof) > (1e-6) % While unconverged on TOF
    %%% Send a warning if it's taking a while
    if (toc(ticLambert) > 20) && (timeFlag == 0)
        warning('This is taking a while...')
        timeFlag = 1;
    end
    
    y = R1 + R2 + A*(psi*c3 - 1)/sqrt(c2);
    
    %%% If necessary, adjusting psi so that y is positive
    if (A > 0.0 && y < 0.0)
        while y < 0.0
            psi = (N/c3)*(1 - (sqrt(c2)/A)*(R1+R2));
            
            if (psi > 1e-6)
                c2 = (1.0 - cos(sqrt(psi)))/psi;
                c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi^3);

            elseif (psi < (-1e-6))
                c2 = (1.0 - cosh(sqrt(-psi)))/psi;
                c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/sqrt((-psi)^3);

            else
                c2 = 1/2;
                c3 = 1/6;
            end
            y = R1 + R2 + A*(psi*c3 - 1)/sqrt(c2);
        end
    end
    
    %%% Form new TOF estimate
    chi = sqrt(y/c2);
    tof_est = ((chi^3)*c3 + A*sqrt(y)) / sqrt(mu);
    
    %%% Use new TOF estimate to set new bounds on psi
    if (rem(transferType,2) == 0) || nRevs == 0
        if (tof_est <= tof)
            psi_lo = psi;
        else
            psi_hi = psi;
        end
    elseif (rem(transferType,2) == 1)
        if (tof_est >= tof)
            psi_lo = psi;
        else
            psi_hi = psi;
        end
    end
    
    %%% Set new psi value
    psi = (psi_hi + psi_lo)/2;
    
    %%% If elliptical
    if (psi > 1e-6)
        c2 = (1.0 - cos(sqrt(psi)))/psi;
        c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi^3);
    
    %%% If hyperbolic
    elseif (psi < (-1e-6))
        c2 = (1.0 - cosh(sqrt(-psi)))/psi;
        c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/sqrt((-psi)^3);
    
    %%% If borderline (parabolic)
    else
        c2 = 1/2;
        c3 = 1/6;
    end
    
end

%%% Calculate final velocities
f = 1 - y/R1;
g_dot = 1 - y/R2;
g = A*sqrt(y/mu);

v1 = (r2 - r1.*f)./g;
v2 = (r2.*g_dot - r1)./g;
    
end



