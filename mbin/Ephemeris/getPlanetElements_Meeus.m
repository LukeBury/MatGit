function [L, a, e, i, Omega, Pi, w, M, ta] = getPlanetElements_Meeus(JDE, bodyName, angleUnit)
%%% Description
%       Returns orbital elements for desired planet at specified time
% --------------------------------------------------------------
%%% Inputs
%       JDE      - Julian date of desired time        ... [1x1]
%       bodyName - Name of desired planet, ex 'Earth' ... [string]
%       angleUnit - Either 'radians' or 'degrees'     ... [string]
% --------------------------------------------------------------
%%% Outputs
%       L - Mean longitude of planet               ... (angleUnit)
%       a - Semimajor axis                         ... (km)
%       e - Eccentricity                           ...
%       i - Inclination                            ... (angleUnit)
%       Omega - Longitude of the ascending node    ... (angleUnit)
%       Pi - Longitude of perihelion (Pi = Om + w) ... (angleUnit)
%       w - Argument of periapsis                  ... (angleUnit)
%       M - Mean anomaly                           ... (angleUnit)
%       v - True anomaly                           ... (angleUnit)
% ===============================================================
%%% Compute T in Julian centuries of 36525 ephemeris days from the epoch
%%% J2000.0 = JDE 2451545
T = (JDE - 2451545)/36525;

%%% AU definition
AU_km = 1.49597870700e8;     % km
% -------------------------------------------------
% Grabbing coefficients for desired planet
% -------------------------------------------------
%%% Coefficients for the orbital elements of Venus
if isequal(bodyName, 'Venus')
    L0_deg     = [181.979801, 58517.8156760, 0.00000165, -0.000000002];
    a0_AU      = [0.72332982, 0, 0, 0];
    e0         = [0.00677188, -0.000047766, 0.0000000975, 0.00000000044];
    i0_deg     = [3.394662, -0.0008568, -0.00003244, 0.000000010];
    Omega0_deg = [76.679920, -0.2780080, -0.00014256, -0.000000198];
    Pi0_deg    = [131.563707, 0.0048646, -0.00138232, -0.000005332];
    
%%% Coefficients for the orbital elements of Earth
elseif isequal(bodyName, 'Earth')
    L0_deg     = [100.466449, 35999.3728519, -0.00000568, 0.0];
    a0_AU      = [1.000001018, 0, 0, 0];
    e0         = [0.01670862, -0.000042037, -0.0000001236, 0.00000000004];
    i0_deg     = [0, 0.0130546, -0.00000931, -0.000000034];
    Omega0_deg = [174.873174, -0.2410908, 0.00004067, -0.000001327];
    Pi0_deg    = [102.937348, 0.3225557, 0.00015026, 0.000000478];
    
%%% Coefficients for the orbital elements of Mars
elseif isequal(bodyName, 'Mars')
    L0_deg     = [355.433275, 19140.2993313, 0.00000261, -0.000000003];
    a0_AU      = [1.523679342, 0, 0, 0];
    e0         = [0.09340062, 0.000090483, -0.0000000806, -0.00000000035];
    i0_deg     = [1.849726, -0.0081479, -0.00002255, -0.000000027];
    Omega0_deg = [49.558093, -0.2949846, -0.00063993, -0.000002143];
    Pi0_deg    = [336.060234, 0.4438898, -0.00017321, 0.000000300];

%%% Coefficients for the orbital elements of Jupiter
elseif isequal(bodyName, 'Jupiter')
    L0_deg     = [34.351484, 3034.9056746, -0.00008501, 0.000000004];
    a0_AU      = [5.202603191, 0.0000001913, 0, 0];
    e0         = [0.04849485, 0.000163244, -0.0000004719, -0.00000000197];
    i0_deg     = [1.303270, -0.0019872, 0.00003318, 0.000000092];
    Omega0_deg = [100.464441, 0.1766828, 0.00090387, -0.000007032];
    Pi0_deg    = [14.331309, 0.2155525, 0.00072252, -0.000004590];

%%% Coefficients for the orbital elements of Saturn
elseif isequal(bodyName, 'Saturn')
    L0_deg     = [50.077471, 1222.1137943, 0.00021004, -0.000000019];
    a0_AU      = [9.554909596, -0.0000021389, 0, 0];
    e0         = [0.05550862, -0.000346818, -0.0000006456, 0.00000000338];
    i0_deg     = [2.488878, 0.0025515, -0.00004903, 0.000000018];
    Omega0_deg = [113.665524, -0.2566649, -0.00018345, 0.000000357];
    Pi0_deg    = [93.056787, 0.5665496, 0.00052809, 0.000004882];

%%% Coefficients for the orbital elements of Uranus
elseif isequal(bodyName, 'Uranus')
    L0_deg     = [314.055005, 429.8640561, 0.00030434, 0.000000026];
    a0_AU      = [19.21844606, -0.0000000372, 0.00000000098, 0.0];
    e0         = [0.04629590, -0.000027337, 0.0000000790, 0.00000000025];
    i0_deg     = [0.773196, 0.0007744, 0.00003749, -0.000000092];
    Omega0_deg = [74.005947, 0.5211258, 0.00133982, 0.000018516];
    Pi0_deg    = [173.005159, 1.4863784, 0.0021450, 0.000000433];

%%% Coefficients for the orbital elements of Neptune
elseif isequal(bodyName, 'Neptune')
    L0_deg     = [304.348665, 219.8833092, 0.00030926, 0.000000018];
    a0_AU      = [30.110386869, -0.0000001663, 0.00000000069, 0.0];
    e0         = [0.00898809, 0.000006408, -0.0000000008, -0.00000000005];
    i0_deg     = [1.769952, -0.0093082, -0.00000708, 0.000000028];
    Omega0_deg = [131.784057, 1.1022057, 0.00026006, -0.000000636];
    Pi0_deg    = [48.123691, 1.4262677, 0.00037918, -0.000000003];

%%% Coefficients for the orbital elements of Pluto
elseif isequal(bodyName, 'Pluto')
    L0_deg     = [238.92903833, 145.20780515, 0.0, 0.0];
    a0_AU      = [39.48211675, -0.00031596, 0.0, 0.0];
    e0         = [0.24882730, 0.00005170, 0.0, 0.0];
    i0_deg     = [17.14001206, 0.00004818, 0.0, 0.0];
    Omega0_deg = [110.30393684, -0.01183482, 0.0, 0.0];
    Pi0_deg    = [224.06891629, -0.04062942, 0.0, 0.0];

%%% Catching bad entry
else
    warning('Incorrect body name provided')
    return
end

% -------------------------------------------------
% Calculating Desired Elements
% -------------------------------------------------
L     = L0_deg(1) + L0_deg(2)*T + L0_deg(3)*T^2 + L0_deg(4)*T^3;                 % deg
a_AU  = a0_AU(1) + a0_AU(2)*T + a0_AU(3)*T^2 + a0_AU(4)*T^3;                     % deg
a     = a_AU * AU_km;                                                            % km
e     = e0(1) + e0(2)*T + e0(3)*T^2 + e0(4)*T^3;                                 %
i     = i0_deg(1) + i0_deg(2)*T + i0_deg(3)*T^2 + i0_deg(4)*T^3;                 % deg
Omega = Omega0_deg(1) + Omega0_deg(2)*T + Omega0_deg(3)*T^2 + Omega0_deg(4)*T^3; % deg
Pi    = Pi0_deg(1) + Pi0_deg(2)*T + Pi0_deg(3)*T^2 + Pi0_deg(4)*T^3;             % deg

w = Pi - Omega;                                                                  % deg
M = L - Pi;                                                                      % deg

Ccen = (2*e - (e^3)/4 + (5/96)*e^5)*sind(M) + ((5/4)*e^2 - (11/24)*e^4)*sind(2*M) + ...
    ((13/12)*e^3 - (43/64)*e^5)*sind(3*M) + (103/96)*(e^4)*sind(4*M) + (1097/960)*(e^5)*sind(5*M); % rad
Ccen = Ccen * 180/pi;                                                            % deg

ta = M + Ccen;                                                                   % deg

% -------------------------------------------------
% Setting angular units
% -------------------------------------------------
if isequal(angleUnit, 'degrees')
    return
elseif isequal(angleUnit,'radians')
    deg2rad = pi/180;
    
    L       = L*deg2rad;
    i       = i*deg2rad;
    Omega   = Omega*deg2rad;
    Pi      = Pi*deg2rad;
    w       = w*deg2rad;
    M       = M*deg2rad;
    ta      = ta*deg2rad;
end
end









