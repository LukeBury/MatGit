function [r,v] = pleph3o(planet,jd)
AU=149597870.691;   %km
MUsun = 1.32712440017987E11;      %km3/s2  
tol=1e-10;
T = (jd - 2451545)/36525;



% -----------------------3rd Order Coefficients--------------------------
% Meeus,Jean, Astronomical Algorithms, 2nd ed.

switch planet
    case 1      % Mercury
        ac = [0.387098310 0 0 0]*AU;         % [km]
        ec = [0.20563175 0.000020407 -0.0000000283 -0.00000000018];
        ic = [7.004986 -0.0059516 0.0000008 0.000000043];       % [deg]
        oc = [48.330893 -0.1254227 -0.00008833 -0.000000200];       % [deg]
        wlc = [77.456119 0.1588643 -0.00001342 -0.000000007];       % [deg]
        Lc = [252.250906 149472.6746358 -0.00000536 0.000000002];        % [deg]


    case 2      % Venus
        ac = [0.72332982 0 0 0] * AU;                          % [km]
        ec = [0.00677192 -0.000047765 0.0000000981 0.00000000046];
        ic = [3.394662 -0.0008568 -0.00003244 0.000000009];             % [deg]
        oc = [76.679920 -0.2780134 -0.00014257 -0.000000164];           % [deg]
        wlc =[131.563703 0.0048646 -0.00138467 -0.000005695];           % [deg]
        Lc = [181.979801 58517.8156760 0.00000165 -0.000000002];        % [deg]

    case 3      % Earth
        ac = [1.000001018 0 0 0] * AU;                         % [km]
        ec = [0.01670863 -0.000042037 -0.0000001267 0.00000000014];
        ic = [0.0 0.0130548 -0.00000931 -0.000000034];                  % [deg]
        oc = [174.873176 -0.2410908 0.00004262 0.000000001];           % [deg]
        wlc = [102.937348 0.3225654 0.00014799 -0.0000000039];             % [deg]
        Lc = [100.466457 35999.3728565 -0.00000568 -0.000000001];                % [deg]
    
    case 4      % Mars
        ac = [1.523679342 0 0 0] * AU;                         % [km]
        ec = [0.09340065 0.000090484 -0.0000000806 -0.00000000025];
        ic = [1.849726 -0.0081477 -0.00002255 -0.000000029];            % [deg]
        oc = [49.558093 -0.2950250 -0.00064048 -0.000001964];           % [deg]
        wlc = [336.060234 0.4439016 -0.00017313 0.000000518];            % [deg]
        Lc = [355.433 19140.2993039 0.00000262 -0.000000003];        % [deg]
    
    case 5      % Jupiter
        ac = [5.202603209 0.0000001913 0 0] * AU;              % [km]
        ec = [0.04849793 0.000163225 -0.0000004714 -0.00000000201];
        ic = [1.303267 -0.0019877 0.0000332 0.000000097];              % [deg]
        oc = [100.464407 0.1767232 0.000907 -0.000007272];            % [deg]
        wlc = [14.331207 0.2155209 0.00072211 -0.000004485];             % [deg]
        Lc = [34.351519 3034.9056606 -0.00008501 0.000000016];          % [deg]

    case 6      % Saturn
        ac = [9.554909192 -0.000002139 0 0] * AU;                 % [km]
        ec = [0.05554814 -0.000346641 -0.0000006436 0.0000000034];
        ic = [2.488879 0.0025514 -0.00004906 0.000000017];              % [deg]
        oc = [113.665503 -0.2566722 -0.00018399 0.00000048];           % [deg]
        wlc = [93.057237 0.5665415 0.0005285 0.000004912];              % [deg]
        Lc = [50.077444 1222.1138488 0.00021004 -0.000000046];          % [deg]

    case 7      % Uranus
        ac = [19.218446062 -0.0000000372 0.00000000098 0]*AU;      % [km]
        ec = [0.04638122 -0.000027293 0.0000000789 0.00000000024];
        ic = [0.773197 -0.0016869 0.00000349 0.000000016];       % [deg]
        oc = [74.005957 0.0741431 0.00040539 0.000000119];       % [deg]
        wlc = [173.005291 0.0893212 -0.00009470 0.000000414];       % [deg]
        Lc = [314.055005 428.4669983 -0.00000486 0.000000006];       % [deg]

    case 8       % Neptune
        ac = [30.110386869 -0.0000001663 0.00000000069 0]*AU;      % [km]
        ec = [0.00945575 0.000006033 0.0 -0.00000000005];
        ic = [1.769953 0.0002256 0.00000023 -0.000000000];       % [deg]
        oc = [131.784057 -0.0061651 -0.00000219 -0.000000078];       % [deg]
        wlc = [48.120276 0.0291866 0.0000761 0.000000000];       % [deg]
        Lc = [304.348665 218.4862002 0.00000059 -0.000000002];       % [deg]
        
    case 9      % Pluto
        ac = [39.48211675 -0.00031596 0 0] * AU;    % [km]
        ec = [0.24882730 0.00005170 0 0];
        ic = [17.14001206 0.00004818 0 0];                   % [deg]
        oc = [110.30393684 -0.01183482 0 0];                 % [deg]
        wlc = [224.06891629 -0.04062942 0 0];                 % [deg]
        Lc = [238.92903833 145.20780515 0 0];                % [deg]

    otherwise
        error('This planet is not defined for the pleph3o function');
end

% --------------------Keplerian elements------------------
a = ac(1) + ac(2)*T + ac(3)*T^2 + ac(4)*T^3;    % Semimajor axis
e = ec(1) + ec(2)*T + ec(3)*T^2 + ec(4)*T^3;    % Eccentricity
i = ic(1) + ic(2)*T + ic(3)*T^2 + ic(4)*T^3;    % Inclination
o = oc(1) + oc(2)*T + oc(3)*T^2 + oc(4)*T^3;    % Long. of asc. node
wl = wlc(1) + wlc(2)*T + wlc(3)*T^2 + wlc(4)*T^3;    % Long. of perihelion
L = Lc(1) + Lc(2)*T + Lc(3)*T^2 + Lc(4)*T^3;    % Mean longitude


L = L-sign(L)*360*floor(abs(L)/360);  %reduce to 0-360 deg range
o = o-sign(o)*360*floor(abs(o)/360);  %reduce to 0-360 deg range
wl = wl-sign(wl)*360*floor(abs(wl)/360);  %reduce to 0-360 deg range

L = L*pi/180;    %conv to rad
o = o*pi/180;    %conv to rad
wl = wl*pi/180;    %conv to rad

w = wl-o;    % Argument of perihelion
M = L-wl;    % Mean anomaly


   %Equation of center
Ccen = (2*e-e^3/4+5/96*e^5)*sin(M) + (5/4*e^2-11/24*e^4)*sin(2*M) + ...
       ((13/12)*e^3-43/64*e^5)*sin(3*M) + (103/96)*e^4*sin(4*M) + ...
       (1097/960)*e^5*sin(5*M); %rad 
   
ta = M+Ccen; % True anomaly


% E0= M+e*sin(M); %initialize E
% dE=Inf; %arbitrary large number
% while abs(dE)>tol
%    E=M+e*sin(E0); %rad
%    dE=E-E0;
%    E0=E;
% end
%  ta=2*atan(sqrt((1+e)/(1-e))*tan(E/2));


 o=o*180/pi;
 w=w*180/pi;
 ta=ta*180/pi;

[r,v] = orbel2rv(MUsun,a,e,i,o,w,ta);



