function [r,v] = pleph2o(planet,jd)
AU=149597870.691;   %km
MUsun = 1.32712440017987E11;      %km3/s2  
tol=1e-8;
T = (jd - 2451545)/36525;
if (jd>2341972 && jd<2524593) %ephemeris valid for 1800 - 2050 but using 1700 -2200
switch planet
    case 1      % Mercury
        ac = [0.38709927 0.00000037 ]*AU;         % [km]
        ec = [0.20563593 0.00001906 ];
        ic = [7.00497902 -0.00594749 ];       % [deg]
        oc = [48.33076593 -0.12534081 ];       % [deg]
        wlc = [77.45779628 0.16047689 ];       % [deg]
        Lc = [252.25032350 149472.67411175];        % [deg]


    case 2      % Venus
        ac = [0.72333566 0.00000390] * AU;                          % [km]
        ec = [0.00677672 -0.00004107];
        ic = [3.39467605 -0.00078890 ];             % [deg]
        oc = [76.67984255 -0.27769418 ];           % [deg]
        wlc =[131.60246718 0.00268329 ];           % [deg]
        Lc = [181.97909950 58517.81538729 ];        % [deg]

    case 3      % Earth
        ac = [1.00000261 0.00000562] * AU;                         % [km]
        ec = [0.01671123 -0.00004392 ];
        ic = [-0.00001531 -0.01294668 ];                  % [deg]
        oc = [    0.0         0.0    ];           % [deg]
        wlc = [102.93768193 0.32327364 ];             % [deg]
        Lc = [100.46457166 35999.37244981 ];                % [deg]
    
    case 4      % Mars
        ac = [1.523671034 0.00001847 ] * AU;                         % [km]
        ec = [0.09339410 0.00007882 ];
        ic = [1.84969142 -0.00813131 ];            % [deg]
        oc = [49.55953891 -0.29257343 ];           % [deg]
        wlc = [-23.94362959 0.44441088 ];            % [deg]
        Lc = [-4.55343205 19140.30268499 ];        % [deg]
    
    case 5      % Jupiter
        ac = [5.20288700 -0.00011607] * AU;              % [km]
        ec = [0.04838624 -0.00013253 ];
        ic = [1.30439695 -0.00183714 ];              % [deg]
        oc = [100.47390909 0.20469106 ];            % [deg]
        wlc = [14.72847983 0.21252668 ];             % [deg]
        Lc = [34.39644051 3034.74612775 ];          % [deg]

    case 6      % Saturn
        ac = [9.53667594 -0.00125060 ] * AU;                 % [km]
        ec = [0.05386179 -0.00050991 ];
        ic = [2.48599187 0.00193609];              % [deg]
        oc = [113.66242448 -0.28867794 ];           % [deg]
        wlc = [92.59887831 -0.41897216 ];              % [deg]
        Lc = [49.95424423 1222.49362201 ];          % [deg]

    case 7      % Uranus
        ac = [19.18916464 -0.00196176]*AU;      % [km]
        ec = [0.04725744 -0.00004397 ];
        ic = [0.77263783 -0.00242939 ];       % [deg]
        oc = [74.01692503 0.04240589];       % [deg]
        wlc = [170.95427630 0.40805281 ];       % [deg]
        Lc = [313.23810451 428.48202785 ];       % [deg]

    case 8       % Neptune
        ac = [30.06992276 0.00026291 ]*AU;      % [km]
        ec = [0.00859048 0.00005105 ];
        ic = [1.77004347 0.00035372 ];       % [deg]
        oc = [131.78422574 -0.00508664 ];       % [deg]
        wlc = [44.96476227 -0.32241464 ];       % [deg]
        Lc = [-55.12002969 218.45945325 ];       % [deg]
        
    case 9      % Pluto
        ac = [39.48211675 -0.00031596 ] * AU;    % [km]
        ec = [0.24882730 0.00005170 ];
        ic = [17.14001206 0.00004818 ];                   % [deg]
        oc = [110.30393684 -0.01183482 ];                 % [deg]
        wlc = [224.06891629 -0.04062942];                 % [deg]
        Lc = [238.92903833 145.20780515 ];                % [deg]

    otherwise
        error('This planet is not defined for the pleph2o function');
end

    a = ac(1) + ac(2)*T ;    % Semimajor axis
    e = ec(1) + ec(2)*T ;    % Eccentricity
    I = ic(1) + ic(2)*T ;    % Inclination
    o = oc(1) + oc(2)*T ;    % Long. of asc. node
    wl = wlc(1) + wlc(2)*T ;    % Long. of perihelion
    L = Lc(1) + Lc(2)*T ;    % Mean longitude

else %ephemeris valid for 3000 BC - 3000 AD
    
switch planet
    case 1      % Mercury
        ac = [0.38709843 0.000000 ]*AU;         % [km]
        ec = [0.20563661 0.00002123 ];
        ic = [7.00559432 -0.00590158 ];       % [deg]
        oc = [48.33961819 -0.12214182 ];      % [deg]
        wlc =[77.45771895 0.15940013 ];       % [deg]
        Lc = [252.25166724 149472.67486623];    % [deg]


    case 2      % Venus
        ac = [0.72332102 -0.00000026] * AU;                          % [km]
        ec = [0.00676399 -0.00005107];
        ic = [3.39777545 0.00043494 ];             % [deg]
        oc = [76.67261496 -0.27274174 ];           % [deg]
        wlc =[131.76755713 0.05679648 ];           % [deg]
        Lc = [181.9797085 58517.8156026 ];        % [deg]

    case 3      % Earth
        ac = [1.00000018 -0.00000003] * AU;                         % [km]
        ec = [0.01673163 -0.00003661 ];
        ic = [-0.00054346 -0.01337178 ];                  % [deg]
        oc = [-5.11260389 -0.24123856 ];           % [deg]
        wlc =[102.93005885 0.3179526 ];             % [deg]
        Lc = [100.46691572 35999.37306329];                % [deg]
    
    case 4      % Mars
        ac = [1.52371243 0.00000097 ] * AU;                         % [km]
        ec = [0.09336511 0.00009149 ];
        ic = [1.85181869 -0.00724757 ];            % [deg]
        oc = [49.71320984 -0.26852431 ];           % [deg]
        wlc =[-23.91744784 0.45223625 ];            % [deg]
        Lc = [-4.56813164 19140.29934243 ];        % [deg]
    
    case 5      % Jupiter
        ac = [5.20248019 -0.00002864] * AU;              % [km]
        ec = [0.0485359 0.00018026 ];
        ic = [1.29861416 -0.00322699 ];              % [deg]
        oc = [100.29282654 0.13024619 ];            % [deg]
        wlc =[14.27495244 0.18199196 ];             % [deg]
        Lc = [34.33479152 3034.90371757 ];          % [deg]
        
        b=-0.00012452;
        c=0.0606406;
        s=-0.35635438;
        f=38.35125;

    case 6      % Saturn
        ac = [9.54149883 -0.00003065 ] * AU;                 % [km]
        ec = [0.05550825 -0.00032044 ];
        ic = [2.49424102 0.00451969];              % [deg]
        oc = [113.63998702 -0.25015002 ];           % [deg]
        wlc =[92.86136063 0.54179478 ];              % [deg]
        Lc = [50.07571329 1222.11494724 ];          % [deg]
        
        b=0.00025899;
        c=-0.13434469;
        s=0.87320147;
        f=38.35125;

    case 7      % Uranus
        ac = [19.18797948 -0.00020455]*AU;      % [km]
        ec = [0.0468574 -0.0000155 ];
        ic = [0.77298127 -0.00180155 ];       % [deg]
        oc = [73.96250215 0.05739699];       % [deg]
        wlc =[172.43404441 0.09266985 ];       % [deg]
        Lc = [314.20276625 428.49512595 ];       % [deg]
        
        b=0.00058331;
        c=-0.97731848;
        s=0.17689245;
        f=7.67025;

    case 8       % Neptune
        ac = [30.06952752 0.00006447 ]*AU;      % [km]
        ec = [0.00895439 0.00000818 ];
        ic = [1.7700552 0.000224 ];       % [deg]
        oc = [131.78635853 -0.00606302 ];       % [deg]
        wlc =[46.68158724 0.01009938 ];       % [deg]
        Lc = [304.22289287 218.46515314];       % [deg]
        
        b=-0.00041348;
        c=0.68346318;
        s=-0.10162547;
        f=7.67025;
        
    case 9      % Pluto
        ac = [39.48686035 0.00449751 ] * AU;    % [km]
        ec = [0.24885238 0.00006016 ];
        ic = [17.1410426 0.00000501 ];                   % [deg]
        oc = [110.30167986 -0.00809981 ];                 % [deg]
        wlc =[224.09702598 -0.00968827];                 % [deg]
        Lc = [238.96535011 145.18042903 ];                % [deg]
        
        b=-0.01262724;
        c=0;
        s=0;
        f=0;
    otherwise
        error('This planet is not defined for the pleph2o function');
end
    if planet<5 %check for inner planet
        a = ac(1) + ac(2)*T ;    % Semimajor axis
        e = ec(1) + ec(2)*T ;    % Eccentricity
        I = ic(1) + ic(2)*T ;    % Inclination
        o = oc(1) + oc(2)*T ;    % Long. of asc. node
        wl = wlc(1) + wlc(2)*T ;    % Long. of perihelion
        L = Lc(1) + Lc(2)*T ;    % Mean longitude
    else
        a = ac(1) + ac(2)*T ;    % Semimajor axis
        e = ec(1) + ec(2)*T ;    % Eccentricity
        I = ic(1) + ic(2)*T ;    % Inclination
        o = oc(1) + oc(2)*T ;    % Long. of asc. node
        wl = wlc(1) + wlc(2)*T ;    % Long. of perihelion
        L = Lc(1) + Lc(2)*T + b*T*T + c*cosd(f*T) + s*sind(f*T);% Mean longitude
    end

end
% --------------------Keplerian elements------------------


L = mod(L,360);  %reduce to 0-360 deg range
o = mod(o,360);  %reduce to 0-360 deg range
wl = mod(wl,360);  %reduce to 0-360 deg range


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

% if M>pi
%    M=M-2*pi; 
% end
% % E0= M+e*sin(M); %initialize E
% E0=atan2(sin(M),cos(M)-e); %initialize E with formula for low eccentricty
% dE=1E10; %arbitrary large number
% 
% while abs(dE)>tol
%    E=M+e*sin(E0); %rad
%    dE=E-E0;
%    E0=E;
% end
%  ta=2*atan(sqrt((1+e)/(1-e))*tan(E/2));

[r,v] = orbel2rv(MUsun,a,e,I,o*180/pi,w*180/pi,ta*180/pi);
end