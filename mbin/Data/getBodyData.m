function [bodies] = getBodyData(mbinPath)
% Values from Wikipedia and Vallado and MONTE and NAIF
% ------------------------------------------------------------------------
%%% Gravitational constant
% ------------------------------------------------------------------------
bodies.constants.G     = 6.6726e-20; % km^3 * kg^-1 * s^-2
bodies.constants.AU_km = 1.495978707e8; % km
% ------------------------------------------------------------------------
%%% System - Sun
% ------------------------------------------------------------------------
%%% Sun
bodies.sun.name  = 'sun';
bodies.sun.title = 'Sun';
bodies.sun.mass  = 1.9891e30; % kg
bodies.sun.u     = bodies.sun.mass * bodies.constants.G; % km^3 * s^-2
bodies.sun.color = [1,1,1]./255;
bodies.sun.R     = 695700; % km
bodies.sun.J2    = 2e-6;

% ------------------------------------------------------------------------
%%% System - Mercury
% ------------------------------------------------------------------------
%%% Mercury
bodies.mercury.name    = 'mercury';
bodies.mercury.title   = 'Mercury';
bodies.mercury.mass    = 3.3011e23; % kg
bodies.mercury.u       = bodies.mercury.mass * bodies.constants.G; % km^3 * s^-2
bodies.mercury.a       = 57909050; % km
bodies.mercury.color   = [1,1,1]./255;
bodies.mercury.R       = 2440; % km
bodies.mercury.meanMot = 8.266781726111516e-07; %2*pi/(87.969 * 86400);
bodies.mercury.MR      = 1.659594516188856e-07; %bodies.mercury.mass / (bodies.mercury.mass + bodies.sun.mass);

% ------------------------------------------------------------------------
%%% System - Venus
% ------------------------------------------------------------------------
%%% Venus
bodies.venus.name    = 'venus';
bodies.venus.title   = 'Venus';
bodies.venus.mass    = 4.8675e24; % kg
bodies.venus.u       = bodies.venus.mass*bodies.constants.G; % km^3 * s^-2
bodies.venus.a       = 108208000; % km
bodies.venus.color   = [1,1,1]./255;
bodies.venus.R       = 6052; % km
bodies.venus.meanMot = 1.245411223565392e-07; %2*pi/(583.92 * 86400);
bodies.venus.Tp      = 2*pi/bodies.venus.meanMot; % sec
bodies.venus.MR      = 2.447080633872110e-06; %bodies.venus.mass / (bodies.venus.mass + bodies.sun.mass);

% ------------------------------------------------------------------------
%%% System - Earth
% ------------------------------------------------------------------------
%%% Earth
bodies.earth.name    = 'earth';
bodies.earth.title   = 'Earth';
bodies.earth.color   = [36, 38, 235]./255;
bodies.earth.img     = imread([mbinPath,'/textures/earthSurfTex.jpg']);
bodies.earth.mass    = 5.9742e24; % kg
bodies.earth.u       = bodies.earth.mass * bodies.constants.G; % km^3 * s^-2
bodies.earth.a       = 149600000; % km
bodies.earth.R       = 6378; % km
bodies.earth.meanMot = 1.991063857259666e-07; %2*pi/(365.242189 * 86400);
bodies.earth.Tp      = 2*pi/bodies.earth.meanMot; %365.242189 * 86400; % sec
bodies.earth.MR      = 3.003459884736794e-06; %bodies.earth.mass / (bodies.earth.mass + bodies.sun.mass);
bodies.earth.J2      = 0.0010826269;
bodies.earth.J3      = -0.0000025323;
bodies.earth.J4      = -0.0000016204;

%%% Moon
bodies.moon.name    = 'moon';
bodies.moon.title   = 'Moon';
bodies.moon.color   = [168, 168, 168]./255;
bodies.moon.img     = imread([mbinPath,'/textures/moonSurfTex.jpg']);
bodies.moon.mass    = 7.347e22; % kg
bodies.moon.u       = bodies.moon.mass*bodies.constants.G; % km^3 * s^-2
bodies.moon.a       = 384748; % km
bodies.moon.R       = 1737; % km
bodies.moon.R_n     = bodies.moon.R / bodies.moon.a;
bodies.moon.meanMot = 2.661699527215069e-06; %2*pi/(27.321661*86400); % rad/s
bodies.moon.Tp      = 2*pi/bodies.moon.meanMot; % sec
bodies.moon.MR      = 0.012148480323827; %bodies.moon.mass / (bodies.moon.mass + bodies.earth.mass); % Mass ratio w/ primary
bodies.moon.J2      = 2.032366226455845e-04; %(9.08901807505999991e-05)*sqrt(5); % from MONTE ... sqrt(5) is normalizer (vallado pg 547)
bodies.moon.LyapTp  = 3.373237818307935;

% ------------------------------------------------------------------------
%%% System - Mars
% ------------------------------------------------------------------------
%%% Mars
bodies.mars.name    = 'mars';
bodies.mars.title   = 'Mars';
bodies.mars.color   = [224, 112, 99]./255;
bodies.mars.img     = imread([mbinPath,'/textures/marsSurfTex.jpg']);
bodies.mars.mass    = 6.39e23; % kg
bodies.mars.u       = bodies.mars.mass*bodies.constants.G; % km^3 * s^-2
bodies.mars.a       = 227939186; % km
bodies.mars.e       = 0.0935; % https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
bodies.mars.R       = 3390; % km
bodies.mars.meanMot = 1.058575972610999e-07; %2*pi/(686.980 * 86400);
bodies.mars.J2      = 0.00196045;
bodies.mars.J3      = 0.000036;

%%% Phobos
bodies.phobos.name    = 'phobos';
bodies.phobos.title   = 'Phobos';
bodies.phobos.mass    = 1.0659e16; % kg
bodies.phobos.u       = bodies.phobos.mass * bodies.constants.G; % km^3 * s^-2
bodies.phobos.color   = [176, 130, 122]./255;
bodies.phobos.a       = 9376; % km
bodies.phobos.R       = 11.2667; % km
bodies.phobos.R_n     = bodies.phobos.R / bodies.phobos.a;
bodies.phobos.meanMot = 2.280329864815889e-04; %2*pi/(0.31891023*86400); % rad/s
bodies.phobos.Tp      = 2*pi/bodies.phobos.meanMot; % sec
bodies.phobos.MR      = 1.668075089546147e-08; %bodies.phobos.mass / (bodies.phobos.mass + bodies.mars.mass);  % Mass ratio w/ primary
bodies.phobos.J2      = 0.10058; %(4.49807434353858004e-02)*sqrt(5); % from MONTE ... sqrt(5) is normalizer (vallado pg 547)

%%% Deimos
bodies.deimos.name    = 'deimos';
bodies.deimos.title   = 'Deimos';
bodies.deimos.mass    = 1.476e15; % kg
bodies.deimos.u       = bodies.deimos.mass * bodies.constants.G; % km^3 * s^-2
bodies.deimos.color   = [158, 129, 122]./255;
bodies.deimos.a       = 23463.2; % km
bodies.deimos.R       = 6.2; % km
bodies.deimos.R_n     = bodies.deimos.R / bodies.deimos.a;
bodies.deimos.meanMot = 5.757882198450546e-05; %2*pi/(1.263*86400); % rad/s
bodies.deimos.Tp      = 2*pi/bodies.deimos.meanMot; % sec
bodies.deimos.MR      = 2.309859149594128e-09; %bodies.deimos.mass / (bodies.deimos.mass + bodies.mars.mass);  % Mass ratio w/ primary

% ------------------------------------------------------------------------
%%% System - Jupiter
% ------------------------------------------------------------------------
%%% Jupiter
bodies.jupiter.name  = 'jupiter';
bodies.jupiter.title = 'Jupiter';
bodies.jupiter.mass  = 1.89819e27; % kg
bodies.jupiter.u     = bodies.jupiter.mass*bodies.constants.G; % km^3 * s^-2
bodies.jupiter.color = [235,214,173]./255;
bodies.jupiter.img   = imread([mbinPath,'/textures/jupiterSurfTex.jpg']);
bodies.jupiter.a     = 778.57e6; % km, https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
bodies.jupiter.R     = 69911; % km
bodies.jupiter.R_ZH  = 71492.0; % km   .... reference radius used for zonal harmonics
bodies.jupiter.J2    = 0.0146956247707; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op / MONTE - Jupiter solution 310
bodies.jupiter.J4    = -5.913138887463e-04; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op / MONTE - Jupiter solution 310
bodies.jupiter.J6    = 2.077510523749e-05; % MONTE - Jupiter solution 310

%%% Io
bodies.io.name    = 'io';
bodies.io.title   = 'Io';
bodies.io.mass    = 8.94e22; % kg
bodies.io.img     = imread([mbinPath,'/textures/ioSurfTex.jpg']);
bodies.io.color   = [1,1,1]./255;
bodies.io.a       = 421800; % km
bodies.io.R       = 1821.5; % km
bodies.io.R_n     = bodies.io.R / bodies.io.a;
bodies.io.meanMot = 4.110592399599715e-05; %2*pi/(1.769138*86400); % rad/s
bodies.io.MR      = 4.709528007310199e-05; %bodies.io.mass / (bodies.io.mass + bodies.jupiter.mass);  % Mass ratio w/ primary
bodies.io.J2      = 8.188668558770e-04; % from MONTE

%%% Europa
bodies.europa.name    = 'europa';
bodies.europa.title   = 'Europa';
bodies.europa.color   = [0, 1, 1];
bodies.europa.img     = imread([mbinPath,'/textures/europaSurfTex.jpg']);
% bodies.europa.mass    = 4.799e22; % kg
bodies.europa.a       = 671100; % km
bodies.europa.R       = 1560.8; % km
bodies.europa.R_n     = bodies.europa.R / bodies.europa.a;
bodies.europa.Tp      = 3.551181*86400; % sec
bodies.europa.meanMot = 2*pi / bodies.europa.Tp; % rad/s
% bodies.europa.MR      = 2.528133998624693e-05; %bodies.europa.mass / (bodies.europa.mass + bodies.jupiter.mass); % Mass ratio w/ primary
% bodies.europa.MR      = bodies.europa.mass / (bodies.europa.mass + bodies.jupiter.mass); %bodies.europa.mass / (bodies.europa.mass + bodies.jupiter.mass); % Mass ratio w/ primary
bodies.europa.MR      = 2.528017528540000e-05; % Value from Mar Vaquero and Javier Roa Vicens at JPL (12/2020)
bodies.europa.J2      = 4.333968885716003e-04; %(1.938209808166e-04)*sqrt(5); % from MONTE ... sqrt(5) is normalizer (vallado pg 547)
bodies.europa.mass    = bodies.europa.MR * bodies.jupiter.mass / (1 - bodies.europa.MR); % kg
% bodies.europa.mass    = 4.798778906471254e+22; warning('Switch this back')
% bodies.europa.u       = bodies.europa.mass*bodies.constants.G; % km^3 / s^2
bodies.europa.LyapTp  = 3.076475265314829;
% bodies.europa.C22     = 1.315e-4; % Bagenal, F., Dowling, T., and McKinnon, W., (eds.), Jupiter, The Planet, Satellites, and Magnetosphere, Cambridge Planetary Society, Cambridge, U.K., 2004, p. 285.
bodies.europa.C22     = bodies.europa.J2 * 0.3; % Hydrostatic equilibrium ... Synchronous Moon Theory. .... Computation of a Science orbit about Europa ..... Bursa, M., "Figure and Dynamic Parameters of Synchronously Orbiting Satellites in the Solar System," Bulletin of the Astronomical Institutes of Czechoslovakia, Vol. 40, No. 2, 1989, pp. 125–130.

%%% Ganymede
bodies.ganymede.name    = 'ganymede';
bodies.ganymede.title   = 'Ganymede';
bodies.ganymede.mass    = 1.48e23; % kg
bodies.ganymede.color   = [1,1,1]./255;
bodies.ganymede.img     = imread([mbinPath,'/textures/ganymedeSurfTex.jpg']);
bodies.ganymede.a       = 1070400; % km
bodies.ganymede.R       = 2631.2; % km
bodies.ganymede.R_n     = bodies.ganymede.R / bodies.ganymede.a;
bodies.ganymede.meanMot = 1.016444383966831e-05; %2*pi/(7.154553*86400); % rad/s
bodies.ganymede.Tp      = 2*pi/bodies.ganymede.meanMot; % sec
bodies.ganymede.MR      = 7.796293389269981e-05; %bodies.ganymede.mass / (bodies.ganymede.mass + bodies.jupiter.mass);  % Mass ratio w/ primary
bodies.ganymede.J2      = 1.254438955403253e-04; %(5.610021555811e-05)*sqrt(5); % from MONTE ... sqrt(5) is normalizer (vallado pg 547)
bodies.ganymede.LyapTp  = 3.096161527716355;

%%% Callisto
bodies.callisto.name    = 'callisto';
bodies.callisto.title   = 'Callisto';
bodies.callisto.mass    = 1.08e23; % kg
bodies.callisto.color   = [68, 36, 82]./255;
bodies.callisto.a       = 1882700; % km
bodies.callisto.R       = 2410.3; % km
bodies.callisto.R_n     = bodies.callisto.R / bodies.callisto.a;
bodies.callisto.meanMot = 4.357479662608672e-06; %2*pi/(16.689017*86400); % rad/s
bodies.callisto.MR      = 5.689306947592001e-05; %bodies.callisto.mass / (bodies.callisto.mass + bodies.jupiter.mass);  % Mass ratio w/ primary
bodies.callisto.J2      = 2.942490965436501e-05; %(1.315921964379e-05)*sqrt(5); % from MONTE ... sqrt(5) is normalizer (vallado pg 547)

% ------------------------------------------------------------------------
%%% System - Saturn
% ------------------------------------------------------------------------
%%% Saturn
bodies.saturn.name  = 'saturn';
bodies.saturn.title = 'Saturn';
bodies.saturn.mass  = 5.683e26; % kg
bodies.saturn.u     = bodies.saturn.mass*bodies.constants.G; % km^3 * s^-2
bodies.saturn.color = [1,1,1]./255;
bodies.saturn.img   = imread([mbinPath,'/textures/saturnSurfTex.jpg']);
bodies.saturn.a     = 1433.53e6; % km, https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html
bodies.saturn.R     = 58232; % km
bodies.saturn.R_ZH  = 60330.0; % km   .... reference radius used for zonal harmonics
bodies.saturn.Tp    = 2*pi*sqrt((bodies.saturn.a^3)/bodies.sun.u);
bodies.saturn.J2    = 1.629069817035e-02; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op / MONTE - Saturn solution 365
bodies.saturn.J4    = -9.343344679305e-04; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op / MONTE - Saturn solution 365
bodies.saturn.J6    = 9.003499615846e-05; % MONTE - Saturn solution 365
bodies.saturn.J8    = -1.000000000000e-05; % MONTE - Saturn solution 365

%%% Enceladus
bodies.enceladus.name    = 'enceladus';
bodies.enceladus.title   = 'Enceladus';
bodies.enceladus.color   = [136, 194, 235]./255;
bodies.enceladus.img     = imread([mbinPath,'/textures/enceladusSurfTex.jpg']);
bodies.enceladus.mass    = 1.08022e20; % kg
bodies.enceladus.u       = bodies.enceladus.mass*bodies.constants.G; % km^3/s^2
bodies.enceladus.a       = 237948; % km
bodies.enceladus.R       = 252; % km
bodies.enceladus.R_n     = bodies.enceladus.R / bodies.enceladus.a;
bodies.enceladus.meanMot = 5.307334465496030e-05; %2*pi/(1.370218 * 86400); % rad/s (from wikipedia)
bodies.enceladus.Tp      = 2*pi/bodies.enceladus.meanMot; % sec
bodies.enceladus.MR      = 1.900404354665370e-07; %bodies.enceladus.mass / (bodies.enceladus.mass + bodies.saturn.mass); % Mass ratio w/ primary
bodies.enceladus.J2      = 5.459785012180e-03; % from MONTE
bodies.enceladus.J3      = -6.672309709362e-05;% from MONTE
bodies.enceladus.LyapTp  = 3.041563160570682;

%%% Mimas
bodies.mimas.name    = 'mimas';
bodies.mimas.title   = 'Mimas';
bodies.mimas.mass    = 3.8e19; % kg
bodies.mimas.color   = [168, 168, 168]./255;
bodies.mimas.a       = 185520; % km
bodies.mimas.R       = 198.2; % km
bodies.mimas.R_n     = bodies.mimas.R / bodies.mimas.a;
bodies.mimas.meanMot = 7.716507848866653e-05; %2*pi/(0.9424218*86400); % rad/s
bodies.mimas.Tp      = 2*pi/bodies.mimas.meanMot; % sec
bodies.mimas.MR      = 6.686608738182065e-08; %bodies.mimas.mass / (bodies.mimas.mass + bodies.saturn.mass); % Mass ratio w/ primary

%%% Rhea
bodies.rhea.name    = 'rhea';
bodies.rhea.title   = 'Rhea';
bodies.rhea.mass    = 2.306e21; % kg
bodies.rhea.color   = [108,145,92]./255;
bodies.rhea.a       = 527040; % km
bodies.rhea.R       = 763.8; % km
bodies.rhea.R_n     = bodies.rhea.R / bodies.rhea.a;
bodies.rhea.meanMot = 1.609785327425133e-05; %2*pi/(4.517500*86400); % rad/s
bodies.rhea.MR      = 4.057699530080738e-06; %bodies.rhea.mass / (bodies.rhea.mass + bodies.saturn.mass);  % Mass ratio w/ primary

%%% Titan
bodies.titan.name    = 'titan';
bodies.titan.title   = 'Titan';
bodies.titan.color   = [57, 181, 134]./255;
bodies.titan.img     = imread([mbinPath,'/textures/titanSurfTex.png']);
bodies.titan.mass    = 1.345e23; % kg
bodies.titan.a       = 1221870; % km
bodies.titan.R       = 2575.5; % km
bodies.titan.R_n     = bodies.titan.R / bodies.titan.a;
bodies.titan.meanMot = 4.560806031133923e-06; %2*pi/(15.945*86400); % rad/s
bodies.titan.Tp      = 2*pi/bodies.titan.meanMot; % sec
bodies.titan.MR      = 2.366147726782945e-04; %bodies.titan.mass / (bodies.titan.mass + bodies.saturn.mass); % Mass ratio w/ primary
bodies.titan.J2      = 3.341074092664e-05; % from MONTE
bodies.titan.LyapTp  = 3.124244231852114;

% ------------------------------------------------------------------------
%%% System - Uranus
% ------------------------------------------------------------------------
%%% Uranus
bodies.uranus.name  = 'uranus';
bodies.uranus.title = 'Uranus';
bodies.uranus.img     = imread([mbinPath,'/textures/uranusSurfTex.jpg']);
bodies.uranus.mass  = 8.681e25; % kg
bodies.uranus.a     = 19.2184*bodies.constants.AU_km; % km (wiki)
bodies.uranus.u     = bodies.uranus.mass * bodies.constants.G; % km^3 * s^-2
bodies.uranus.color = [1,1,1]./255;
bodies.uranus.R     = 25362; % km
bodies.uranus.J2    = 0.00351068; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op
bodies.uranus.J4    = -.00003417; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op

%%% Cordelia
bodies.cordelia.name    = 'cordelia';
bodies.cordelia.title   = 'Cordelia';
bodies.cordelia.mass    = 4.4e16; % kg
bodies.cordelia.a       = 49751.722; % km
bodies.cordelia.R       = 20.1; % km
bodies.cordelia.R_n     = bodies.cordelia.R / bodies.cordelia.a;
bodies.cordelia.e       = 0.00026;
bodies.cordelia.i_deg   = 0.08479; % deg, to Uranus equator
bodies.cordelia.meanMot = 2.170587228950805e-04; %2*pi/ (0.335034 * 86400);
bodies.cordelia.Tp      = 2*pi/bodies.cordelia.meanMot;
bodies.cordelia.MR      = 5.068540488157864e-10; %bodies.cordelia.mass / (bodies.cordelia.mass + bodies.uranus.mass);  % Mass ratio w/ primary

%%% Ophelia
bodies.ophelia.name    = 'ophelia';
bodies.ophelia.title   = 'Ophelia';
bodies.ophelia.mass    = 5.3e16; % kg
bodies.ophelia.a       = 53763.39; % km
bodies.ophelia.R       = 21.4; % km
bodies.ophelia.R_n     = bodies.ophelia.R / bodies.ophelia.a;
bodies.ophelia.e       = 0.00992;
bodies.ophelia.i_deg   = 0.10362; % deg, to Uranus equator
bodies.ophelia.meanMot = 1.932039766654610e-04; %2*pi/ (0.37640039 * 86400);
bodies.ophelia.Tp      = 2*pi/bodies.ophelia.meanMot;
bodies.ophelia.MR      = 6.105287405557191e-10; %bodies.ophelia.mass / (bodies.ophelia.mass + bodies.uranus.mass);  % Mass ratio w/ primary

%%% Titania
bodies.titania.name    = 'titania';
bodies.titania.title   = 'Titania';
bodies.titania.mass    = 3.42e21; % kg
bodies.titania.color   = [1,1,1]./255;
bodies.titania.a       = 436300; % km
bodies.titania.R       = 788.9; % km
bodies.titania.R_n     = bodies.titania.R / bodies.titania.a;
bodies.titania.meanMot = 8.353223425815073e-06; %2*pi/(8.705867*86400); % rad/s
bodies.titania.MR      = 3.939483089135297e-05; %bodies.titania.mass / (bodies.titania.mass + bodies.uranus.mass);  % Mass ratio w/ primary

%%% Oberon
bodies.oberon.name    = 'oberon';
bodies.oberon.title   = 'Oberon';
bodies.oberon.mass    = 2.883e21; % kg
bodies.oberon.color   = [176, 113, 105]./255;
bodies.oberon.a       = 583500; % km
bodies.oberon.R       = 761.4; % km
bodies.oberon.R_n     = bodies.oberon.R / bodies.oberon.a;
bodies.oberon.meanMot = 5.401529243748597e-06; %2*pi/(13.463234*86400); % rad/s
bodies.oberon.MR      = 3.320935672646651e-05; %bodies.oberon.mass / (bodies.oberon.mass + bodies.uranus.mass);  % Mass ratio w/ primary

% ------------------------------------------------------------------------
%%% System - Neptune
% ------------------------------------------------------------------------
%%% Neptune
bodies.neptune.name  = 'neptune';
bodies.neptune.title = 'Neptune';
bodies.neptune.img     = imread([mbinPath,'/textures/neptuneSurfTex.jpg']);
bodies.neptune.mass  = 1.024e26; % kg
bodies.neptune.a     = 30.07*bodies.constants.AU_km;
bodies.neptune.u     = bodies.neptune.mass * bodies.constants.G; % km^3 * s^-2
bodies.neptune.color = [0.1961,0.3922,0.7843];
bodies.neptune.R     = 24622; % km
bodies.neptune.J2    = 0.00340843; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op
bodies.neptune.J4    = -0.0000334; % NAIF 05/18 https://ssd.jpl.nasa.gov/?gravity_fields_op

%%% Triton
bodies.triton.name    = 'triton';
bodies.triton.title   = 'Triton';
bodies.triton.mass    = 2.14e22; % kg
bodies.triton.color   = [0.1961,0.3922,0.7843];
bodies.triton.img     = imread([mbinPath,'/textures/tritonSurfTex.jpg']);
bodies.triton.a       = 354759; % km
bodies.triton.R       = 1353; % km
bodies.triton.R_n     = bodies.triton.R / bodies.triton.a;
bodies.triton.meanMot = 1.237431662696238e-05; %2*pi/(5.876854*86400); % rad/s
bodies.triton.Tp      = 2*pi/bodies.triton.meanMot;
bodies.triton.MR      = 2.089407096563804e-04; %bodies.triton.mass / (bodies.triton.mass + bodies.neptune.mass); % Mass ratio w/ primary
bodies.triton.J2      = 4.333968885716003e-04; %bodies.europa.J2 ; % Referencing "?Secular Evolution of a Satellite by Tidal Effect: Application to Triton" by A Correia that mentions using Europa J2 as estimate for Triton's
bodies.triton.LyapTp  = 3.120561393869679;

% ------------------------------------------------------------------------
%%% System - Pluto
% ------------------------------------------------------------------------
%%% Pluto
bodies.pluto.name  = 'pluto';
bodies.pluto.title = 'Pluto';
bodies.pluto.mass  = 1.309e22; % kg
bodies.pluto.u     = bodies.pluto.mass * bodies.constants.G; % km^3 * s^-2
bodies.pluto.color = [1,1,1]./255;
bodies.pluto.R     = 1187; % km

%%% Charon
bodies.charon.name    = 'charon';
bodies.charon.title   = 'Charon';
bodies.charon.mass    = 1.586e21; % kg
bodies.charon.color   = [1,1,1]./255;
bodies.charon.a       = 19591; % km
bodies.charon.R       = 606; % km
bodies.charon.R_n     = bodies.charon.R / bodies.charon.a;
bodies.charon.meanMot = 1.138559183467410e-05; %2*pi/(6.3872*86400); % rad/s
bodies.charon.MR      = 0.108067593349687; %bodies.charon.mass / (bodies.charon.mass + bodies.pluto.mass);  % Mass ratio w/ primary

end
