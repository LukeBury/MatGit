clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))

% ========================================================================
%%% Importing Data
% ========================================================================
% %%% General data on solar system bodies
bodies = getBodyData(mbinPath);
% 
% %%% Color options/schemes
colors = get_colors();



secondary = bodies.europa;
primary   = bodies.jupiter;
prms.u = secondary.MR;
prms.n = 1;

L123 = EquilibriumPoints(secondary.MR, 1,1:3);

rNorm = secondary.a;         % n <-> km
tNorm_old = 1/secondary.meanMot; % n <-> sec
vNorm_old = rNorm / tNorm_old;       % n <-> km/sec

[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

JC_50 = L2FlyoverVelocity_2_JC(300, secondary.MR, L123(2,:), vNorm_old, 1)

l2_mps     = JC_2_L2FlyoverVelocity(JC_50, prms, L123(2,:), vNorm)
l2_mps_old = JC_2_L2FlyoverVelocity(JC_50, prms, L123(2,:), vNorm_old)
% ========================================================================
%%% Testing something 
% ========================================================================
% r0 = [0.003; 0.002; 0.002];
% v0 = [-0.04; 0.05; -0.05];
% 
% X0 = [r0; v0];
% 
% 
% tol = 1e-13;
% options = odeset('RelTol',tol,'AbsTol',tol);
% 
% 
% [~, X] = ode113(@Int_2BI, [0 0.5*pi], X0, options, bodies.europa.MR);
% 
% figure; hold all
% PlotBoi3(' ',' ',' ',13)
% plot3(X(:,1),X(:,2),X(:,3),'linewidth',4,'color',colors.cyan)
% plot3(X(1,1),X(1,2),X(1,3),'bo')
% plotBodyTexture3(bodies.europa.R_n, [0,0,0], bodies.europa.img)
% axis equal
% grid off
% axis off
% set(gcf,'color',colors.orccaPPT)
% % ========================================================================
% %%% Looking at enceladus polar landings
% % ========================================================================
% 
% dataPath       = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/ManifoldShallowImpactDataSets/shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_All.nodes2000.txt';
% enceladusTrajs = dlmread(dataPath,',',1,0);
% polarLandings  = enceladusTrajs(abs(enceladusTrajs(:,8))>80, :);
% 
% 
% primary = bodies.saturn;
% secondary = bodies.enceladus;

% tol = 1e-13;
% options = odeset('RelTol',tol,'AbsTol',tol);
% 
% prms.u = secondary.MR;
% 
% 
% 
% 
% 
% %%% 6/1/20 testing my C_J2 against the C_J2 presented in a paper I'm
% %%% reviewing
% xd = rand(1);
% yd = rand(1);
% zd = rand(1);
% x = rand(1);
% y = rand(1);
% z = rand(1);
% u = rand(1);
% r1 = sqrt((u+x)^2 + y^2 + z^2);
% r2 = sqrt((x-(1-u))^2 + y^2 + z^2);
% R1 = rand(1);
% R2 = rand(1);
% J2p = rand(1);
% J2s = rand(1);
% n = 1;
% 
% Cj2_me = -(xd^2 + yd^2 + zd^2) + x^2 + y^2 + 2*(1-u)/r1 + 2*u/r2 - (1-u)*(R1^2)*J2p*(3*z*z-r1^2)/(r1^5) - u*(R2^2)*J2s*(3*z*z-r2^2)/(r2^5)
% 
% Cj2_them = -(xd^2 + yd^2 + zd^2) + 2*(1-u)*(1/r1 + (3*J2p*R1*R1/(2*(r1^3)))*(1/3 - z*z/(r1^2))) + 2*u*(1/r2 + (3*J2s*R2*R2/(2*(r2^3)))*(1/3 - z*z/(r2^2))) + (n^2)*((1-u)*r1*r1 + u*r2*r2 - z*z)




% % n_trajs = size(polarLandings,1);
% n_trajs = size(polarLandings,1)-2000;
% 
% trajs = [];
% for kk = 1:n_trajs
%     [T, X] = ode113(@Int_CR3Bn, [0, polarLandings(kk, 7)], polarLandings(kk,1:6)', options, prms);
%     
%     trajs = [trajs; X; NaN(1,6)];
%     
%     if mod(kk,500) == 0
%         fprintf('%1.1f\n',kk*100/n_trajs)
%     end
% end
% figure; hold all
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(20)
% axis equal
% plot3(trajs(:,1),trajs(:,2),trajs(:,3),'r')
% 
% figure; hold all
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(20)
% axis equal
% plot3(X(:,1),X(:,2),X(:,3),'r')






% %%% Filename
% fileName = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/ManifoldShallowImpactDataSets/PolarShallowImpacts_SaturnEnceladus.txt';
% 
% %%% Open File
% datafile = fopen(fileName,'wt');
% 
% %%% Write header
% headerString = 'x0,y0,z0,xd0,yd0,zd0,Tf,impact_lat,impact_lon,JC\n';
% fprintf(datafile,headerString);
% 
% %%% Write data
% dataStart = 1;
% dataStop = n_trajs;
% 
% for kk = dataStart:dataStop
%     fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.2f,%1.2f,%1.16f\n',...
%         polarLandings(kk,1), polarLandings(kk,2), polarLandings(kk,3),...
%         polarLandings(kk,4), polarLandings(kk,5), polarLandings(kk,6),...
%         polarLandings(kk,7), polarLandings(kk,8), polarLandings(kk,9),...
%         polarLandings(kk,13));
% 
% end
% 
% %%% Close file
% fclose(datafile);
% ========================================================================
%%% Checking dot(r,v) for 2 body orbit
% ========================================================================
% earth = bodies.earth;
% tol = 1e-13;
% options = odeset('RelTol',tol,'AbsTol',tol);
% 
% tf = 3600*8;
% t_vec = linspace(0,tf,1000);
% 
% X0 = [12000;0;0;0;7;0];
% 
% [T, X] = ode113(@Int_2BI, t_vec, X0, options, earth.u);
% 
% figure; hold all
% plot(X(:,1),X(:,2),'r')
% plotBody2(earth.R,[0,0,0],colors.blue,colors.black,1,1)
% PlotBoi2('x','y',20,'LaTex')
% axis equal
% 
% 
% RVs = NaN(size(X,1),1);
% for kk = 1:length(RVs)
%     RVs(kk) = dot(X(kk,1:3),X(kk,4:6));
% end
% 
% figure; hold all
% plot(T, RVs)
% PlotBoi2('t','r dot v',20,'LaTex')
% ========================================================================
%%% Checking stability indicies for Jake Vendl
% ========================================================================
% trajFile = load('/Users/lukebury/Downloads/Distant_Prograde_ICs.mat');
% trajData = trajFile.Distant_Prograde_ICs;
% 
% n_trajs = size(trajData,1);
% 
% primary   = bodies.earth;
% secondary = bodies.moon;
% 
% stm0_vec = reshape(eye(6),36,1);
% 
% tol = 1e-13;
% options = odeset('RelTol',tol,'AbsTol',tol);
% 
% prms.u    = secondary.MR;
% prms.R2_n = secondary.R_n;
% 
% stabilityIndices = NaN(n_trajs,2);
% JCs              = NaN(n_trajs,1);
% for kk = 1:n_trajs
%     X0n_i = trajData(kk,1:6)';
%     Tpn_i = trajData(kk,7);
%     
%     [~, Xn_i] = ode113(@Int_CR3BnSTM, [0 Tpn_i], [X0n_i; stm0_vec], options, prms);
%     stm_tf_t0                           = reshape(Xn_i(end,7:42),6,6);
%     monodromy                           = stm_tf_t0;
%     [eigenVectors_new, eigenValues_new] = eig(monodromy);
%     [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
%     stabilityIndices(kk,:)              = [S1, S2];
%     
%     JCs(kk) = getJacobiConstant_ZH(X0n_i', prms);
% end
% 
% figure;  hold all
% s1 = plot(JCs, stabilityIndices(:,1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
% s2 = plot(JCs, stabilityIndices(:,2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
% PlotBoi2('Jacobi Constant','Stability Indices',18,'LaTex')
% plot(unique([min(JCs) max(JCs)]),[2 2],'k','linewidth',1)
% xlim(unique([min(JCs) max(JCs)]))
% legend([s1, s2],'S_1','S_2')
% 
% for kk = 1:n_trajs
%     fprintf('%1.9f,%1.9f,%1.9f\n',stabilityIndices(kk,1),stabilityIndices(kk,2),JCs(kk))
% end

% ========================================================================
%%% 02/13/20 testing 
% ========================================================================
% 
% % trajs_Low = dlmread("/Users/lukebury/CU_Google_Drive/Documents/JulGit/Research/HighLatitudeLanding/equatorialImpactTrajs_noImpact6pi.txt",',',1,0);
% % trajs_High = dlmread("/Users/lukebury/CU_Google_Drive/Documents/JulGit/Research/HighLatitudeLanding/highLatImpactTrajs.txt",',',1,0);
% 
% trajs_Low = dlmread("/Users/lukebury/CU_Google_Drive/Documents/JulGit/Research/HighLatitudeLanding/equatorialImpactTrajs_150mps.txt",',',1,0);
% trajs_High = dlmread("/Users/lukebury/CU_Google_Drive/Documents/JulGit/Research/HighLatitudeLanding/highLatImpactTrajs_150mps_70lat.txt",',',1,0);
% 
% figure; hold all
% p_low  = plot3(trajs_Low(:,1),trajs_Low(:,2),trajs_Low(:,3),'b','linewidth',2);
% p_high = plot3(trajs_High(:,1),trajs_High(:,2),trajs_High(:,3),'r','linewidth',2);
% 
% PlotBoi3_CR3Bn(20)
% axis equal
% 
% prms.u = bodies.europa.MR;
% prms.R2_n = bodies.europa.R / bodies.europa.a;
% JC = getJacobiConstant_ZH(trajs_Low(2,1:6),prms);
% plotCR3BP_YZNeck(JC, prms.u, 2, 0, prms, colors.black, 2)
% plotCR3BP_YZNeck(JC, prms.u, 1, 0, prms, colors.black, 2)
% 
% L123 = EquilibriumPoints(prms.u,[1,2,3]);
% plotCR3BP_Neck(prms, L123, JC, 600, 200, colors.black, 2)
% 
% legend('"Low" Trajectories','"High" Trajectories', 'Zero-Velocity Curves')
% 
% % figure; hold all
% % plot3(trajs_High(:,1),trajs_High(:,2),trajs_High(:,3),'r')
% % PlotBoi3('x','y','z',20)
% % axis equal

% ========================================================================
%%% 02/11/20 Testing (rotation matrix between two vectors
% ========================================================================
% % a = [1; 3; 7]; b = [-2; 5; -8];
% % 
% % a = a ./ norm(a);
% % b = b ./ norm(b); 
% % 
% % v = cross(a,b);
% % 
% % ss = [0, -v(3), v(2);...
% %       v(3), 0, -v(1);...
% %       -v(2), v(1), 0];
% % 
% % R = eye(3) + ss + ss*ss*(1 / (1 + dot(a,b)));
% % 
% % R * b
% n = 1000;
% hem = vHatHemisphere(n, 'x');
% 
% % a = [1; 0; 0];
% a = [1; 0; 0];
% b = [-1; 0; 0.01];
% 
% a = a ./ norm(a);
% b = b ./ norm(b); 
% 
% v = cross(a,b);
% ss = [0, -v(3), v(2);...
%       v(3), 0, -v(1);...
%       -v(2), v(1), 0];
% R = eye(3) + ss + ss*ss*(1 / (1 + dot(a,b)));
% 
% newHem = zeros(size(hem));
% 
% for kk = 1:n
%     newHem(kk,:) = (R * hem(kk,:)')';
% end
% 
% figure; hold all
% plot3(hem(:,1),hem(:,2),hem(:,3),'.')
% quiver3(0,0,0,a(1),a(2),a(3),'b')
% quiver3(0,0,0,b(1),b(2),b(3),'r')
% axis equal
% PlotBoi3('x','y','z',10)
% 
% figure; hold all
% plot3(newHem(:,1),newHem(:,2),newHem(:,3),'.')
% quiver3(0,0,0,a(1),a(2),a(3),'b')
% quiver3(0,0,0,b(1),b(2),b(3),'r')
% axis equal
% PlotBoi3('x','y','z',10)
% ========================================================================
%%% 11/14/19 Testing
% % ========================================================================
% tocWhole = 12;
% % receiver = 'broncosferever@gmail.com';
% receiver = '6052547400@vtext.com';
% 
% subject = sprintf('Fortuna run complete!');
% message = sprintf('Ellapsed time: %1.0f seconds\n\n Great job, buddy',tocWhole);
%     
% sendEmailFromMatlab(receiver, subject, message)


% % ========================================================================
% %%% 05/28/19 Testing
% % ========================================================================
% %%% Systems
% 
% primary1 = bodies.jupiter;
% primary2 = bodies.saturn;
% 
% secondary1 = bodies.europa;
% secondary2 = bodies.ganymede;
% secondary3 = bodies.enceladus;
% secondary4 = bodies.titan;
% secondary5 = bodies.moon;
% 
% %%% Normalizing constants
% rNorm1 = secondary1.a;         % n <-> km
% tNorm1 = 1/secondary1.meanMot; % n <-> sec
% vNorm1 = rNorm1 / tNorm1;       % n <-> km/sec
% 
% rNorm2 = secondary2.a;         % n <-> km
% tNorm2 = 1/secondary2.meanMot; % n <-> sec
% vNorm2 = rNorm2 / tNorm2;       % n <-> km/sec
% 
% rNorm3 = secondary3.a;         % n <-> km
% tNorm3 = 1/secondary3.meanMot; % n <-> sec
% vNorm3 = rNorm3 / tNorm3;       % n <-> km/sec
% 
% rNorm4 = secondary4.a;         % n <-> km
% tNorm4 = 1/secondary4.meanMot; % n <-> sec
% vNorm4 = rNorm4 / tNorm4;       % n <-> km/sec
% 
% rNorm5 = secondary5.a;         % n <-> km
% tNorm5 = 1/secondary5.meanMot; % n <-> sec
% vNorm5 = rNorm5 / tNorm5;       % n <-> km/sec
% 
% 
% tBurn_hrs = 12; % hours
% tBurn_sec = tBurn_hrs * 3600; % sec
% tBurn_n1 = tBurn_hrs / tNorm1
% tBurn_n2 = tBurn_hrs / tNorm2
% tBurn_n3 = tBurn_hrs / tNorm3
% tBurn_n4 = tBurn_hrs / tNorm4
% tBurn_n5 = tBurn_hrs / tNorm5
% 
% primary = bodies.jupiter;
% secondary = bodies.europa;
% 
% L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
% 
% %%% Normalizing constants
% rNorm = secondary.a;         % n <-> km
% tNorm = 1/secondary.meanMot; % n <-> sec
% vNorm = rNorm / tNorm;       % n <-> km/sec
% 
% 
% 
% %%% Choosing ode tolerance
% tol = 1e-13;
% 
% %%% Setting integrator options
% options = odeset('Event',@event_ImpactEscape_CR3Bn, 'RelTol',tol,'AbsTol',tol);
% 
% %%% Setting necessary parameters for integration
% prms.u = secondary.MR;
% prms.R2_n = secondary.R_n;
% prms.L1x = L123(1,1);
% prms.L2x = L123(2,1);
% 
% % X0 = [1.0204617015266166;-0.0012084061620286;0.0000000000000000;-0.0028596241729016;0.0009291482175995;0.0000000000000000]';
% % X0 = [1.0204617015266166,-0.0019850725823437,0.0000000000000000,-0.0010976413130697,-0.0030157447223195,0.0098771743854318];
% % X0 = [1.020461701526617; 0.000037454199667; 0.000074923811876; -0.010853703861100; 0.001140770244121; 0]';
% X0 = [1.020461701526617; 0.000410765286112; 0.003440050933435; -0.003449156947614; -0.001991371692182; 0.012257623746852];
% t_f = 12.49307;
% 
% t_i = 0; % sec
% % t_f = 4*pi; % Long bc events are watching for impact or escape
% n_dt = 10000;
% time0_n = linspace(t_i,t_f,n_dt);
% 
% stm0 = reshape(eye(6,6),36,1);
% X0 = [X0; stm0];
% 
% [time_n, X_BCR_n] = ode113(@Int_CR3BnSTM, time0_n, X0, options, prms);
% 
% r_SCR_n = X_BCR_n(:,1:3) - [1-secondary.MR, 0, 0];
% R_SCR_n = rowNorm(r_SCR_n);
% 
% figure(1); hold all
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'r','linewidth',1.5)
% [JC_scInitial] = JacobiConstantCalculator(secondary.MR,X0(1:3)',X0(4:6)');
% plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 2, 0, prms, colors.std.black, 1.5)
% plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 1, 0, prms, colors.std.black, 1.5)
% plotCR3BP_Neck(secondary,L123,JC_scInitial,600,200,colors.std.black,1.5)
% % plotBody2(secondary.R_n,[1-secondary.MR,0,0],colors.std.blue,colors.std.black,1,0)
% plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
% view(0,90)
% axis equal
% PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
% % importantIndices = [4991, 6737, 8271, 9081]; % spikes in STM determinant
% importantIndices = [1235, 3066, 5004, 6733, 8284]; % sharp peaks in z-dot
% 
% plot3(X_BCR_n(importantIndices,1),X_BCR_n(importantIndices,2),X_BCR_n(importantIndices,3),'m.','markersize',15)
% % plot3(X_BCR_n(importantIndices(3),1),X_BCR_n(importantIndices(3),2),X_BCR_n(importantIndices(3),3),'k.','markersize',16)
% 
% determinants = NaN(size(X_BCR_n,1),1);
% determinants_pos = NaN(size(X_BCR_n,1),1);
% determinants_vel = NaN(size(X_BCR_n,1),1);
% for kk = 1:size(X_BCR_n,1)
%     stm = reshape(X_BCR_n(kk,7:42),6,6);
%     stm_pos = stm(1:3,1:3);
%     stm_vel = stm(4:6,4:6);
%     
%     determinants(kk)     = det(stm);
%     determinants_pos(kk) = det(stm_pos);
%     determinants_vel(kk) = det(stm_vel);
% end
% figure
% plot(percentchange(determinants),'linewidth',2)
% title('|\Phi|')
% PlotBoi2('Index','\% Change in $|\Phi|$',18,'LaTex')
% figure
% plot(percentchange(determinants_pos),'linewidth',2)
% title('|\Phi_{pos}|')
% PlotBoi2('Index','\% Change in $|\Phi|$',18,'LaTex')
% figure
% plot(percentchange(determinants_vel),'linewidth',2)
% title('|\Phi_{vel}|')
% PlotBoi2('Index','\% Change in $|\Phi|$',18,'LaTex')
% 
% figure
% plot(R_SCR_n,X_BCR_n(:,3),'linewidth',2)
% PlotBoi2('R_{SCR}','Z Amplitude',18,'LaTex')
% 
% 
% figure
% plot(X_BCR_n(:,3),'linewidth',2)
% PlotBoi2('Index','Z Amplitude',18,'LaTex')
% 
% figure
% plot(X_BCR_n(:,6),'linewidth',2)
% PlotBoi2('Index','Z-dot Amplitude',18,'LaTex')

% ========================================================================
%%% Testing
% ========================================================================

% %%% Choosing ode tolerance
% tol = 1e-11;
% %%% Setting integrator options
% options = odeset('RelTol',tol,'AbsTol',tol);
% t_i = 0; % sec
% t_f = 2*365.242189*86400; % Long bc events are watching for impact or escape
% n_dt = 10000;
% time0 = linspace(t_i,t_f,n_dt);
% 
% X0 = [-1.116247042938662e+08;...
%     9.651787718164873e+07;...
%     -5.115716280996529e+03;...
%     3.056174806192775e+05;...
%     -2.642564712324188e+05;...
%     14.006328899915344];
% [~, X_out] = ode113(@Int_2BI, time0, X0, options, Sun.u);




% JD_departure = 2460545;
% JD_arrival   = 2460919;
% transferTime_sec = (JD_arrival - JD_departure)*86400;
% 
% [L_rad, a, e, i_rad, Omega_rad, Pi_rad, w_rad, M_rad, ta_rad] = ...
%     getPlanetElements_Meeus(JD_departure, 'Earth', 'radians');
% [r_departure, v_departure] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
% 
% %%% Find Jupiter position at arrival time
% [L_rad, a, e, i_rad, Omega_rad, Pi_rad, w_rad, M_rad, ta_rad] = ...
%     getPlanetElements_Meeus(JD_arrival, 'Venus', 'radians');
% [r_arrival, v_arrival] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
% 
% % r_departure
% % r_arrival
% 
% nOrbits = 1;
% [V1, V2,exitflag] = lambertSolver(r_departure, r_arrival, transferTime_sec, nOrbits, 0, bodies.sun.u);
% 
% transferType = 3;
% [v1,v2] = lambertTargeting(r_departure, r_arrival, transferTime_sec, transferType, 1, bodies.sun.u, 0);
% 
% %%% Creating time vector
% t0 = 0;             % sec
% tf = transferTime_sec; % sec
% time = linspace(t0,tf,10000);
% 
% 
% %%% Choosing ode tolerance
% tol = 2.22045e-14;
% 
% %%% Setting integrator options
% options = odeset('RelTol',tol,'AbsTol',tol);
% 
% 
% [time, X_Earth] = ode113(@Int_2BI,time, [r_departure; v_departure], options, bodies.sun.u);
% [time, X_Venus] = ode113(@Int_2BI,time, [r_arrival; v_arrival], options, bodies.sun.u);
% 
% [time_sc, X_sc] = ode113(@Int_2BI,time, [r_departure; V1], options, bodies.sun.u);
% 
% [time_sc, X_sc_m] = ode113(@Int_2BI,time, [r_departure; v1], options, bodies.sun.u);
% 
% figure; hold all
% p1 = plot3(X_Earth(:,1),X_Earth(:,2),X_Earth(:,3),'color',colors.std.blue,'linewidth',2);
% plot3(X_Earth(1,1),X_Earth(1,2),X_Earth(1,3),'ro','linewidth',2)
% p2 = plot3(X_Venus(:,1),X_Venus(:,2),X_Venus(:,3),'color',colors.std.orange,'linewidth',2);
% plot3(X_Venus(1,1),X_Venus(1,2),X_Venus(1,3),'go','linewidth',2)
% p3 = plot3(X_sc(:,1),X_sc(:,2),X_sc(:,3),'--k','linewidth',2);
% p4 = plot3(X_sc_m(:,1),X_sc_m(:,2),X_sc_m(:,3),'--c','linewidth',2);
% PlotBoi3('x','y','z',15)
% axis equal
% legend([p1 p2 p3 p4],{'Earth','Venus','Old Lambert','My Lambert'})



% primary = bodies.jupiter;
% secondary = bodies.europa;
% 
% L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
% 
% %%% Normalizing constants
% rNorm = secondary.a;         % n <-> km
% tNorm = 1/secondary.meanMot; % n <-> sec
% vNorm = rNorm / tNorm;       % n <-> km/sec
% 
% 
% 
% %%% Choosing ode tolerance
% tol = 1e-13;
% 
% %%% Setting integrator options
% options = odeset('Event',@event_ImpactEscape_CR3Bn, 'RelTol',tol,'AbsTol',tol);
% 
% %%% Setting necessary parameters for integration
% prms.u = secondary.MR;
% prms.R2_n = secondary.R_n;
% prms.L1x = L123(1,1);
% prms.L2x = L123(2,1);
% 
% % X0 = [1.0204617015266166;-0.0012084061620286;0.0000000000000000;-0.0028596241729016;0.0009291482175995;0.0000000000000000]';
% % X0 = [1.0204617015266166,-0.0019850725823437,0.0000000000000000,-0.0010976413130697,-0.0030157447223195,0.0098771743854318];
% % X0 = [1.020461701526617; 0.000037454199667; 0.000074923811876; -0.010853703861100; 0.001140770244121; 0]';
% X0 = [1.020461701526617; 0.000410765286112; 0.003440050933435; -0.003449156947614; -0.001991371692182; 0.012257623746852];
% t_f = 12.49307;
% 
% t_i = 0; % sec
% % t_f = 4*pi; % Long bc events are watching for impact or escape
% n_dt = 10000;
% time0_n = linspace(t_i,t_f,n_dt);
% 
% [time_n, X_BCR_n] = ode113(@Int_CR3Bn, time0_n, X0, options, prms);
% 
% figure; hold all
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'r','linewidth',1.5)
% [JC_scInitial] = JacobiConstantCalculator(secondary.MR,X0(1:3)',X0(4:6)');
% plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 2, 0, prms, colors.std.black, 1.5)
% plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 1, 0, prms, colors.std.black, 1.5)
% plotCR3BP_Neck(secondary,L123,JC_scInitial,600,200,colors.std.black,1.5)
% % plotBody2(secondary.R_n,[1-secondary.MR,0,0],colors.std.blue,colors.std.black,1,0)
% plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
% view(0,90)
% axis equal
% PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')






% f_me = @(x) x - (1-u)*(x+u)/(abs(x+u)^3) - u*(x-1+u)/(abs(x-1+u)^3) - 3*(1-u)*(prms.R1^2)*prms.J2p*(x+u)/(2*(abs(x+u)^5)) - 3*u*(prms.R2^2)*prms.J2s*(x-1+u)/(2*(abs(x-1+u)^5))
% 
% 
% f = @(x) x*(1 + 1.5*(A1 + A2)) - (1-u)*(x-u)*(1/(abs(x-u)^3) + 3*A1/(2*(abs(x-u)^5))) - u*(x-u+1)*(1/(abs(x-u+1)^3) + 3*A2/(2*(abs(x-u+1)^5)))
% f_test = @(x) x - (1-u)*(x-u)*(1/(abs(x-u)^3) + 3*A1/(2*(abs(x-u)^5))) - u*(x-u+1)*(1/(abs(x-u+1)^3) + 3*A2/(2*(abs(x-u+1)^5)))









% 
% % s_max      = 1;
% % trajID_max = 1;
% % for kk = 1:990
% %     s = size(data_L2_4pi(data_L2_4pi(:,1)==kk, c_x:c_zd),1);
% %     if s > s_max
% %         s_max = s;
% %         trajID_max = kk;
% %     end
% % end
% 
% trajID_low = 209;
% testTraj_low = data_L2_4pi(data_L2_4pi(:,1)==trajID_low, [c_x:c_zd,c_t]);
% figure; hold all
% plot3(testTraj_low(:,1),testTraj_low(:,2),testTraj_low(:,3),'b')
% plotSecondary(secondary)
% axis equal
% PlotBoi3_CR3Bn(20)
% 
% % trajID_high = [1];
% % testTraj_high =  data_high_1pi(data_high_1pi(:,1)==trajID_high, c_x:c_zd);
% figure; hold all
% plotSecondary(secondary)
% axis equal
% PlotBoi3_CR3Bn(20)
% 
% for kk = 1:250:50879
%     testTraj_high =  data_high_1pi(data_high_1pi(:,1)==kk, c_x:c_zd);
%     plot3(testTraj_high(:,1),testTraj_high(:,2),testTraj_high(:,3),'r')
% end
% 
% 
% 
% savepath_test = '/Volumes/LB_External_Drive/Research/High_Latitude_Landing_Study/Data/testData';
% 
% 
% 
% filename_testTraj_low = fullfile(savepath_test, sprintf('testData_lowLatTraj.txt'));
% f_testTraj_low = fopen(filename_testTraj_low, 'wt');
% fprintf(f_testTraj_low,'trajID,x_n,y_n,z_n,xd_n,yd_n,zd_n,t_n\n');  % header
% for kk = 1:size(testTraj_low,1)
%     fprintf(f_testTraj_low,'%1d,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f\n',209, testTraj_low(kk,1), testTraj_low(kk,2), testTraj_low(kk,3), testTraj_low(kk,4), testTraj_low(kk,5), testTraj_low(kk,6), testTraj_low(kk,7));
% end
% fclose(f_testTraj_low);
% 
% 
% 
% filename_testTraj_high = fullfile(savepath_test, sprintf('testData_highLatTraj.txt'));
% f_testTraj_high = fopen(filename_testTraj_high, 'wt');
% fprintf(f_testTraj_high,'trajID,x_n,y_n,z_n,xd_n,yd_n,zd_n,t_n\n');  % header
% for trajID_kk = 1:250:50879
%     testTraj_high =  data_high_1pi(data_high_1pi(:,1)==trajID_kk, [c_x:c_zd,c_t]);
%     
%     for kk = 1:size(testTraj_high,1)
%         fprintf(f_testTraj_high,'%1d,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f,%1.13f\n',trajID_kk, testTraj_high(kk,1), testTraj_high(kk,2), testTraj_high(kk,3), testTraj_high(kk,4), testTraj_high(kk,5), testTraj_high(kk,6), testTraj_high(kk,7));
%     end
% end
% fclose(f_testTraj_high);




