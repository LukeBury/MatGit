clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================
run_p3 = 0;
run_p4 = 1;

% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% Bodies
% ------------------------------------
Sun     = bodies.sun;
Mercury = bodies.mercury;
Venus   = bodies.venus;
Earth   = bodies.earth;
Mars    = bodies.mars;
Jupiter = bodies.jupiter;
Saturn  = bodies.saturn;
Uranus  = bodies.uranus;
Neptune = bodies.neptune;
Pluto   = bodies.pluto;

% ------------------------------------
%%% Constants given
% ------------------------------------
au_km = 1.49597870691e8; % km

Mercury.a = 0.387  * au_km; % km
Venus.a   = 0.723  * au_km; % km
Earth.a   = 1      * au_km; % km
Mars.a    = 1.524  * au_km; % km
Jupiter.a = 5.203  * au_km; % km
Saturn.a  = 9.537  * au_km; % km
Uranus.a  = 19.191 * au_km; % km
Neptune.a = 30.069 * au_km; % km
Pluto.a   = 39.482 * au_km; % km

% ------------------------------------
%%% vInfinity values
% ------------------------------------
VInfs = [2, 4, 6, 8, 10]; % km/s

% ------------------------------------
%%% Alpha values
% ------------------------------------
alphas_rad = linspace(0,2*pi,1000); % rad

% ========================================================================
%%% P1 - Inner-Planet Tisserand Plot
% ========================================================================
% ------------------------------------
%%% Calculating circular velocities of planet
% ------------------------------------
[ V_Mercury ] = visviva_v( Mercury.a, Mercury.a, Sun.u); % km/s
[ V_Venus ]   = visviva_v( Venus.a,   Venus.a,   Sun.u); % km/s
[ V_Earth ]   = visviva_v( Earth.a,   Earth.a,   Sun.u); % km/s
[ V_Mars ]    = visviva_v( Mars.a,    Mars.a,    Sun.u); % km/s
[ V_Jupiter ] = visviva_v( Jupiter.a, Jupiter.a, Sun.u); % km/s
[ V_Saturn ]  = visviva_v( Saturn.a,  Saturn.a,  Sun.u); % km/s
[ V_Uranus ]  = visviva_v( Uranus.a,  Uranus.a,  Sun.u); % km/s
[ V_Neptune ] = visviva_v( Neptune.a, Neptune.a, Sun.u); % km/s
[ V_Pluto ]   = visviva_v( Pluto.a,   Pluto.a,   Sun.u); % km/s

% ------------------------------------
%%% Calculate Tisserand-Plot data
% ------------------------------------
% vSC_p = vP + vInf_p
rps_ras_Mercury = cell(length(VInfs),length(alphas_rad));
rps_ras_Venus   = cell(length(VInfs),length(alphas_rad));
rps_ras_Earth   = cell(length(VInfs),length(alphas_rad));
rps_ras_Mars    = cell(length(VInfs),length(alphas_rad));

rps_ras_Mercury_p2 = cell(1,length(alphas_rad));
rps_ras_Venus_p2   = cell(1,length(alphas_rad));
rps_ras_Earth_p2   = cell(1,length(alphas_rad));
rps_ras_Mars_p2    = cell(1,length(alphas_rad));
rps_ras_Jupiter_p2 = cell(1,length(alphas_rad));
rps_ras_Saturn_p2  = cell(1,length(alphas_rad));
rps_ras_Uranus_p2  = cell(1,length(alphas_rad));
rps_ras_Neptune_p2 = cell(1,length(alphas_rad));
rps_ras_Pluto_p2   = cell(1,length(alphas_rad));


for k = 1:length(VInfs)
    for j = 1:length(alphas_rad)
        %%% Find vInf_p vector at this alpha
        vInf_p = R3([VInfs(k), 0, 0], alphas_rad(j)); % km/s
        
        %%% Calculate post-flyby heliocentric velocity
        vSc_p_Mercury = vInf_p + [V_Mercury, 0, 0]; % km/s
        vSc_p_Venus   = vInf_p + [V_Venus, 0, 0];   % km/s
        vSc_p_Earth   = vInf_p + [V_Earth, 0, 0];   % km/s
        vSc_p_Mars    = vInf_p + [V_Mars, 0, 0];    % km/s
        
        %%% Find a and e of new orbit
        [aSc_Mercury, eSc_Mercury,~,~,~,~] = RV2COE([0, Mercury.a, 0], vSc_p_Mercury, Sun.u);
        [aSc_Venus,   eSc_Venus,~,~,~,~]   = RV2COE([0, Venus.a, 0],   vSc_p_Venus,   Sun.u);
        [aSc_Earth,   eSc_Earth,~,~,~,~]   = RV2COE([0, Earth.a, 0],   vSc_p_Earth,   Sun.u);
        [aSc_Mars,    eSc_Mars,~,~,~,~]    = RV2COE([0, Mars.a, 0],    vSc_p_Mars,    Sun.u);

        %%% Calculate post-flyby rp and ra, and store results
        rps_ras_Mercury{k,j} = [aSc_Mercury*(1-eSc_Mercury), aSc_Mercury*(1+eSc_Mercury)]./au_km;
        rps_ras_Venus{k,j}   = [aSc_Venus*(1-eSc_Venus),     aSc_Venus*(1+eSc_Venus)]./au_km;
        rps_ras_Earth{k,j}   = [aSc_Earth*(1-eSc_Earth),     aSc_Earth*(1+eSc_Earth)]./au_km;
        rps_ras_Mars{k,j}    = [aSc_Mars*(1-eSc_Mars),       aSc_Mars*(1+eSc_Mars)]./au_km;

    end
    
    
end


% ------------------------------------
%%% Re-structuring data for P1 plot
% ------------------------------------
rps_ras_Mercury_mat = padcat(rps_ras_Mercury{:,:});
rps_ras_Venus_mat   = padcat(rps_ras_Venus{:,:});
rps_ras_Earth_mat   = padcat(rps_ras_Earth{:,:});
rps_ras_Mars_mat    = padcat(rps_ras_Mars{:,:});

% ------------------------------------
%%% Plotting
% ------------------------------------
figure; hold all
plot(rps_ras_Mercury_mat(:,1),rps_ras_Mercury_mat(:,2),'b.')
plot(rps_ras_Venus_mat(:,1),rps_ras_Venus_mat(:,2),'r.')
plot(rps_ras_Earth_mat(:,1),rps_ras_Earth_mat(:,2),'k.')
plot(rps_ras_Mars_mat(:,1),rps_ras_Mars_mat(:,2),'g.')

pMercury = plot(Mercury.a/au_km, Mercury.a/au_km,'b^','markersize',12,'markerfacecolor','b','markeredgecolor','k');
pVenus   = plot(Venus.a/au_km, Venus.a/au_km,'r^','markersize',12,'markerfacecolor','r','markeredgecolor','k');
pEarth   = plot(Earth.a/au_km, Earth.a/au_km,'k^','markersize',12,'markerfacecolor','k','markeredgecolor','k');
pMars    = plot(Mars.a/au_km, Mars.a/au_km,'g^','markersize',12,'markerfacecolor','g','markeredgecolor','k');

pJupiter = plot([0.19, 1.6],[Jupiter.a/au_km, Jupiter.a/au_km],'--','color',colors.std.cyan,'linewidth',2);
pSaturn  = plot([0.19, 1.6],[Saturn.a/au_km, Saturn.a/au_km],'--','color',colors.std.brown,'linewidth',2);
pUranus  = plot([0.19, 1.6],[Uranus.a/au_km, Uranus.a/au_km],'--','color',colors.std.ltorange,'linewidth',2);
pNeptune = plot([0.19, 1.6],[Neptune.a/au_km, Neptune.a/au_km],'--','color',colors.std.ltblue,'linewidth',2);
pPluto   = plot([0.19, 1.6],[Pluto.a/au_km, Pluto.a/au_km],'--','color',colors.std.ltpurp,'linewidth',2);


legend([pMercury,pVenus,pEarth,pMars,pJupiter,pSaturn,pUranus,pNeptune,pPluto],...
    'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto','location','southeast')
xlim([0.19 1.6])
ylim([0.32 40])
setLogPlot('y')
PlotBoi2('Radius of Periapsis, $AU$','Radius of Apoapsis, $AU$',18,'LaTex')




% ========================================================================
%%% P2 - Route to beyond Pluto
% ========================================================================
%%% Grabbing necessary data
rpa_ras_Earth_mat_p2 = padcat(rps_ras_Earth{2,:});
rpa_ras_Venus_mat_p2 = padcat(rps_ras_Venus{4,:});
rpa_ras_Mars_mat_p2  = padcat(rps_ras_Mars{5,:});

%%% Plotting
figure; hold all
plot(rpa_ras_Earth_mat_p2(:,1),rpa_ras_Earth_mat_p2(:,2),'r.')
plot(rpa_ras_Venus_mat_p2(:,1),rpa_ras_Venus_mat_p2(:,2),'k.')
plot(rpa_ras_Mars_mat_p2(:,1),rpa_ras_Mars_mat_p2(:,2),'g.')

pMercury = plot(Mercury.a/au_km, Mercury.a/au_km,'b^','markersize',12,'markerfacecolor','b','markeredgecolor','k');
pVenus   = plot(Venus.a/au_km, Venus.a/au_km,'r^','markersize',12,'markerfacecolor','r','markeredgecolor','k');
pEarth   = plot(Earth.a/au_km, Earth.a/au_km,'k^','markersize',12,'markerfacecolor','k','markeredgecolor','k');
pMars    = plot(Mars.a/au_km, Mars.a/au_km,'g^','markersize',12,'markerfacecolor','g','markeredgecolor','k');

pJupiter = plot([0.19, 1.6],[Jupiter.a/au_km, Jupiter.a/au_km],'--','color',colors.std.cyan,'linewidth',2);
pSaturn  = plot([0.19, 1.6],[Saturn.a/au_km, Saturn.a/au_km],'--','color',colors.std.brown,'linewidth',2);
pUranus  = plot([0.19, 1.6],[Uranus.a/au_km, Uranus.a/au_km],'--','color',colors.std.ltorange,'linewidth',2);
pNeptune = plot([0.19, 1.6],[Neptune.a/au_km, Neptune.a/au_km],'--','color',colors.std.ltblue,'linewidth',2);
pPluto   = plot([0.19, 1.6],[Pluto.a/au_km, Pluto.a/au_km],'--','color',colors.std.ltpurp,'linewidth',2);

legend([pMercury,pVenus,pEarth,pMars,pJupiter,pSaturn,pUranus,pNeptune,pPluto],...
    'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto','location','southeast')
xlim([0.19 1.6])
ylim([0.32 40])
setLogPlot('y')
PlotBoi2('Radius of Periapsis, $AU$','Radius of Apoapsis, $AU$',18,'LaTex')



% ========================================================================
%%% P3
% ========================================================================

if run_p3 == 1

% ------------------------------------
%%% Calculate Tisserand-Plot data
% ------------------------------------
VInfs_p3 = linspace(8,10,1000);
minVInf_p3 = [];
for k = 1:length(VInfs_p3)
    for j = 1:length(alphas_rad)
        %%% Find vInf_p vector at this alpha
        vInf_p = R3([VInfs_p3(k), 0, 0], alphas_rad(j)); % km/s
        
        %%% Calculate post-flyby heliocentric velocity
        vSc_p_Mars    = vInf_p + [V_Mars, 0, 0];    % km/s


        
        %%% Find a and e of new orbit
        [aSc_Mars,    eSc_Mars,~,~,~,~]    = RV2COE([0, Mars.a, 0],    vSc_p_Mars,    Sun.u);

        %%% Calculate post-flyby rp and ra, and store results
        ra_postMars    = aSc_Mars*(1+eSc_Mars);
        
        if ra_postMars > Pluto.a
            minVInf_p3 = VInfs_p3(k)
            break
        end
        
    end
    
    if isempty(minVInf_p3) == 0
        break
    end
    
end


end % if run_p3 == 1





% ========================================================================
%%% P4
% ========================================================================

if run_p4 == 1
VInfs_p4 = linspace(13,15,1000);
minVInf_p4 = [];
for k = 1:length(VInfs_p4)
    for j = 1:length(alphas_rad)
        %%% Find vInf_p vector at this alpha
        vInf_p = R3([VInfs_p4(k), 0, 0], alphas_rad(j)); % km/s
        
        %%% Calculate post-flyby heliocentric velocity
        vSc_p_Venus    = vInf_p + [V_Venus, 0, 0];    % km/s
        
        %%% Find a and e of new orbit
        [aSc_Venus,    eSc_Venus,~,~,~,~]    = RV2COE([0, Venus.a, 0],    vSc_p_Venus,    Sun.u);

        %%% Calculate post-flyby rp and ra, and store results
        ra_postVenus    = aSc_Venus*(1+eSc_Venus);
        
        if ra_postVenus > Neptune.a
            minVInf_p4 = VInfs_p4(k)
            break
        end
        
    end
    
    if isempty(minVInf_p4) == 0
        break
    end
    
end





end % run_p4 == 1

