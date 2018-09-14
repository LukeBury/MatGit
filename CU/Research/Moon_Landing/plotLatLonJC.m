function plotLatLonJC( u, radius )
%%% Plots lat/lon map of JC values on surface of secondary
%%% Inputs
%         u - mass ratio of system
%         radius - radius of secondary body
% ========================================================================
%%% Create field of lat/lon points
lats = -90:5:90;
lons = -180:5:180;

%%% Set up mesh grid
nlats = length(lats);
nlons = length(lons);
JCs = zeros(nlats, nlons);

%%% Set body positions
rP = [-u, 0, 0];  % Primary
rS = [1-u, 0, 0]; % Secondary

%%% Determine JCs across secondary surface
for klons = 1:nlons
    for klats = 1:nlats
        %%% Turn lat/lon into SCR position
        [rSCR] = latlon2surfECEF(lats(klats), lons(klons), radius);

        %%% rSCR --> rBCR
        rBCR = rSCR + rS;
        
        %%% Calculate distance to bodies
        r1 = norm(rP - rBCR);
        r2 = norm(rBCR);
        
        %%% Calculate Jacobi constant
        JCs(klats,klons) = rBCR(1)^2 + rBCR(2)^2 + 2*(1-u)/r1 + 2*u/r2;
        
    end
end

%%% Plot JCs across secondary surface
figure
surf(lons,lats,JCs)
PlotBoi3('Longitude, °','Latitude, °','Jacobi Constant Value', 14)


end



