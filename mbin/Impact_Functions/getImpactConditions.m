function [latLon_deg, impactAngle_deg, latLonHeadingHat, impactColor] = getImpactConditions(X_impact, prms, impactAngleBins_deg, impactAngleColors)
%%% Description
%       Returns impact location in latitude/longitude, impact angle,
%       a unit vector corresponding to direction of motion at time of
%       impact, and a color code corresponding to the impact angle bin.
%       Built for trajectories impacting the surface of the secondary body
%       
% ------------------------------------------------------------------------
%%% Inputs
%       X_impact - [1x6] State at impact of secondary surface
%
%       prms     - [struct] system parameters - fields (u - mass ratio, 
%                  R2 - normalized radius of secondary body, rNorm -
%                  normalizing constant for distances (equal to semimajor
%                  axis of system))
%
%       impactAngleBins_deg - [1xn] bin values for impact angle (bin is equal
%                             to or great than this number)
%
%       impactAngleColors   - [nx3] matrix of color codes for impact angle
%                             bins
% ------------------------------------------------------------------------
%%% Outputs
%       latLon_deg        - [1x2] Latitude and longitude of impact site
%                           (degrees)
%
%       impactAngle_deg   - [1x1] Impact angle, where 0 degrees is tangent
%                           to the surface and 90 degrees is perpendicular
%                           (degrees)
%
%       latLonHeading_hat - [1x2] Unit vector in Latitude-Longitude space
%                           pointing in direction of motion at time of
%                           impact
%
%       impactColor       - [1x3] Color code corresponding to the impact
%                           angle bin for the trajectory
% ------------------------------------------------------------------------
% Created: 10/22/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Impact Latitude/Longitude
% -------------------------------------------------
[lat_deg, lon_deg] = BCR2latlon(X_impact(1,1:3), 'secondary', prms.u);
latLon_deg = [lat_deg, lon_deg];


% -------------------------------------------------
%%% Impact Angle
% -------------------------------------------------
%%% Creating SCR position vector
rImpact_SCR_n = X_impact(1,1:3) - [1-prms.u,0,0];
rHatImpact_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);

%%% Angle between velocity and surface at impact
[impactAngle_deg] = calcImpactAngle(rImpact_SCR_n,X_impact(1,4:6),'degrees');

% -------------------------------------------------
%%% Heading at impact
% -------------------------------------------------
%%% Heading direction in state-space
headingVec = X_impact(1,4:6) - rHatImpact_SCR_n * dot(X_impact(1,4:6), rHatImpact_SCR_n);
headingHat = headingVec ./ norm(headingVec);

%%% Extend a vector 1 meter from impact site in direction of heading at
%%% impact
oneNormalizedMeter = ((1e-3)/prms.rNorm);
refPoint_SCR = rImpact_SCR_n + (headingHat .* oneNormalizedMeter);

%%% Find the latitude/longitude of the heading-extention, and difference it
%%% with the latitude/longitude of the impact site. This will give a
%%% heading at impact in latitude-longitude space instead of state space
refPoint_BCR = refPoint_SCR + [1-prms.u, 0, 0];
[latHeadingRef_deg, lonHeadingRef_deg] = BCR2latlon(refPoint_BCR, 'secondary', prms.u);
latLonHeading = [latHeadingRef_deg, lonHeadingRef_deg] - latLon_deg;
latLonHeadingHat = latLonHeading ./ norm(latLonHeading);

% -------------------------------------------------
%%% Color for impact angle by bin
% -------------------------------------------------
%%% Color for impact angle
bin_impactAngle = discretize(impactAngle_deg, impactAngleBins_deg);
impactColor = impactAngleColors(bin_impactAngle, :);

end % function