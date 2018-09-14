% Convert ECEF position (km) to ECI (km) given:
    % 1) ECEF position vector (km)
    % 2) Greenwich Siderial Time (rads)
function [r_ECI] = ECEF2ECI(r_ECEF, gst)
% Creating reverse R3 rotation matrix
R3 = [cos(gst) -sin(gst) 0; sin(gst) cos(gst) 0; 0 0 1];
% Rotating r_ECEF into r_ECI
r_ECI = R3 * r_ECEF; % km
end