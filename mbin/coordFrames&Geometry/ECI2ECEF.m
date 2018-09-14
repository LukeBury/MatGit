% Convert ECI (km) position to ECEF (km) given:
    % 1) ECI position vector (km)
    % 2) Greenwich Siderial Time (rads)
function [r_ECEF] = ECI2ECEF(r_ECI, gst)
% Creating R3 rotation matrix
R3 = [cos(gst) sin(gst) 0; -sin(gst) cos(gst) 0; 0 0 1];
% Rotating r_ECI into r_ECEF
r_ECEF = R3 * r_ECI; % km
end
