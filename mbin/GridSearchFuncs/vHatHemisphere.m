function [vHats] = vHatHemisphere(n_target,centerAxis)
%%% Description
%       Returns a set of vHats that somewhat evenly populate a hemisphere
%       
% --------------------------------------------------------------
%%% Inputs
%       n_target   - [1x1] desired number of vHat vectors 
%       centerAxis - Axis for planar vHats to be centered on ('x','-x','y','-y','z',or '-z')
% --------------------------------------------------------------
%%% Outputs
%       vHats - [nx3] vHat vectors in hemisphere 
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% Rough testing
% ------------------------------------
% n_target = 1297; % ~ 5 deg spacing
% n_target = 352;  % ~ 10 deg spacing
% n_target = 145;  % ~ 15 deg spacing - should be used
% n_target = 37;   % ~ 30 deg spacing - decent
% n_target = 17;   % ~ 45 deg spacing
% n_target = 5;    % ~ 90 deg spacing

% ------------------------------------
%%% Obtaining vHat Hemisphere
% ------------------------------------
%%% Get points for evenly spaced sphere
[x,y,z] = mySphere(n_target*2);

%%% Find hemisphere of points
points = find(y<=0);
x_hem  = x(points);
y_hem  = y(points);
z_hem  = z(points);

if     isequal(centerAxis,'x')
    vHats  = [-y_hem',x_hem',z_hem'];
elseif     isequal(centerAxis,'-x')
    vHats  = [y_hem',x_hem',z_hem'];
elseif isequal(centerAxis,'y')
    vHats  = [x_hem',-y_hem',z_hem'];
elseif isequal(centerAxis,'-y')
    vHats  = [x_hem',y_hem',z_hem'];
elseif isequal(centerAxis,'z')
    vHats  = [x_hem',z_hem',-y_hem'];
elseif isequal(centerAxis,'-z')
    vHats  = [x_hem',z_hem',y_hem'];
else
    warning('Input Error')
end

end
