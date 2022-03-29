function [lons_new] = convert_lon180_to_lon360(lons_old)
%%% Description
% Converts a longitude system ranging from [-180, 180] to [0 360]
%       
% ------------------------------------------------------------------------
%%% Inputs
% lons_old - [nx1] column vectors of latitudes and longitudes where the
%                longitudes range from [-180, 180] (degrees)
% ------------------------------------------------------------------------
%%% Outputs
% lons_new  - [nx2] column vectors of latitudes and longitudes where the
%                 longitudes range from [0, 360] (degrees)
% ------------------------------------------------------------------------
% Created: 6/23/21
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
lons_new = zeros(size(lons_old));

for kk = 1:length(lons_old)
    if lons_old(kk) < 0
        lons_new(kk) = 360 + lons_old(kk);
    else
        lons_new(kk) = lons_old(kk);
end


end % function