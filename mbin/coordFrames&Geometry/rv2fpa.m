function [ fpa_rad ] = rv2fpa( r, v)
%%% Inputs:
% 1) r - [nx3] Position vectors 
% 2) v - [nx3] Velocity vectors
%%% Outputs:
% 1) fpa - [nx1] (rad) flight path angle at each time
% =========================================================================
%%% Dimension Check
if size(r) ~= size(v)
    warning('Dimension mismatch with r and v')
    return
elseif size(r,2) ~=3
    warning('Bad input dimensions')
    return
end

%%% Preallocating fpa
fpa_rad = zeros(size(r,1),1);

%%% Calculating fpa
for kk = 1:size(r,1)
    rhat = r(kk,:)./norm(r(kk,:)); % unit vector of r
    v_r_scalar = dot(v(kk,:),rhat); % magnitude of radial component of v
    v_r = v_r_scalar.*rhat; % radial component of v
    fpa_rad(kk,1) = asin(norm(v_r)/norm(v(kk,:))); % rad
    
    %%% Alternative method for fpa
%     fpa(kk,1) = acos(norm(cross(r,v))/(norm(r)*norm(v))); % This is another method for fpa
end


end

