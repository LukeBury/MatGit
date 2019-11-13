function [periapses] = get_periapses(positionVectors)
%%% Description
%       Finds 3D positions of periapses based on distance minima
%       
%%% Example use with CR3BP trajectories for finding periapses about the
%   secondary:
%       %%% Secondary-centered position
%       r_SCR_n = X_BCR_n(:,1:3) - [1-secondary.MR, 0, 0];
%
%       [r_SCR_periapses] = get_periapses(r_SCR_n);
% 
%       %%% Converting SCR back to BCR
%       r_BCR_periapses = r_SCR_periapses + [1-secondary.MR, 0, 0];
% --------------------------------------------------------------
%%% Inputs
%       positionVectors - [nx3] Position vectors with respect to body of
%                         interest
% --------------------------------------------------------------
%%% Outputs
%       periapses - [nx3] Position vectors of periapses with respect to
%                   body of interest
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% Loop through positions and find periapses
% ------------------------------------
periapses = [];
for kk = 2:size(positionVectors,1)-1
    lastMagnitude    = norm(positionVectors(kk-1,:));
    currentMagnitude = norm(positionVectors(kk,:));
    nextMagnitude    = norm(positionVectors(kk+1,:));
    
    if (nextMagnitude > currentMagnitude) && (lastMagnitude > currentMagnitude)
        periapses = [periapses; positionVectors(kk,:)];
    end
end

end