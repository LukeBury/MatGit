function [discretizedTrajectory] = discretize_by_path_distance(trajectory, nPoints)
%%% Description
%       Used to discretize a trajectory in n states which are approximately
%       evenly spaced in path distance. This function DOES NOT interpolate,
%       but rather finds points in the input trajectory which are closest
%       to the desired discretizations. Assumes that the input trajectory
%       is sampled at a fine enough resolution so that the motion is well
%       modeled
%       
% --------------------------------------------------------------
%%% Inputs
%       trajectory - [nx6] Time history of state vector of trajectory
%       nPoints    - [1x1] Number of points to discretize the trajectory
%                    into
% --------------------------------------------------------------
%%% Outputs
%       discretizedTrajectory - [nPointsx6] Time history of input
%               trajectory discretized approximately evenly in path
%               distance
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% Calculate path distance
% ------------------------------------
%%% Number of states in input trajectory
nStates = size(trajectory,1); 

%%% Preallocate variable for storing path distance
pathDistance = zeros(nStates,1);

%%% Loop through input trajectory and calculate cumulative path distances
for kk = 2:nStates
    pathDistance(kk) = pathDistance(kk-1) + norm(trajectory(kk,1:3)-trajectory(kk-1,1:3));
end

%%% Create vector of evenly space path distances based on found range
pathDistance_discretized = linspace(0,pathDistance(end),nPoints);

%%% Preallocate output
discretizedTrajectory = zeros(nPoints,6);

%%% Loop through discretized path distance and find corresponding states
for kk = 1:nPoints
    closestIndex = abs(pathDistance - pathDistance_discretized(kk)) == min(abs(pathDistance - pathDistance_discretized(kk)));
    discretizedTrajectory(kk,:) = trajectory(closestIndex,:);
end

end