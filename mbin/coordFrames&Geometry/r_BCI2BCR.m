function [ r_BCR ] = r_BCI2BCR( r_BCI, times, w )
%%% Inputs
%         r_BCI - [nx3] matrix of position values in Body-Centered-Inertial
%                 frame
%         times - (sec) vector of times associated with positions 
%         w     - (rad/s) angular velocity of the system 
%
%%% Outputs
%         r_BCR - [nx3] matrix of position values in Body-Centered-Rotating
%                 frame 
%
% ========================================================================
%%% Dimension checks
if size(r_BCI,2) ~= 3
    warning('Wrong dimensions for input matrix')
    return
elseif size(r_BCI,1) ~= length(times)
    warning('Dimensions of position matrix and times don''t agree')
    return
end

%%% Preallocate output vector
r_BCR = zeros(size(r_BCI));

%%% Rotate BCI into BCR
for kk = 1:length(times)
    % Create rotation angle (rad)
    theta = times(kk)*w;
    
    % Rotate position
    r_BCR(kk,:) = R3(r_BCI(kk,:),-theta);
end

end
