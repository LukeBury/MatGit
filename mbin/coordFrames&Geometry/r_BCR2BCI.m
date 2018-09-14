function [ r_BCI ] = r_BCR2BCI( r_BCR, times, w )
%%% Inputs
%         r_BCR - [nx3] matrix of position values in Body-Centered-Rotating
%                 frame
%         times - (sec) vector of times associated with positions 
%         w     - (rad/s) angular velocity of the system 
%
%%% Outputs
%         r_BCI - [nx3] matrix of position values in Body-Centered-Inertial
%                 frame 
%
%%% Note
%         - input a negative w value to change order to BCI->BCR
% ========================================================================
%%% Dimension checks
if size(r_BCR,2) ~= 3
    warning('Wrong dimensions for input matrix')
    return
elseif size(r_BCR,1) ~= length(times)
    warning('Dimensions of position matrix and times don''t agree')
    return
end

%%% Preallocate output vector
r_BCI = zeros(size(r_BCR));

%%% Rotate BCR into BCI
for kk = 1:length(times)
    % Create rotation angle (rad)
    theta = times(kk)*w;
    
    % Rotate position
    r_BCI(kk,:) = R3(r_BCR(kk,:),theta);
end

end

