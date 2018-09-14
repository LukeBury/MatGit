function [ v_BCI ] = v_BCR2BCI( v_BCR, r_BCI, times, w )
%%% Inputs
%         v_BCR - [nx3] matrix of velocity values in Body-Centered-Rotating
%                 frame
%         r_BCI - [nx3] matrix of position values in Body-Centered-Inertial
%                 frame
%         times - (sec) vector of times associated with positions 
%         w     - (rad/s) angular velocity of the system 
%
%%% Outputs
%         v_BCI - [nx3] matrix of position values in Body-Centered-Inertial
%                 frame 
%
%%% Note
%         - input a negative w value to change order to BCI->BCR
% ========================================================================
%%% Dimension checks
if size(v_BCR,2) ~= 3
    warning('Wrong dimensions for input matrix')
    return
elseif size(v_BCR,1) ~= length(times)
    warning('Dimensions of position matrix and times don''t agree')
    return
elseif size(v_BCR) ~= size(r_BCI)
    warning('Dimensions of v_BCR and r_BCI don''t agree')
    return
end

%%% Preallocate output vector
v_BCI = zeros(size(v_BCR));

%%% Rotate BCR into BCI
for kk = 1:length(times)
    % Create rotation angle (rad)
    theta = times(kk)*w;
    
    % Rotate position
    v_BCI(kk,:) = R3(v_BCR(kk,:),theta) + cross([0,0,w],r_BCI(kk,:));
end

end

