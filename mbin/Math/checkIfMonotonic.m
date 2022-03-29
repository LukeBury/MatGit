function [isMonotonic, turnPoints] = checkIfMonotonic(input_vector)
%%% Description
% Check if a vector is monotonic
%       
% ------------------------------------------------------------------------
%%% Inputs
% input_vector - [nx1 or 1xn] Vector to be checked
% ------------------------------------------------------------------------
%%% Outputs
% isMonotonic - [logical] 'true' or 'false' for whether or not vector is
%                monotonic
% turnPoints  - [1xn] indices where monotonicity is broken
% ------------------------------------------------------------------------
% Created: 01/12/21
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------
% % % % % % % % %%% Establish vector for turnpoints
% % % % % % % % turnPoints = [];
% % % % % % % % 
% % % % % % % % %%% Calculate signs of differences between elements
% % % % % % % % signs = sign(diff(input_vector));
% % % % % % % % 
% % % % % % % % %%% Set the first non-zero sign (in case first sign is zero)
% % % % % % % % signed_indices = find(signs ~= 0);
% % % % % % % % firstNonZeroSign = signs(signed_indices(1));
% % % % % % % % 
% % % % % % % % %%% If the length of unique signs is 1, it's a monotonic vector
% % % % % % % % if length(unique(signs)) == 1
% % % % % % % %     isMonotonic = true;
% % % % % % % % 
% % % % % % % % %%% If the length of unique signs is 2, it's only monotonic if one of the
% % % % % % % % %%% unique signs is 0
% % % % % % % % elseif length(unique(signs)) == 2
% % % % % % % %     if ismember(0, signs)
% % % % % % % %         isMonotonic = true;
% % % % % % % %     else
% % % % % % % %         isMonotonic = false;
% % % % % % % %         if firstNonZeroSign == 1
% % % % % % % %             turnPoints = find(signs == -1) + 1;
% % % % % % % %         elseif firstNonZeroSign == -1
% % % % % % % %             turnPoints = find(signs == 1) + 1;
% % % % % % % %         else
% % % % % % % %             warning('Not sure what happened here1')
% % % % % % % %         end
% % % % % % % %         
% % % % % % % %         turnPoints2 = [];
% % % % % % % %         for kk = 1:length(turnPoints)
% % % % % % % %             if ~isequal(signs(turnPoints(kk)-1), signs(turnPoints(kk)-2))
% % % % % % % %                 turnPoints2 = [turnPoints2; turnPoints(kk)];
% % % % % % % %             end
% % % % % % % %         end
% % % % % % % %     end
% % % % % % % % 
% % % % % % % % %%% If the length of unique signs is 3, it's not a monotonic vector
% % % % % % % % elseif length(unique(signs)) == 3
% % % % % % % %     isMonotonic = false;
% % % % % % % %     if firstNonZeroSign == 1
% % % % % % % %         turnPoints = find(signs == -1) + 1;
% % % % % % % %     elseif firstNonZeroSign == -1
% % % % % % % %         turnPoints = find(signs == 1) + 1;
% % % % % % % %     else
% % % % % % % %         warning('Not sure what happened here2')
% % % % % % % %     end
% % % % % % % %     
% % % % % % % %     turnPoints2 = [];
% % % % % % % %     for kk = 1:length(turnPoints)
% % % % % % % %         if ~isequal(signs(turnPoints(kk)-1), signs(turnPoints(kk)-2))
% % % % % % % %             turnPoints2 = [turnPoints2; turnPoints(kk)];
% % % % % % % %         end
% % % % % % % %     end
% % % % % % % % end
    
turnPoints = find(diff(sign(diff(input_vector)))) + 1;

if isempty(turnPoints)
    isMonotonic = true;
else
    isMonotonic = false;
end
end % function