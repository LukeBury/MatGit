function [ percChange ] = percentchange( vec )
%%% Use:
%   Returns percent-change vector for input vector
%
%%% Inputs:
%   1) row [1xn] or column [nx1] vector
%
%%% Outputs:
%   1) percent change vector in same shape as input
% ========================================================================
if vec(1) == 0
    warning('percentchange.m was passed a vector starting with 0')
end
percChange = 100*(vec - vec(1))./vec(1);

end

