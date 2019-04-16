function [ MJD ] = julian2modifiedJulian(JD)
%%% Description
%       Converts Julian Date to Modified Julian Date
% --------------------------------------------------------------
%%% Inputs
%       JD - Julian Date
% --------------------------------------------------------------
%%% Outputs
%       MJD - Modified Julian Date
% ===============================================================
MJD = JD - 2400000.5;
end









