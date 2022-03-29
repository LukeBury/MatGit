function [ JD ] = modifiedJulian2Julian(MJD)
%%% Description
%       Converts Julian Date to Modified Julian Date
% --------------------------------------------------------------
%%% Inputs
%       MJD - Modified Julian Date
% --------------------------------------------------------------
%%% Outputs
%       JD - Julian Date
% ===============================================================
JD = MJD + 2400000.5;
end









