function [uniqueFilename] = get_uniqueFilename(filename, savePath, fileType)
%%% Description
% Generates a unique filename by placing an integer at the end
%       
% ------------------------------------------------------------------------
%%% Inputs
% filename - [string] Desired name of file
% savePath - [string] Path to save directory (ex: '~/Documents/data/')
% fileType - [string] Type of file (ex: 'txt' or 'csv')
% ------------------------------------------------------------------------
%%% Outputs
% uniqueFilename - [string] Unique filename with path
% ------------------------------------------------------------------------
% Created: 12/15/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialize file version and updated filename
fileVersion  = 0;
filename_new = filename;

%%% Add integers to the end of the filename until it's unique
while isfile([savePath,filename_new, '.', fileType])
    fileVersion = fileVersion + 1;
    filename_new = sprintf('%s_(%1d)', filename, fileVersion);
end

%%% Create 
uniqueFilename = [savePath, filename_new, '.', fileType];

end % function