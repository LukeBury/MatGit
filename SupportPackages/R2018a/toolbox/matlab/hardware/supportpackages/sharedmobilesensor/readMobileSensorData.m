function data = readMobileSensorData(fname, varargin)
%readMobileSensorData Imports sensor data from file previously collected by
%MATLAB Mobile
%   S = readMobileSensorData(FILENAME) will read file FILENAME, extract sensor
%   data from it and create structure s, where each field in this structure
%   is a timetable with sensor data.
%   If corresponding sensor data is present in the file structure s will
%   have the following fields:
%   Acceleration - timetable with acceleration data
%   AngularVelocity - timetable with angular velocity data
%   MagneticField - timetable with magnetic field data
%   Orientation - timetable with orientation data
%   Position - timetable with position data
%   
%   For example:
%   S = readMobileSensorData('sensorlog.zip');
%   acceleration = S.Acceleration;
%   plot(acceleration.Timestamp, acceleration{:, 2:end});

%   Copyright 2017-2018 MathWorks, Inc.

try
    unzipFolder = tempdir;
    if nargin > 1
        validateattributes(varargin{1}, {'char', 'string'}, {}, 2)
        unzipFolder = varargin{1};
    end
    
    % if .zip is missing - add it
    if ~endsWith(fname, ".zip")
        fname = strcat(fname, ".zip");
    end
    
    epochTime = datenum(1970,1,1,0,0,0);
    secondsInDay = 86400;
    
    % unzip archive
    dataFiles = unzip(fname, unzipFolder);
    
    % go over extracted data files and for each file extract the protion 
    % of its name which tells what sensor it is.
    % then read data into the table and save this table as a struct field 
    data = struct();
    for idx=1:length(dataFiles)
        dataFile = dataFiles{idx};
        [~, fn] = fileparts(dataFile);
        nameParts = split(fn, '_');
        sensorName = nameParts{2};
        tbl = readtable(dataFile, 'ReadVariableNames', true);
        tbl.Timestamp = datetime(tbl.timestamp / 1000 / secondsInDay + epochTime, 'ConvertFrom', 'datenum');
        tbl.Timestamp.Format = 'dd-MMM-uuuu HH:mm:ss.SSS';
        tbl.timestamp = [];
        fieldName = 'unknown';
        switch (sensorName)
            case 'accel'
                fieldName = 'Acceleration';
            case 'angvel'
                fieldName = 'AngularVelocity';
            case 'magfield'
                fieldName = 'MagneticField';
            case 'orient'
                fieldName = 'Orientation';
            case 'pos'
                fieldName = 'Position';
        end
        data.(fieldName) = table2timetable(tbl);
        delete(dataFiles{idx});
    end
    
catch e
    throwAsCaller(e);
end

end
