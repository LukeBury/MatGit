classdef (Sealed, CaseInsensitiveProperties=true, TruncatedProperties=true) mobiledev < hgsetget
%MOBILEDEV Read sensor data from mobile device running MATLAB Mobile.
%
%   Supported OS:
%   Android, Apple iOS
%
%   m = mobiledev() creates an object that reads sensor data from a mobile device
%   connected to the same network as the computer running MATLAB or MathWorks Cloud.
%
%
%   mobiledev methods:
%
%       Accessing logged data.
%           accellog - returns logged acceleration data
%           angvellog - returns logged angular velocity data
%           magfieldlog - returns logged magnetic field data
%           orientlog - returns logged orientation data       
%           poslog - returns logged position data
%
%       Discarding logged data.
%           discardlogs - discard all logged data
%
%   mobiledev properties:
%       Connected - Shows status of connection between MATLAB Mobile and mobiledev object in MATLAB
%       Logging - Shows and controls status of data transfer from device to MATLAB
%       InitialTimestamp - Time when first data point was transferred from 
%                          device to mobiledev in date format dd-mmm-yyyy HH:MM:SS.FFF.
%
%       Acceleration - Current acceleration reading: X, Y, Z in m/s^2
%       AngularVelocity - Current angular velocity reading: X, Y, Z in radians per second
%       Orientation - Current orientation reading: Azimuth, Pitch and Roll in degrees
%       MagneticField - Current magnetic field reading:  X, Y, Z in microtesla
%
%       Latitude - Current latitude reading in degrees
%       Longitude - Current longitude reading in degrees
%       Speed - Current speed reading in meters per second
%       Course - Current course reading in degrees relative to true north
%       Altitude - Current altitude reading in meters
%       HorizontalAccuracy - Current horizontal accuracy reading in meters
%
%       AccelerationSensorEnabled - Turns on/off accelerometer
%       AngularVelocitySensorEnabled - Turns on/off gyroscope
%       MagneticSensorEnabled - Turns on/off magnetometer
%       OrientationSensorEnabled - Turns on/off orientation sensor
%       PositionSensorEnabled - Turns on/off position sensor
%       SampleRate - Sets sample rate at which device will acquire the data
%
%   Usage
%
%   Before starting, if you are planning to use a computer to run MATLAB, 
%   connect your device to the same network as this computer. You may use Wi-Fi
%   or the cellular network. If you are planning to use MathWorks Cloud 
%   your device should have access to internet.
%
%   1. Start MATLAB Mobile.
%   2. Connect MATLAB Mobile to your computer or MathWorks Cloud. Refer to MATLAB
%      Mobile documentation for help. You need to do this step only
%      once.
%   3. In MATLAB, enter:  m = mobiledev() to create mobiledev object.
%
%   Access Data
%
%   You can get the latest value of a specific measurement by
%   querying the corresponding property. For example:
%
%       m.Acceleration
%
%   You can use mobiledev methods to access the logged measurement values.
%   For example, to get logged acceleration values:
%
%       [a, t] = accellog(m)
%  

% Copyright 2014-2018 The MathWorks, Inc.
    
    properties(GetAccess = public, SetAccess = private, Dependent = true)
        
        % InitialTimestamp - Time when first data point was transferred from 
        %                    device to mobiledev in date format dd-mmm-yyyy HH:MM:SS.FFF.
        InitialTimestamp
        
        % Acceleration - Current acceleration reading: X, Y, Z in m/s^2
        %
        % Acceleration is defined in relation to the X, Y and Z axes.
        %
        % If you set the device down face-up on a table, the positive X-axis
        % extends out of the right side of the device, positive Y-axis
        % extends out of the top side, and the positive Z-axis extends out
        % of the front face of the device. This is independent of the
        % orientation of the device.
        %
        % See also accellog
        Acceleration
        
        % AngularVelocity - Current angular velocity reading: X, Y, Z in 
        % radians per second.
        %
        % Angular velocity is defined in relation to the X, Y and Z axes and
        % in standard right-hand rotational vector notation.
        %
        % If you set the device down face-up on a table, the positive X-axis
        % extends out of the right side of the device, positive Y-axis
        % extends out of the top side, and the positive Z-axis extends out
        % of the front face of the device. This is independent of the
        % orientation of the device.
        %
        % See also angvellog
        AngularVelocity
        
        % MagneticField - Current magnetic field reading:  X, Y, Z in
        % microtesla.
        %
        % Magnetic field is defined in relation to the X, Y and Z axes.
        %
        % If you set the device down face-up on a table, the positive X-axis
        % extends out of the right side of the device, positive Y-axis
        % extends out of the top side, and the positive Z-axis extends out
        % of the front face of the device. This is independent of the
        % orientation of the device.
        %
        % See also magfieldlog
        MagneticField
        
        % Orientation -  Current orientation reading: Azimuth, Pitch and 
        % Roll in degrees.
        %
        % Orientation is defined in relation to the X, Y and Z axes.
        %
        % If you set the device down face-up on a table, the positive X-axis
        % extends out of the right side of the device, positive Y-axis
        % extends out of the top side, and the positive Z-axis extends out
        % of the front face of the device. This is independent of the
        % orientation of the device.
        %
        % Azimuth is the angle between the positive Y-axis and magnetic
        % north.
        %
        % Positive pitch is defined when the device starts by laying flat on
        % a surface and the positive Z-axis begins to tilt towards the
        % positive Y-axis.
        %
        % Positive roll is defined when the device starts by laying flat on
        % a surface and the positive Z-axis begins to tilt towards the
        % positive X-axis.
        %
        % Azimuth and roll range from -180 to 180 degrees.
        % Pitch ranges from -90 to 90 degrees.        
        %
        % See also orientlog
        Orientation
        
        % Latitude - Current latitude in degrees.
        %
        % Position data is obtained from GPS, Wi-Fi or cellular network,
        % whichever is available.
        %
        % See also poslog
        Latitude
        
        % Longitude - Current longitude in degrees.
        %
        % Position data is obtained from GPS, Wi-Fi or cellular network,
        % whichever is available.
        %
        % See also poslog
        Longitude

        % Speed - Current speed reading, in meters per second.
        Speed

        % Course - Current course reading, in degrees relative to true
        % north.
        Course
        
        % Altitude - Current altitude, in meters.
        %
        % Position data is obtained from GPS, Wi-Fi or cellular network,
        % whichever is available.
        %
        % See also poslog
        Altitude
        
        % HorizontalAccuracy - Current horizontal accuracy in meters.
        %
        % Position data is obtained from GPS, Wi-Fi or cellular network,
        % whichever is available.
        %
        % See also poslog
        HorizontalAccuracy
        
        % Connected - Shows status of connection between MATLAB Mobile and
        % mobiledev object in MATLAB.
        Connected        
    end
    
    properties(Access = public, Dependent = true, Transient = true)
        
        % AccelerationSensorEnabled - Turns on/off accelerometer
        AccelerationSensorEnabled
        
        % AngularVelocitySensorEnabled - Turns on/off gyroscope
        AngularVelocitySensorEnabled
        
        % MagneticSensorEnabled - Turns on/off magnetometer
        MagneticSensorEnabled
        
        % OrientationSensorEnabled - Turns on/off orientation sensor (synthetic sensor)
        OrientationSensorEnabled
        
        % PositionSensorEnabled - Turns on/off position sensor
        PositionSensorEnabled
        
        % Logging - Shows and controls status of data transfer from device to MATLAB
        Logging
        
        % Sample rate in Hz at which device attempts to acquire sensor data 
        SampleRate
    end
    
    properties(Access = private, Hidden, Transient = true)
        Controller        
    end
    
    properties(Constant, Access = private, Hidden)
        ConnectionTimeout = 5; % 5 sec
        
        ConnectionPollingPeriod  = 0.25; %0.25 sec
    end        
    
    properties(Access = public, Dependent = true, Hidden)
        MaxLogSize 
        Size
    end

    methods
        function obj = mobiledev(varargin)
            m = message('MATLAB:mobilesensor:mobiledev:ConnectorNotStarted');
            try
                m.getString();
            catch
                testSPPKGRootDirs = {};
                
                % registerrealtimecataloglocation
                vendorMFilePath = fileparts(mfilename('fullpath'));
                toolboxIndex = strfind(vendorMFilePath, [filesep, 'toolbox', filesep]);
                supportPackageBasePath = vendorMFilePath(1:toolboxIndex);
                
                testSPPKGRootDirs = [testSPPKGRootDirs, cellstr(supportPackageBasePath)];
                
                for idx = 1:numel(testSPPKGRootDirs)
                    mobilesensor.internal.Utility.registercataloglocation(testSPPKGRootDirs{idx});
                    try
                        m.getString();
                        break;
                    catch
                    end
                end
            end
            
            if ~mobilesensor.internal.MobileDevController.isConnectorOn()
                error(message('MATLAB:mobilesensor:mobiledev:ConnectorNotStarted'));
            end
                        
            try
                dataStruct = [];
                if nargin > 0
                    dataStruct = varargin{1};
                end
                obj.Controller = mobilesensor.internal.MobileDevController.getInstance(dataStruct);
                obj.Controller.register();
                
                for t = 0:obj.ConnectionPollingPeriod:obj.ConnectionTimeout
                    if obj.Controller.DeviceInitFinished
                        return;
                    end
                    pause(obj.ConnectionPollingPeriod);
                end
                                
            catch e
                if strcmp(e.identifier, 'MATLAB:mobilesensor:mobiledev:MultipleInstancesNotAllowed')
                    obj.Controller.warningNoTrace(message(e.identifier));
                else
                    throwAsCaller(e);
                end
            end
        end
                
        function discardlogs(obj)
            % DISCARDLOGS - Discards all logged measurements and InitialTimestamp
            %    discardlogs(obj)
            %
            % See also ACCELLOG, ANGVELLOG, MAGFIELDLOG, ORIENTLOG, POSLOG
            try
                obj.Controller.discardLogs();
            catch e
                throwAsCaller(e);
            end
        end
        
        function [a, t] = accellog(obj)
            % ACCELLOG - Returns logged acceleration data
            %    [a, t] = accellog(obj)
            %    a - [m x 3] matrix containing acceleration data points
            %    t - [m x 1] vector of timestamps, time in seconds 
            %        relative to InitialTimestamp
            %
            % See also ANGVELLOG, DISCARDLOGS, MAGFIELDLOG, ORIENTLOG, POSLOG
            [a, t] = obj.Controller.accellog();
        end
        
        
        function [a, t] = angvellog(obj)
            % ANGVELLOG - Returns logged angular velocity data
            %    [a, t] = angvellog(obj)
            %    a - [m x 3] matrix containing angular velocity data points
            %    t - [m x 1] vector of timestamps, time in seconds 
            %        relative to InitialTimestamp
            %
            % See also ACCELLOG, DISCARDLOGS, MAGFIELDLOG, ORIENTLOG, POSLOG
            [a, t] = obj.Controller.angvellog();
        end
        
        function [m, t] = magfieldlog(obj)
            % MAGFIELDLOG - Returns logged magnetic field data
            %    [m, t] = magfieldlog(obj)
            %    m - [m x 3] matrix containing magnetic field data points
            %    t - [m x 1] vector of timestamps, time in seconds 
            %        relative to InitialTimestamp
            %
            % See also ACCELLOG, ANGVELLOG, DISCARDLOGS, ORIENTLOG, POSLOG
            [m, t] = obj.Controller.magfieldlog();
        end
        
        function [o, t] = orientlog(obj)
            % ORIENTLOG - Returns logged orientation data
            %    [o, t] = orientlog(obj)
            %    o - [m x 3] matrix containing orientation data points
            %        Each row in matrix o represents azimuth, pitch and roll.
            %        If some values are unavailable, they are represented as NaN
            %    t - [m x 1] vector of timestamps, time in seconds 
            %        relative to InitialTimestamp
            %
            % See also ACCELLOG, ANGVELLOG, DISCARDLOGS, MAGFIELDLOG, POSLOG
            [o, t] = obj.Controller.orientlog();
        end
        
        
        function [lat, lon, t, speed, course, alt, horizacc] = poslog(obj)
            % POSLOG � Returns logged position data
            %    [lat, lon, t, speed, course, alt, horizacc] = poslog(obj)
            %    lat - [m x 1] vector of latitude values
            %    lon - [m x 1] vector of longitude values
            %    t - [m x 1] vector of timestamps, 
            %        time in seconds relative to InitialTimestamp
            %    speed - [m x 1] vector of speed values
            %    course - [m x 1] vector of course values
            %    alt - [m x 1] vector of altitude values
            %    horizacc - [m x 1] vector of horizontal accuracy values
            %    Position data is obtained from GPS, Wi-Fi or cellular
            %    network, which ever is available
            %
            % See also ACCELLOG, ANGVELLOG, DISCARDLOGS, MAGFIELDLOG, ORIENTLOG
            [lat, lon, t, speed, course, alt, horizacc] = obj.Controller.poslog();
        end                              
        
        function value = get.InitialTimestamp(obj)
            value = obj.Controller.getInitialTimestamp();
            
            if ~isempty(value)
                value = datestr(value, 'dd-mmm-yyyy HH:MM:SS.FFF');
            else 
                value = '';
            end
        end
        
        function value = get.Acceleration(obj)
            if obj.AccelerationSensorEnabled
                value = obj.Controller.getCurrentValueFor('Acceleration');
            else
                value = [];
            end
            
            if isempty(value)
                value = zeros(0, 3);
            end
            
        end
        
        function value = get.AngularVelocity(obj)
            if obj.AngularVelocitySensorEnabled
                value = obj.Controller.getCurrentValueFor('AngularVelocity');
            else
                value = [];
            end
            
            if isempty(value)
                value = zeros(0, 3);
            end
            
        end
        
        function value = get.Orientation(obj)
            if obj.OrientationSensorEnabled
                value = obj.Controller.getCurrentValueFor('Orientation');
            else
                value = [];
            end
            
            if isempty(value)
                value = zeros(0, 3);
            end
            
        end
        
        function value = get.MagneticField(obj)
            if obj.MagneticSensorEnabled
                value = obj.Controller.getCurrentValueFor('MagneticField');
            else
                value = [];
            end
            
            if isempty(value)
                value = zeros(0, 3);
            end
            
        end
        
        function value = get.HorizontalAccuracy(obj)
            if obj.PositionSensorEnabled
                value = obj.Controller.getCurrentValueFor('HorizontalAccuracy');
            else
                value = [];
            end
            
            if isempty(value)
                value = zeros(0, 1);
            end
        end
        
        function value = get.Latitude(obj)
            
            if obj.PositionSensorEnabled
                value = obj.Controller.getCurrentValueFor('Latitude');
            else
                value = [];
            end
            if isempty(value)
                value = zeros(0, 1);
            end
            
            
        end
        
        function value = get.Longitude(obj)
            
            if obj.PositionSensorEnabled
                value = obj.Controller.getCurrentValueFor('Longitude');
            else
                value = [];
            end
            if isempty(value)
                value = zeros(0, 1);
            end
            
        end
        
        function value = get.Altitude(obj)
            if obj.PositionSensorEnabled
                value = obj.Controller.getCurrentValueFor('Altitude');
            else
                value = [];
            end
            
            if isempty(value)
                value = zeros(0, 1);
            end
        end
        
        function value = get.Speed(obj)
            if obj.PositionSensorEnabled
                value = obj.Controller.getCurrentValueFor('Speed');
            else
                value = [];
            end
            if isempty(value)
                value = zeros(0, 1);
            end
        end
        
        function value = get.Course(obj)
            if obj.PositionSensorEnabled
                value = obj.Controller.getCurrentValueFor('Course');
            else
                value = [];
            end
            
            if isempty(value)
                value = zeros(0, 1);
            end
        end
        
        
        function value = get.Connected(obj)
            try
                value = obj.Controller.ConnectionState;
            catch e
                throwAsCaller(e);
            end
        end
        
        function value = get.Logging(obj)
            try
                value = obj.Controller.AcquisitionState;
            catch e
                throwAsCaller(e);
            end
        end
        
        function set.Logging(obj, isLogging)
            try
                isLogging = obj.getLogical(isLogging);
                if isLogging
                    obj.Controller.sendStart();
                else
                    obj.Controller.sendStop();
                end
            catch e
                throwAsCaller(e);
            end
        end

        function value = get.AccelerationSensorEnabled(obj)
            try
                value = obj.Controller.getSensorState('Acceleration');
            catch e
                throwAsCaller(e);
            end
        end
        
        function set.AccelerationSensorEnabled(obj, state)
            try
                obj.setSensorState('Acceleration', state);
            catch e
                throwAsCaller(e);
            end
        end
        
        function value = get.AngularVelocitySensorEnabled(obj)
            try
                value = obj.Controller.getSensorState('AngularVelocity');
            catch e
                throwAsCaller(e);
            end
        end
        
        function set.AngularVelocitySensorEnabled(obj, state)
            try
                obj.setSensorState('AngularVelocity', state);
            catch e
                throwAsCaller(e);
            end
        end
        
        function value = get.MagneticSensorEnabled(obj)
            try
                value = obj.Controller.getSensorState('MagneticField');
            catch e
                throwAsCaller(e);
            end
        end
        
        function set.MagneticSensorEnabled(obj, state)
            try
                obj.setSensorState('MagneticField', state);
            catch e
                throwAsCaller(e);
            end    
        end
        
        function value = get.OrientationSensorEnabled(obj)
            try
                value = obj.Controller.getSensorState('Orientation');
            catch e
                throwAsCaller(e);
            end
        end
        
        function set.OrientationSensorEnabled(obj, state)
            try
                obj.setSensorState('Orientation', state);
            catch e
                throwAsCaller(e);
            end    
        end
        
        function value = get.PositionSensorEnabled(obj)
            try
                value = obj.Controller.getSensorState('Position');
            catch e
                throwAsCaller(e);
            end
        end
        
        function set.PositionSensorEnabled(obj, state)
            try
                obj.setSensorState('Position', state);
            catch e
                throwAsCaller(e);
            end    
        end
        
        function value = get.SampleRate(obj)
            value = obj.Controller.SampleRate;
        end
        
        function set.SampleRate(obj, rate)
            try
                obj.Controller.updateDeviceSampleRate(rate);
            catch e
                throwAsCaller(e);
            end
            
        end        
       
        function value = get.Size(obj)
            value = obj.Controller.getTotalLogSize();
        end


        function value = get.MaxLogSize(obj)
            value = obj.Controller.MaxLogSize;
        end
        
        function set.MaxLogSize(obj, size)
            if (isnumeric(size))
                obj.Controller.MaxLogSize = size;
            end
        end
        
        function disp(obj)
            
            if ~isvalid(obj)
                error(message('MATLAB:class:InvalidHandle'));
            end

            doubleXYZFormat = '[%.4f %.4f %.4f] %s';
            doubleFormat = '%.4f %s';
            positionFormat = '%.6f %s';
            emptyXYZFormat = '[0x3 double] (Waiting for data)';
            emptyDoubleFormat = '[0x1 double] (Waiting for data)';
            
            [accelVal, accelLogSize] = obj.Controller.getCurrentValueFor('Acceleration');            
            accel = emptyXYZFormat;
            if ~isempty(accelVal)
                accel = sprintf(doubleXYZFormat, accelVal, '(m/s^2)');
            else
                accelLogSize = obj.Controller.getLogSizeFor('Acceleration');
            end
            
            [angvelVal, angvelLogSize] = obj.Controller.getCurrentValueFor('AngularVelocity'); 
            angvel = emptyXYZFormat;
            if ~isempty(angvelVal)
                angvel = sprintf(doubleXYZFormat, angvelVal, '(rad/s)');
            else
                angvelLogSize = obj.Controller.getLogSizeFor('AngularVelocity');
            end
            
            [magfieldVal, magfieldLogSize] = obj.Controller.getCurrentValueFor('MagneticField'); 
            magfield = emptyXYZFormat;
            if ~isempty(magfieldVal)
                magfield = sprintf(doubleXYZFormat, magfieldVal, '(microtesla)');
            else
                magfieldLogSize = obj.Controller.getLogSizeFor('MagneticField');
            end
            
            [orientVal, orientLogSize] = obj.Controller.getCurrentValueFor('Orientation');
            orient = emptyXYZFormat;
            if ~isempty(orientVal)
                orient = sprintf(doubleXYZFormat, orientVal, '(degrees)');
            else
                orientLogSize = obj.Controller.getLogSizeFor('Orientation');
            end
            
            [latitudeVal, posLogSize] = obj.Controller.getCurrentValueFor('Latitude');
            latitude = emptyDoubleFormat;
            if ~isempty(latitudeVal)
                latitude = sprintf(positionFormat, latitudeVal, '(degrees)');
            else
                posLogSize = obj.Controller.getLogSizeFor('Position');
            end
            
            longitudeVal = obj.Controller.getCurrentValueFor('Longitude');
            longitude = emptyDoubleFormat;
            if ~isempty(longitudeVal)
                longitude = sprintf(positionFormat, longitudeVal, '(degrees)');
            end
            
            altitudeVal = obj.Controller.getCurrentValueFor('Altitude');
            altitude = emptyDoubleFormat;
            if ~isempty(altitudeVal)
                altitude = sprintf(doubleFormat, altitudeVal, '(m)');
            end
            
            horizontalAccuracyVal = obj.Controller.getCurrentValueFor('HorizontalAccuracy');
            horizontalAccuracy = emptyDoubleFormat;
            if ~isempty(horizontalAccuracyVal)
                horizontalAccuracy = sprintf(doubleFormat, horizontalAccuracyVal, '(m)');
            end
            
            speedVal = obj.Controller.getCurrentValueFor('Speed');
            speed = emptyDoubleFormat;
            if ~isempty(speedVal)
                speed = sprintf(doubleFormat, speedVal, '(m/s)');
            end
            
            courseVal = obj.Controller.getCurrentValueFor('Course');
            course = emptyDoubleFormat;
            if ~isempty(courseVal)
                course = sprintf(doubleFormat, courseVal, '(degrees)');
            end                       
            
            loggedValuesStr = '    (%d Logged values)';
            
            mobiledevHelpStr = 'mobiledev';
            if matlab.internal.display.isHot()
                mobiledevHelpStr = '<a href="matlab:helpPopup mobiledev">mobiledev</a>';
            end
            
            loggingOnDeviceStr = ''; % If device "Log" option is selected on the device
                                     % show this info in the disp
            if obj.Logging && ~obj.Controller.isStreamToMATLABMode
                loggingOnDeviceStr = '    (Acquiring on device)';
            end
            
            fprintf('%s with properties:\n\n', mobiledevHelpStr);
            fprintf('                   Connected: %d\n', obj.Connected);
            fprintf('                     Logging: %d%s\n', obj.Logging, loggingOnDeviceStr);      
            
            if obj.Controller.isStreamToMATLABMode
                fprintf('            InitialTimestamp: ''%s''\n', obj.InitialTimestamp);
                fprintf('\n');
            end
            
            fprintf('   AccelerationSensorEnabled: %d', obj.AccelerationSensorEnabled);
            if accelLogSize > 0
                fprintf(loggedValuesStr, accelLogSize);
            end
            fprintf('\n');           
            
            fprintf('AngularVelocitySensorEnabled: %d', obj.AngularVelocitySensorEnabled);
            if angvelLogSize > 0
                fprintf(loggedValuesStr, angvelLogSize);
            end
            fprintf('\n');
                        
            fprintf('       MagneticSensorEnabled: %d', obj.MagneticSensorEnabled);
            if magfieldLogSize > 0
                fprintf(loggedValuesStr, magfieldLogSize);
            end
            fprintf('\n');
            
            fprintf('    OrientationSensorEnabled: %d', obj.OrientationSensorEnabled);
            if orientLogSize > 0
                fprintf(loggedValuesStr, orientLogSize);
            end            
            fprintf('\n');
            
            fprintf('       PositionSensorEnabled: %d', obj.PositionSensorEnabled);            
            if posLogSize > 0
                fprintf(loggedValuesStr, posLogSize);
            end
            fprintf('\n');
            if obj.Logging && obj.Controller.isStreamToMATLABMode
                fprintf('\n');
                fprintf('Current Sensor Values:\n');
                if obj.AccelerationSensorEnabled
                    fprintf('                Acceleration: %s\n', accel);
                end
                
                if obj.AngularVelocitySensorEnabled
                    fprintf('             AngularVelocity: %s\n', angvel);
                end
                
                if obj.MagneticSensorEnabled
                    fprintf('               MagneticField: %s\n', magfield);
                end
                
                if obj.OrientationSensorEnabled
                    fprintf('                 Orientation: %s\n', orient);
                end
                
                if obj.PositionSensorEnabled
                    fprintf('\n');
                    fprintf('         Position Data:\n');
                    fprintf('                    Latitude: %s\n', latitude);
                    fprintf('                   Longitude: %s\n', longitude);
                    fprintf('                       Speed: %s\n', speed);
                    fprintf('                      Course: %s\n', course);
                    fprintf('                    Altitude: %s\n', altitude);
                    fprintf('          HorizontalAccuracy: %s\n', horizontalAccuracy);
                end
            end
            
            if matlab.internal.display.isHot()
                supportedFunctionsStr = '<a href="matlab:methods mobiledev">Supported functions</a>';
                fprintf('\n%s\n', supportedFunctionsStr);
            end
        end

    end
    
    methods(Access = private)
 
        function logicalValue = getLogical(obj, aValue)
            if length(aValue) == 1 && (islogical(aValue) || isnumeric(aValue))
                logicalValue = logical(aValue);
            else
                error(message('MATLAB:mobilesensor:mobiledev:WrongLogicalInputValue'));
            end
        end
        
        function setSensorState(obj, sensorName, state)                        
            state = obj.getLogical(state);
            obj.Controller.updateDeviceSensorState(sensorName, state);
        end 
        
        function value = read(obj, name, n)
            if nargin < 3
                value = obj.Controller.readProperty(name);
            else
                value = obj.Controller.readProperty(name, n);
            end
        end
        
    end
    
    methods (Hidden)
        
        function delete(obj)
            % delete - Stops listening and frees all associated resources
            % delete(obj)
            try
                if (isvalid(obj.Controller))
                    obj.Controller.unregister();
                end
            catch e
                throwAsCaller(e);
            end
        end

        function s = saveobj(obj)
            s = obj.Controller.getAllValuesAsStruct();
        end
        
        % Disable and hide these methods.
        function c = horzcat(varargin)
            %HORZCAT Horizontal concatenation of mobiledev objects.
            
            if (nargin == 1)
                c = varargin{1};
            else
                error(message('MATLAB:mobilesensor:mobiledev:NoConcatenation'));
            end
        end
        function c = vertcat(varargin)
            %VERTCAT Vertical concatenation of mobiledev objects.
            
            if (nargin == 1)
                c = varargin{1};
            else
                error(message('MATLAB:mobilesensor:mobiledev:NoConcatenation'));
            end
        end
        function c = cat(varargin)
            %CAT Concatenation of mobiledev objects.
            if (nargin > 2)
                error(message('MATLAB:mobilesensor:mobiledev:NoConcatenation'));
            else
                c = varargin{2};
            end
        end

        % Hidden methods from the hgsetget super class.
        function res = eq(obj, varargin)
            res = eq@hgsetget(obj, varargin{:});
        end
        function res =  fieldnames(obj, varargin)
            res = fieldnames@hgsetget(obj,varargin{:});
        end
        function res = ge(obj, varargin)
            res = ge@hgsetget(obj, varargin{:});
        end
        function res = gt(obj, varargin)
            res = gt@hgsetget(obj, varargin{:});
        end
        function res = le(obj, varargin)
            res = le@hgsetget(obj, varargin{:});
        end
        function res = lt(obj, varargin)
            res = lt@hgsetget(obj, varargin{:});
        end
        function res = ne(obj, varargin)
            res = ne@hgsetget(obj, varargin{:});
        end
        function res = findobj(obj, varargin)
            res = findobj@hgsetget(obj, varargin{:});
        end
        function res = findprop(obj, varargin)
            res = findprop@hgsetget(obj, varargin{:});
        end
        function res = addlistener(obj, varargin)
            res = addlistener@hgsetget(obj, varargin{:});
        end
        function res = notify(obj, varargin)
            res = notify@hgsetget(obj, varargin{:});
        end        
        
        function getdisp(obj, varargin)
            getdisp@hgsetget(obj, varargin{:});
        end
        
        function setdisp(obj, varargin)
            setdisp@hgsetget(obj, varargin{:});
        end
        
        function get(obj, varargin)
            get@hgsetget(obj, varargin{:});
        end
        
        function set(obj, varargin)
            set@hgsetget(obj, varargin{:});
        end
       
    end
    
    methods (Static = true, Hidden = true)
       function out = loadobj(s)
           if isstruct(s)
              out = mobiledev(s);
           end
       end
       
       function ver = compatibilityVersion()
           ver = '1.0.0';
       end
    end
end

