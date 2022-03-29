% ========================================================================
%%% Description
% ========================================================================
% Template for use of parfor in Matlab

% Created: 07/12/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Template
% ========================================================================

%%% Variables to be used
var1 = 9;

%%% Loop options/parameters
n_parforLoops = 10;
nOtherConditions = 6;

%%% Conditions to be looped through
X0_data = zeros(n_parforLoops,nOtherConditions);
otherConditions = linspace(1,13,nOtherConditions);


%%% Create data structure
parfor_data = {};

warning('Something like this ... Not verified if data is actually stored properly here... ')
for parforLoopIndex = 1:n_parforLoops
    % -------------------------------------------------
    % Reducing broadcast variables
    % -------------------------------------------------
    var1_parfor            = var1;
    otherConditions_parfor = otherConditions;
    
    % -------------------------------------------------
    % Grab outer loop data
    % -------------------------------------------------
    X0_i = X0_data(parforLoopIndex,:);
        
    % -------------------------------------------------
    % Preallocate data structure for each loop of otherConditions
    % -------------------------------------------------
    resultsFromX0_i = zeros(nOtherConditions,6);
    
    % -------------------------------------------------
    % Inner loop
    % -------------------------------------------------
    for innerLoopIndex = 1:nOtherConditions
        %%% Preallocating temporary variables
        result_i = [];
        
        %%% Calculate results
        result_i = otherConditions_parfor(innerLoopIndex) + X0_i(2) + var1_parfor*parforLoopIndex;
        
        %%% Store results
        resultsFromX0_i(innerLoopIndex,:) = result_i;

    end % burnMag_kk = 1:n_DVs_per_location
    
    % -------------------------------------------------
    % Storing data from all inner loops
    % -------------------------------------------------
    parfor_data{parforLoopIndex}.resultsFromX0 = resultsFromX0_i;
    
end % parfor 

% ========================================================================
%%% Post-Processing
% ========================================================================
% -------------------------------------------------
% Combine Data
% -------------------------------------------------
%%% Turning data in one large struct array
parfor_data_structArray = [parfor_data{:}];

%%% Seperating fields into vectors & matrices
resultsFromX0 = [parfor_data_structArray(:).resultsFromX0];


989

resultsFromX0          = reshape(resultsFromX0,6,100);
resultsFromX0          = resultsFromX0';








