% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
savePath = '~/CU_Google_Drive/Documents/MatGit/mbin/Data/InitialConditions/PO_Families/';
addpath(genpath(mbinPath))
ticWhole = tic;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

%%% Periodic orbit ICs
PO_ICs = get_PO_ICs();


% -------------------------------------------------
%%% Choose data
% -------------------------------------------------
% famName = 'Jupiter_Europa.CR3BP.L2_Lyapunov';
% famName = 'Jupiter_Europa.CR3BP.L2_NHalo';
% famName = 'Jupiter_Europa.CR3BP.L2_SHalo';
% famName = 'Jupiter_Europa.CR3BP.L2_Vertical';
% famName = 'Jupiter_Europa.CR3BP.L2_L_2T'; % Full Axial family (Also L2_V_T)
% famName = 'Jupiter_Europa.CR3BP.L2_L_T_P2'; % Butterly?
famName = 'Jupiter_Europa.CR3BP.L2_L_P2'; % 



family = [famName,'.txt'];


% --------------------------
% Actually load family
% --------------------------
%%% Path from mbin to data
dataPathFromMBin = '/Data/InitialConditions/PO_Families/';

%%% PO data file
PO_datafile = [mbinPath, dataPathFromMBin, family];

% -------------------------------------------------
%%% Set up parameters
% -------------------------------------------------
%%% Set primary and secondary bodies
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(family, bodies);

%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% prms for integration
prms.u    = secondary.MR;
prms.R2 = secondary.R_n;
prms.n    = 1;

rLPs_n = EquilibriumPoints(prms.u, prms.n);

%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Load data
% -------------------------------------------------
%%% Load the data file
PO_data = dlmread(PO_datafile,',',1,0);

%%% Grab header line
fid = fopen(PO_datafile, 'rt');  %the 't' is important!
header = fgetl(fid);
fclose(fid);

%%% Number of ICs
n_POs = size(PO_data,1);

% --------------------------
% Create column specifiers
% --------------------------
PO_header_2020 = 'x0,y0,z0,xd0,yd0,zd0,Tp,JC,stabilityIndex1,stabilityIndex2,alpha,beta,impactFlag';

if contains(header, PO_header_2020)
    c_x0 = 1;   c_y0 = 2;   c_z0 = 3;
    c_xd0 = 4;  x_yd0 = 5;  c_zd0 = 6;
    c_Tp = 7;   c_JC = 8;   c_S1 = 9;   c_S2 = 10;
    c_alpha = 11;   c_beta = 12;    c_impactFlag = 13;
end












% -------------------------------------------------
%%% Preparing Save File
% -------------------------------------------------
%%% Create unique filename
[uniqueFilename] = get_uniqueFilename(famName, savePath, 'txt');

%%% Open File
datafile = fopen(uniqueFilename,'wt');

% -------------------------------------------------
%%% Writing data
% -------------------------------------------------
%%% Write header
headerString = ['x0,y0,z0,xd0,yd0,zd0,Tp,JC,stabilityIndex1,stabilityIndex2,alpha,beta,impactFlag,error\n'];
fprintf(datafile,headerString);
    
actualErrors = zeros(n_POs,1);    
actualErrors2 = zeros(n_POs,1);    
for index = 1:n_POs
    stm0_colVec = reshape(eye(6),36,1);
    [T_PO, X_PO] = ode113(@Int_CR3BnSTM, [0, PO_data(index,c_Tp)], [PO_data(index,c_x0:c_zd0)'; stm0_colVec], options, prms);
%     [T_PO, X_PO] = ode113(@Int_CR3Bn, [0, PO_data(index,c_Tp)], PO_data(index,c_x0:c_zd0)', options, prms);
    
%     actualError = norm(X_PO(end,1:6)' - X_PO(1,1:6)');
%     actualErrors(index) = actualError;
    
    actualError2 = norm(X_PO(end,1:6)' - X_PO(1,1:6)') / norm(X_PO(1,1:6));
    actualErrors2(index) = actualError2;
    
    if PO_data(index,c_y0)==0
        fprintf(datafile,'%1.16f,%1.1f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.4f,%1.4f,%1.4f,%1.4f,%1d,%1.1e\n',...
            PO_data(index,1), PO_data(index,2), PO_data(index,3), PO_data(index,4), PO_data(index,5), PO_data(index,6),...
            PO_data(index,7), PO_data(index,8), PO_data(index,9), PO_data(index,10), PO_data(index,11), PO_data(index,12),...
            PO_data(index,13), actualError2);
    else
        fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.4f,%1.4f,%1.4f,%1.4f,%1d,%1.1e\n',...
            PO_data(index,1), PO_data(index,2), PO_data(index,3), PO_data(index,4), PO_data(index,5), PO_data(index,6),...
            PO_data(index,7), PO_data(index,8), PO_data(index,9), PO_data(index,10), PO_data(index,11), PO_data(index,12),...
            PO_data(index,13), actualError2);
    end
end
%%% Close file
fclose(datafile);

% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\nElapsed time: %1.4f seconds\n',tocWhole)
















