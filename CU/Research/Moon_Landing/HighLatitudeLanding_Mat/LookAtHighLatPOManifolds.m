% ========================================================================
%%% Description
% ========================================================================
% For Ploltting the stable and unstable manifolds of periodic orbits in the
% CR3BP

% Created: 11/04/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
POfamilyPath = '~/CU_Google_Drive/Documents/MatGit/mbin/Data/InitialConditions/PO_Families/';
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

% ========================================================================
%%% Run Switches
% ========================================================================
plot_POs_alone = true;
    print_stability = true;

plot_unstable_manifolds = true;
plot_stable_manifolds   = true;

% ========================================================================
%%% Setup
% ========================================================================


% PO_0 = [1.0188945432097165; % 50 mps, from Hg2_2T (1.2e-4 to keep EV<1e-2 with 500 nodes)
%  0.0000084264050928;
%  0.0000248138790548;
%  0.0000820188150458;
%  -0.0019277484317759;
%  -0.0056774556726644;
%  8.5379815507225345];
% 
% PO_0 = [1.0188947224382787; % new 50mps, Hg2_2T %%%%%%%%%%%%%%%%%%%%%%
%  -0.0000000000000000;
%  -0.0000000000001439;
%  -0.0000000000004953;
%  -0.0019281746165965;
%  -0.0056777093361476;
%  8.5379815423875147];

% PO_0 = [1.0186054299379799;% 100 mps, from Hg2_2T (4.3e-4 to keep EV<1e-2 with 500 nodes)
%  -0.0000005334787993;
%  -0.0000036928818854;
%  -0.0000083258198245;
%  -0.0013196928865993;
%  -0.0091352799493384;
%  8.5904353683049912];

% PO_0 = [1.0180777740785210;% 150 mps, from Hg2_2T (4.3e-4 to keep EV<1e-2 with 500 nodes)
%  -0.0000000861992053;
%  -0.0000041953313344;
%  -0.0000077533285057;
%  -0.0002705149954614;
%  -0.0131661271511462;
%  8.6836127475510114];

% % 300mps, Hg2_2T
% PO_0 = [1.0152698354469363;
%  0.0000000000000000;
%  0.0000000000002841;
%  0.0000000000003682;
%  0.0040879775689879;
%  -0.0273925852552547;
%  9.2737612302749870];


% PO_0 = [1.0190308182870114; % 50 mps, from Hg2_3T (4.3e-4 to keep <10mps with 500 nodes)
%  0.0000000968206097;
%  0.0025020221301097;
%  0.0000007933183840;
%  -0.0018955785424653;
%  0.0000005758076360;
%  8.7158528096876111];
% 
% PO_0 = [1.0190308183070651; % new 50mps, Hg2_3T %%%%%%%%%%%%%%%%%%%%%%
%  0.0000000000000000;
%  0.0025020221451141;
%  0.0000000000003619;
%  -0.0018955785913221;
%  0.0000000000002428;
%  8.7158528198120102];

% PO_0 = [1.0188046165784910;% 100 mps, from Hg2_3T (4.3e-4 to keep <10mps with 500 nodes)
%  0.0000004688769755;
%  0.0040741741479739;
%  0.0000045794561843;
%  -0.0014286142403243;
%  0.0000059916355900;
%  8.7761778175414857];

% PO_0 = [1.0184135914987997;% 150 mps, from Hg2_3T (4.3e-4 to keep EV<1e-2 with 500 nodes)
%  -0.0000002981364929;
%  0.0058029263480480;
%  -0.0000056175795963;
%  -0.0006036488712580;
%  -0.0000127436302248;
%  8.8802317229718888];

% % 300mps, Hg2_3T
% PO_0 = [1.0159655791985782;
%  0.0000000000000000;
%  0.0108955925318208;
%  -0.0000000000000071;
%  0.0050466697632364;
%  -0.0000000000000398;
%  9.5238277843384402];


% PO_0 = [1.0180936955576518;% 50 mps, from Hg2_1P3
%  -0.0000004625151271;
%  -0.0000009222812682;
%  -0.0000039071534677;
%  -0.0036617222260194;
%  -0.0073000076291021;
%  22.1149413247433380];


% PO_0 = [1.0154991317649156;% 50 mps, from Se7_2T (4.3e-4 to keep EV<1e-2 with 500 nodes)
%  0.0000080370822492;
%  0.0000122484534463;
%  -0.0000360252574998;
%  0.0091061417750163;
%  0.0138776563068140;
%  9.9398777701788941];
% 
% PO_0 = [1.0154991476627844; % new 50mps, Se7_2T %%%%%%%%%%%%%%%%%%%%%%
%  0.0000000000000000;
%  -0.0000000000000354;
%  -0.0000000000000887;
%  0.0091061337833824;
%  0.0138776979896017;
%  9.9398777698981551];

% PO_0 = [1.0152368326685817;%  100mps, Se7_2T 
%  0.0000000000000000;
%  -0.0000000000010431;
%  0.0000000000020702;
%  0.0096884725865535;
%  0.0159216216057773;
%  10.0278144143953831];

% PO_0 = [1.0147925578457493;%  150mps, Se7_2T 
%  0.0000000000000000;
%  0.0000000000001684;
%  -0.0000000000011037;
%  0.0106579367211800;
%  0.0189691534357703;
%  10.1800949484084668];

% % 300mps, Se7_2T
% PO_0 = [1.0122062673132128;
%  0.0000000000000000;
%  -0.0000000000036441;
%  0.0000000000072322;
%  0.0159032005896660;
%  0.0330537405532260;
%  11.1339009665925044];

% PO_0 = [1.0153633231136696;% 50 mps, from Se7_3T (4.3e-4 to keep EV<1e-2 with 500 nodes)
%  0.0000012040901437;
%  0.0051228699076742;
%  -0.0000031310895275;
%  0.0099290101214616;
%  -0.0000042745282244;
%  9.9682370804001010];
% 
PO_0 = [1.0153633233036210; % new 50mps, Se7_3T %%%%%%%%%%%%%%%%%%%%%%
 0.0000000000000000;
 0.0051228701667772;
 -0.0000000000000656;
 0.0099290101709634;
 0.0000000000000456;
 9.9682370803616287];

% PO_0 = [1.0146667909509495; % 100 mps, Se7_3T
%  0.0018236536231968;
%  0.0052188464858317;
%  -0.0044983852675845;
%  0.0106373652823456;
%  -0.0067559902752775;
%  10.0503018109275875];
% 
% PO_0 = [1.0142780217657212;% 150 mps, Se7_3T(4.3e-4 to keep EV<1e-2 with 500 nodes)
%  0.0015707956689025;
%  0.0063214802403087;
%  -0.0029585949079481;
%  0.0122743273037210;
%  -0.0060312893819830;
%  10.1927243636450662];

% % 300mps, Se7_3T
% PO_0 = [1.0107139446857507;
%  0.0000000000000000;
%  0.0089879608530492;
%  0.0000000000018992;
%  0.0253447776784976;
%  -0.0000000000038923;
%  11.1086027519725139];

% 50 mps, Se7_5T
% PO_0 = [1.0191445186061201;
%  0.0000005196723063;
%  -0.0000011311151303;
%  -0.0000018104388537;
%  0.0022538324851705;
%  -0.0049056749220322;
%  13.0489905359020586];

% 50 mps, Se7_6T
% PO_0 = [1.0191958383652204;
%  0.0000006326690041;
%  -0.0022547431599923;
%  -0.0000016369625531;
%  0.0022777483886735;
%  0.0000027757749988;
%  13.1636496248989872];

% 50 mps, Se7_5P2
% PO_0 = [1.0175907727953430;
%  0.0000028203131541;
%  -0.0034358867501494;
%  -0.0000074666632023;
%  0.0053948530134847;
%  0.0000095584311252;
%  23.1020212034909278];
% 
% % % 50 mps, Se7_6P2
% PO_0 = [1.0175993645812185;
%  0.0000010398743246;
%  -0.0000016206236921;
%  -0.0000037985497322;
%  0.0052046621053846;
%  -0.0081113861158770;
%  23.1788228884912400];


% PO_0 = [0.9880134301736027;
%  0.0000000000000000;
%  -0.0168339215530371;
%  0.0000000000000000;
%  -0.0040733760074276;
%  0.0000000000000000;
%  4.2714646720470677];

%%% Manifold Colors
% color_stable   = colors.red;
% color_unstable = colors.grn;

color_stable   = colors.blue2;
color_unstable = colors.red2;

% -------------------------------------------------
%%% Options
% -------------------------------------------------
%%% Number of manifolds per PO
n_nodes = 500;
% n_nodes = 1200;

%%% Scale the perturbation of the manifold node in the unstable direction
% pertScale = 1.2e-4;
% pertScale = 4.3e-4;
pertScale = 2.1e-3;
% pertScale = 3e-3;
% pertScale = 4.3e-3;
% pertScale = 2e-2;

%%% Set propagation time for unstable manifolds
% Tf_manifolds_n = pi/2;
% Tf_manifolds_n = 1.965*pi;
% Tf_manifolds_n = pi*3.5;
% Tf_manifolds_n = 6*pi;
% Tf_manifolds_n = 8*pi;
% Tf_manifolds_n = 18*pi;
% Tf_manifolds_n = 25*pi;
% Tf_manifolds_n = 30*pi;
% Tf_manifolds_n = 40*pi;
Tf_manifolds_n = 50*pi;

% 300 nodes, 1e-4, 50pi
%%% 3B System
famName = 'Jupiter_Europa';

% -------------------------------------------------
%%% Set up System 
% -------------------------------------------------
%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);


%%% Normalizing constants
% rNorm:  n <-> km
% tNorm:  n <-> sec
% vNorm:  n <-> km/sec
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);


%%% Setting parameter structure
prms.u     = secondary.MR;
prms.rNorm = rNorm;
prms.R1    = primary.R / rNorm;
prms.R2    = secondary.R_n;
prms.n     = 1;

if contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
    prms.J2p = primary.J2; prms.J4p = primary.J4; prms.J6p = primary.J6; prms.J2s = secondary.J2;
    warning('Add something for prms.n')
    return
end

%%% Equillibrium Points
rLPs_n = collinearEquilibriumPoints_ZH(prms);

% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol                      = 1e-13;
options                  = odeset('RelTol',tol,'AbsTol',tol);
% options_XZStop           = odeset('Event',@event_yEqualsZeroPastL2,'RelTol',tol,'AbsTol',tol);
% options_impactOrL1Escape = odeset('Event',@event_ImpactorL1Escape_CR3Bn,'RelTol',tol,'AbsTol',tol);
% options_impact           = odeset('Event',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_zEquals0 = odeset('Event',@event_zEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
options_zEquals0_inNeck = odeset('Event',@event_zEqualsZero_inNeck, 'RelTol',tol,'AbsTol',tol);
    prms.xStop = 0.97;


% ========================================================================
%%% Find the manifolds and propagate them
% ========================================================================
% -------------------------------------------------
%%% Setup
% -------------------------------------------------
%%% STM vector
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);

% -------------------------------------------------
%%% Plot the original PO by itself
% -------------------------------------------------
if plot_POs_alone
    [T_PO, X_PO] = ode113(@Int_CR3BnSTM, [0, PO_0(7)], [PO_0(1:6); stm0_colVec], options, prms);
    
    figure; hold all
    plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'linewidth',2,'color',colors.mag)
%     plotTrajShadows(X_PO, 2, colors.grey, 'x', 0.965, 'y', 0.025, 'z', -4.5e-3)
    plotTrajShadows(X_PO, 2, colors.grey, 'x', 0.975, 'y', 0.025, 'z', -1e-2)
    axis equal
    ax = gca;
    ax.FontSize = 14;
    PlotBoi3_CR3Bn(26)
    plotSecondary(secondary)
    view(48,16)
    
    stm_tf_t0                           = reshape(X_PO(end,7:42),6,6);
    monodromy                           = stm_tf_t0;
    [eigenVectors_new, eigenValues_new] = eig(monodromy);
    [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
    fprintf('Stability Indices: %1.3f, %1.3f\n', S1, S2)
        
        
    if print_stability
        a = [real(diag(eigenValues_new)), imag(diag(eigenValues_new))];
        figure; hold all
        axis equal
        plotBody2(1, [0,0,0], colors.white, colors.black, 1, 0)
        plot(a(:,1),a(:,2),'o', 'markerfacecolor', colors.blue2, 'markeredgecolor', colors.black)
        PlotBoi2('Real', 'Imaginary', 26, 'LaTex')
    end
    
    latlons = zeros(size(X_PO,1),2);
    for kk = 1:size(X_PO,1)
        [lat_deg, lon_deg] = BCR2latlon(X_PO(kk,1:3), 'secondary', prms.u);
        latlons(kk,:) = [lat_deg, lon_deg];
    end
    % plot(latlons(:,2),latlons(:,1),'.','color',colors.mag,'markersize',7)
    if (abs(S1) < 2) && (abs(S2) < 2)
        fprintf('Stable - No manifolds\n')
        return
    end
end


% 989
% return

% -------------------------------------------------
%%% Loop through chosen indicies and compute manifolds at each
% -------------------------------------------------
% --------------------------
%%% Propagate initial conditions to get manifold node points
% --------------------------
%%% Integrate to get manifold nodes
[T_nodes, XSTM_nodes] = ode113(@Int_CR3BnSTM, linspace(0, PO_0(7), n_nodes+1), [PO_0(1:6); stm0_colVec], options, prms);

%%% Get the unstable eigenvector at the first node (initial condition)
monodromy_N0 = reshape(XSTM_nodes(end,7:42),6,6);
[eVec_stable_N0, eVec_unstable_N0] = getStableAndUnstableEigenvectors(monodromy_N0);

%%% Use state transition matrix to get stable and unstable eigenvectors 
%%% at each node
unstableEigenvectors_unit = NaN(6, n_nodes);
stableEigenvectors_unit   = NaN(6,n_nodes);
for node_i = 1:n_nodes
    %%% Grab STM from T0 to current node
    stm_Ni_T0                      = reshape(XSTM_nodes(node_i, 7:42),6,6);

    %%% Use STM to propagate the unstable eigenvector from N0 to the
    %%% current node
    unstableEigenvectors_unit(:,node_i) = stm_Ni_T0 * eVec_unstable_N0;
    stableEigenvectors_unit(:,node_i)   = stm_Ni_T0 * eVec_stable_N0;

    %%% Turn unstable eigenvector into a unit vector
    unstableEigenvectors_unit(:,node_i) = unstableEigenvectors_unit(:,node_i) ./ norm(unstableEigenvectors_unit(:,node_i));
    stableEigenvectors_unit(:,node_i)   = stableEigenvectors_unit(:,node_i) ./ norm(stableEigenvectors_unit(:,node_i));
end 


%%% Get rid of the repeat final/first point
X_nodes = XSTM_nodes(1:end-1,1:6);
T_nodes = T_nodes(1:end-1);

% --------------------------
%%% Preallocating for parfor
% --------------------------
man_unstable_p_POi = cell(n_nodes,1);
man_unstable_m_POi = cell(n_nodes,1);
man_stable_p_POi   = cell(n_nodes,1);
man_stable_m_POi   = cell(n_nodes,1);

man_unstable_p_zEquals0_POi = cell(n_nodes,1);
man_unstable_m_zEquals0_POi = cell(n_nodes,1);
man_stable_p_zEquals0_POi   = cell(n_nodes,1);
man_stable_m_zEquals0_POi   = cell(n_nodes,1);

% --------------------------
%%% Looping through manifolds ICs in parallel, computing the manifold
%%% direction, and propagating forward in time
% --------------------------
JCs = cell(n_nodes,1);
parfor node_i = 1:n_nodes
% for node_i = 1:n_nodes
    % --------------------------
    %%% Pick node and integrate PO with STM from it
    % --------------------------        
    %%% Set IC of current manifold node
    X0_manNode_i_n = X_nodes(node_i,:)';

    %%% Create the perturbed initial condition and integrate the manifold
    X0_man_unstable_p_i = X0_manNode_i_n + unstableEigenvectors_unit(:,node_i).*pertScale;
    X0_man_unstable_m_i = X0_manNode_i_n - unstableEigenvectors_unit(:,node_i).*pertScale;
    X0_man_stable_p_i   = X0_manNode_i_n + stableEigenvectors_unit(:,node_i).*pertScale;
    X0_man_stable_m_i   = X0_manNode_i_n - stableEigenvectors_unit(:,node_i).*pertScale;
    
    JCs{node_i}.un_p = getJacobiConstant_ZH(X0_man_unstable_p_i', prms);
    JCs{node_i}.un_m = getJacobiConstant_ZH(X0_man_unstable_m_i', prms);
    JCs{node_i}.st_p = getJacobiConstant_ZH(X0_man_stable_p_i', prms);
    JCs{node_i}.st_m = getJacobiConstant_ZH(X0_man_stable_m_i', prms);

    %%% Integrate the manifolds
%     [~, X_man_unstable_p] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_p_i, options, prms);
%     [~, X_man_unstable_m] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_m_i, options, prms);
% 
%     [~, X_man_stable_p] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_p_i, options, prms);
%     [~, X_man_stable_m] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_m_i, options, prms);
%     [~, X_man_unstable_p, time_event_unstable_p, X_event_unstable_p, index_event_unstable_p] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_p_i, options_zEquals0, prms);
%     [~, X_man_unstable_m, time_event_unstable_m, X_event_unstable_m, index_event_unstable_m] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_m_i, options_zEquals0, prms);
% 
%     [~, X_man_stable_p, time_event_stable_p, X_event_stable_p, index_event_stable_p] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_p_i, options_zEquals0, prms);
%     [~, X_man_stable_m, time_event_stable_m, X_event_stable_m, index_event_stable_m] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_m_i, options_zEquals0, prms);
%     
    [~, X_man_unstable_p, time_event_unstable_p, X_event_unstable_p, index_event_unstable_p] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_p_i, options_zEquals0_inNeck, prms);
    [~, X_man_unstable_m, time_event_unstable_m, X_event_unstable_m, index_event_unstable_m] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_m_i, options_zEquals0_inNeck, prms);

    [~, X_man_stable_p, time_event_stable_p, X_event_stable_p, index_event_stable_p] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_p_i, options_zEquals0_inNeck, prms);
    [~, X_man_stable_m, time_event_stable_m, X_event_stable_m, index_event_stable_m] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_m_i, options_zEquals0_inNeck, prms);
    
    
    %%% The integrator events include points where trajectories escape the
    %%% neck. Here, we discard those events, so the only events are z=0
    X_event_unstable_p    = X_event_unstable_p(abs(X_event_unstable_p(:,3)) < 1e-10, :);
    time_event_unstable_p = time_event_unstable_p(abs(X_event_unstable_p(:,3)) < 1e-10, :);
    
    X_event_unstable_m    = X_event_unstable_m(abs(X_event_unstable_m(:,3)) < 1e-10, :);
    time_event_unstable_m = time_event_unstable_m(abs(X_event_unstable_m(:,3)) < 1e-10, :);
    
    X_event_stable_p    = X_event_stable_p(abs(X_event_stable_p(:,3)) < 1e-10, :);
    time_event_stable_p = time_event_stable_p(abs(X_event_stable_p(:,3)) < 1e-10, :);
    
    X_event_stable_m    = X_event_stable_m(abs(X_event_stable_m(:,3)) < 1e-10, :);
    time_event_stable_m = time_event_stable_m(abs(X_event_stable_m(:,3)) < 1e-10, :);
    
    % --------------------------
    %%% Store manifolds for this PO
    % --------------------------
    man_unstable_p_POi{node_i} = X_man_unstable_p;
    man_unstable_m_POi{node_i} = X_man_unstable_m;
    man_stable_p_POi{node_i}   = X_man_stable_p;
    man_stable_m_POi{node_i}   = X_man_stable_m;
    
    man_unstable_p_zEquals0_POi{node_i} = [X_event_unstable_p(:,1:6), time_event_unstable_p];
    man_unstable_m_zEquals0_POi{node_i} = [X_event_unstable_m(:,1:6), time_event_unstable_m];
    man_stable_p_zEquals0_POi{node_i}   = [X_event_stable_p(:,1:6), time_event_stable_p];
    man_stable_m_zEquals0_POi{node_i}   = [X_event_stable_m(:,1:6), time_event_stable_m];
end

JCs_mat = zeros(n_nodes, 4);
JC_diffs_vec = zeros(n_nodes * 4, 1);
L2_excessVel_mat = zeros(n_nodes, 4);
DVs_vec = zeros(n_nodes*4,1);
Ev_vec = zeros(n_nodes*4,1);
JC_0 = getJacobiConstant_ZH(PO_0(1:6)', prms);

pos_un = zeros(n_nodes,1);
vel_un = zeros(n_nodes,1);
pos_st = zeros(n_nodes,1);
vel_st = zeros(n_nodes,1);

dJC_pos_un = zeros(n_nodes,1);
dJC_pos_st = zeros(n_nodes,1);
dJC_vel_un = zeros(n_nodes,1);
dJC_vel_st = zeros(n_nodes,1);

vel_angle_un = zeros(n_nodes,1);
vel_angle_st = zeros(n_nodes,1);

for node_i = 1:n_nodes
    JCs_mat(node_i, :) = [JCs{node_i}.un_p, JCs{node_i}.un_m, JCs{node_i}.st_p, JCs{node_i}.st_m];
    
    JC_diffs_vec(node_i*4-3: node_i*4) = JCs_mat(node_i, :)' - JC_0;
    
    L2_excessVel_mat(node_i, 1) = JC_2_L2FlyoverVelocity(JCs_mat(node_i, 1), prms, rLPs_n(2,:), vNorm);
    L2_excessVel_mat(node_i, 2) = JC_2_L2FlyoverVelocity(JCs_mat(node_i, 2), prms, rLPs_n(2,:), vNorm);
    L2_excessVel_mat(node_i, 3) = JC_2_L2FlyoverVelocity(JCs_mat(node_i, 3), prms, rLPs_n(2,:), vNorm);
    L2_excessVel_mat(node_i, 4) = JC_2_L2FlyoverVelocity(JCs_mat(node_i, 4), prms, rLPs_n(2,:), vNorm);
    
    DVs_vec(node_i*4-3: node_i*4)    = sqrt(abs(JCs_mat(node_i, :)' - JC_0)).*(vNorm*1000);
    Ev_vec(node_i*4-3: node_i*4) = sqrt(abs(JCs_mat(node_i, :)' - JC_0)).*vNorm;
    
    %%% For looking at relative sizes of position and velocity jumps for
    %%% manifolds
    X0_manNode_i_n = X_nodes(node_i,:)';
    
    X0_man_pos_un = X0_manNode_i_n + [unstableEigenvectors_unit(1:3,node_i); 0; 0; 0].*pertScale;
    X0_man_pos_st = X0_manNode_i_n + [stableEigenvectors_unit(1:3,node_i); 0; 0; 0].*pertScale;
    X0_man_vel_un = X0_manNode_i_n + [0; 0; 0; unstableEigenvectors_unit(4:6,node_i)].*pertScale;
    X0_man_vel_st = X0_manNode_i_n + [0; 0; 0; stableEigenvectors_unit(4:6,node_i)].*pertScale;
    
    JC_new_pos_un = getJacobiConstant_ZH(X0_man_pos_un', prms);
    JC_new_pos_st = getJacobiConstant_ZH(X0_man_pos_st', prms);
    JC_new_vel_un = getJacobiConstant_ZH(X0_man_vel_un', prms);
    JC_new_vel_st = getJacobiConstant_ZH(X0_man_vel_st', prms);
    
    dJC_pos_un(node_i) = JC_new_pos_un - JC_0;
    dJC_pos_st(node_i) = JC_new_pos_st - JC_0;
    dJC_vel_un(node_i) = JC_new_vel_un - JC_0;
    dJC_vel_st(node_i) = JC_new_vel_st - JC_0;
    
    pos_un(node_i) = norm(unstableEigenvectors_unit(1:3,node_i).*pertScale);
    vel_un(node_i) = norm(unstableEigenvectors_unit(4:6,node_i).*pertScale);
    pos_st(node_i) = norm(stableEigenvectors_unit(1:3,node_i).*pertScale);
    vel_st(node_i) = norm(stableEigenvectors_unit(4:6,node_i).*pertScale);
    
    
    vel_angle_un(node_i) = acos(dot(X0_manNode_i_n(4:6),X0_man_vel_un(4:6)) / (norm(X0_manNode_i_n(4:6))*norm(X0_man_vel_un(4:6))))*180/pi;
    vel_angle_st(node_i) = acos(dot(X0_manNode_i_n(4:6),X0_man_vel_st(4:6)) / (norm(X0_manNode_i_n(4:6))*norm(X0_man_vel_st(4:6))))*180/pi;

end
  
figure; hold all
p_pos = plot(pos_un, 'b', 'linewidth',2);
plot(pos_st, 'b', 'linewidth',2)
p_vel = plot(vel_un, 'r', 'linewidth',2);
plot(vel_st, 'r', 'linewidth',2)
legend([p_pos, p_vel],'Position ','Velocity','FontSize',14, 'location', 'best')
PlotBoi2('Node Number', 'Norm of step for manifold', 26, 'LaTex')

figure; hold all
histogram(JC_diffs_vec, 50,'Normalization','probability')
PlotBoi2('$\Delta JC$','Proportion',26, 'LaTex')

figure; hold all
histogram(Ev_vec, 50,'Normalization','probability')
PlotBoi2('$v_{N}\sqrt{|\Delta JC|}$','Proportion',26, 'LaTex')

figure; hold all
p_pos = plot(dJC_pos_un, 'b', 'linewidth',2);
plot(dJC_pos_st, 'b', 'linewidth',2)
p_vel = plot(dJC_vel_un, 'r', 'linewidth',2);
plot(dJC_vel_st, 'r', 'linewidth',2)
legend([p_pos, p_vel],'Position ','Velocity','FontSize',14, 'location', 'best')
PlotBoi2('Node Number', '$\Delta JC$ Due to component of manifold step', 26, 'LaTex')


figure; hold all
p_un = plot(vel_angle_un, 'color', colors.red2);
p_st = plot(vel_angle_st, 'color', colors.blue2);
PlotBoi2('Node Number', 'Angle between old and new velocities, $^\circ$',26,'LaTex')
legend([p_un, p_st],'Unstable Manifolds', 'Stable Manifolds')

fprintf('Minimum L2 excess velocity: %1.1f mps\n', min(L2_excessVel_mat,[],'all'))
fprintf('Maximum L2 excess velocity: %1.1f mps\n', max(L2_excessVel_mat,[],'all'))
% fprintf('Minimum DV: %1.1f mps\n', min(DVs_vec,[],'all'))
% fprintf('Maximum DV: %1.1f mps\n', max(DVs_vec,[],'all'))  
fprintf('Minimum Ev: %1.3e\n', min(Ev_vec))
fprintf('Maximum Ev: %1.3e\n', max(Ev_vec))  
    

    
figure; hold all
axis equal
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
view(0,90)
for man_i = 1:n_nodes
    if plot_unstable_manifolds
        p_u = plot3(man_unstable_p_POi{man_i}(:,1), man_unstable_p_POi{man_i}(:,2), man_unstable_p_POi{man_i}(:,3),'color', color_unstable,'linewidth',1);
        plot3(man_unstable_m_POi{man_i}(:,1), man_unstable_m_POi{man_i}(:,2), man_unstable_m_POi{man_i}(:,3),'color', color_unstable,'linewidth',1);
    end
    if plot_stable_manifolds
        p_s = plot3(man_stable_p_POi{man_i}(:,1), man_stable_p_POi{man_i}(:,2), man_stable_p_POi{man_i}(:,3),'color', color_stable,'linewidth',1);
        plot3(man_stable_m_POi{man_i}(:,1), man_stable_m_POi{man_i}(:,2), man_stable_m_POi{man_i}(:,3),'color', color_stable,'linewidth',1);
    end
end

if (plot_unstable_manifolds) && (~plot_stable_manifolds)
    legend(p_u, 'Unstable manifold','FontSize',14, 'location', 'best')
elseif (~plot_unstable_manifolds) && (plot_stable_manifolds)
    legend(p_s, 'Stable manifold','FontSize',14, 'location', 'best')
elseif (plot_unstable_manifolds) && (plot_stable_manifolds)
    legend([p_s, p_u], 'Stable manifold', 'Unstable manifold','FontSize',14, 'location', 'best')
end



figure; hold all
PlotBoi2('$r_2$','Inclination, $^\circ$', 26, 'LaTex')
xlim([0 0.016])


% figure(6395); hold all
% PlotBoi2('$r_2$','Inclination, $^\circ$', 26, 'LaTex')
% xlim([0 0.016])
R_i_unstable_m_stack = [];
R_i_unstable_p_stack = [];
R_i_stable_m_stack = [];
R_i_stable_p_stack = [];

for man_i = 1:n_nodes
    if plot_unstable_manifolds
        [X_SCI_unstable_m] = X_BaCR2SCI(man_unstable_m_zEquals0_POi{man_i}(:,1:6), man_unstable_m_zEquals0_POi{man_i}(:,7), prms);
        [X_SCI_unstable_p] = X_BaCR2SCI(man_unstable_p_zEquals0_POi{man_i}(:,1:6), man_unstable_p_zEquals0_POi{man_i}(:,7), prms);

        R_i_unstable_m = zeros(size(X_SCI_unstable_m,1),2);
        R_i_unstable_p = zeros(size(X_SCI_unstable_p,1),2);

        for kk = 1:size(X_SCI_unstable_m,1)
            [a,e,i,raan,w,ta] = ECI2OE(X_SCI_unstable_m(kk,1:3), X_SCI_unstable_m(kk,4:6),prms.u);
            R_i_unstable_m(kk,:) = [norm(X_SCI_unstable_m(kk,1:3)), i*180/pi];
        end

        for kk = 1:size(X_SCI_unstable_p,1)
            [a,e,i,raan,w,ta] = ECI2OE(X_SCI_unstable_p(kk,1:3), X_SCI_unstable_p(kk,4:6),prms.u);
            R_i_unstable_p(kk,:) = [norm(X_SCI_unstable_p(kk,1:3)), i*180/pi];
        end
        R_i_unstable_m_stack = [R_i_unstable_m_stack; R_i_unstable_m];
        R_i_unstable_p_stack = [R_i_unstable_p_stack; R_i_unstable_p];


        plot(R_i_unstable_m(:,1), R_i_unstable_m(:,2), 'r.')
        plot(R_i_unstable_p(:,1), R_i_unstable_p(:,2), 'r.')
        
    end
    
    
    if plot_stable_manifolds
        [X_SCI_stable_m] = X_BaCR2SCI(man_stable_m_zEquals0_POi{man_i}(:,1:6), man_stable_m_zEquals0_POi{man_i}(:,7), prms);
        [X_SCI_stable_p] = X_BaCR2SCI(man_stable_p_zEquals0_POi{man_i}(:,1:6), man_stable_p_zEquals0_POi{man_i}(:,7), prms);
        
        R_i_stable_m = zeros(size(X_SCI_stable_m,1),2);
        R_i_stable_p = zeros(size(X_SCI_stable_p,1),2);
        
        for kk = 1:size(X_SCI_stable_m,1)
            [a,e,i,raan,w,ta] = ECI2OE(X_SCI_stable_m(kk,1:3), X_SCI_stable_m(kk,4:6),prms.u);
            R_i_stable_m(kk,:) = [norm(X_SCI_stable_m(kk,1:3)), i*180/pi];
        end

        for kk = 1:size(X_SCI_stable_p,1)
            [a,e,i,raan,w,ta] = ECI2OE(X_SCI_stable_p(kk,1:3), X_SCI_stable_p(kk,4:6),prms.u);
            R_i_stable_p(kk,:) = [norm(X_SCI_stable_p(kk,1:3)), i*180/pi];
        end
        R_i_stable_m_stack = [R_i_stable_m_stack; R_i_stable_m];
        R_i_stable_p_stack = [R_i_stable_p_stack; R_i_stable_p];
%         figure(6395)
        plot(R_i_stable_m(:,1), R_i_stable_m(:,2), 'b.')
        plot(R_i_stable_p(:,1), R_i_stable_p(:,2), 'b.')
    end
    

end


figure; hold all
PlotBoi2('$r_2$','$|\dot{z}|$', 26, 'LaTex')
xlim([0 0.016])


for man_i = 1:n_nodes
    if plot_unstable_manifolds
        plot(rowNorm(man_unstable_m_zEquals0_POi{man_i}(:,1:3) - [1-prms.u,0,0]), abs(man_unstable_m_zEquals0_POi{man_i}(:,6)), 'r.')
        plot(rowNorm(man_unstable_p_zEquals0_POi{man_i}(:,1:3) - [1-prms.u,0,0]), abs(man_unstable_p_zEquals0_POi{man_i}(:,6)), 'r.')
    end

    if plot_stable_manifolds
        plot(rowNorm(man_stable_m_zEquals0_POi{man_i}(:,1:3) - [1-prms.u,0,0]), abs(man_stable_m_zEquals0_POi{man_i}(:,6)), 'b.')
        plot(rowNorm(man_stable_p_zEquals0_POi{man_i}(:,1:3) - [1-prms.u,0,0]), abs(man_stable_p_zEquals0_POi{man_i}(:,6)), 'b.')
    end

end





if 1+1==1
    for man_i = 1:n_nodes
        if plot_unstable_manifolds
            plot(rowNorm(man_unstable_m_zEquals0_POi{man_i}(:,1:3) - [1-prms.u,0,0]), abs(man_unstable_m_zEquals0_POi{man_i}(:,6)), '.', 'markersize', 12, 'color', colors.red2)
            plot(rowNorm(man_unstable_p_zEquals0_POi{man_i}(:,1:3) - [1-prms.u,0,0]), abs(man_unstable_p_zEquals0_POi{man_i}(:,6)), '.', 'markersize', 12, 'color', colors.red2)
        end

        if plot_stable_manifolds
            plot(rowNorm(man_stable_m_zEquals0_POi{man_i}(:,1:3) - [1-prms.u,0,0]), abs(man_stable_m_zEquals0_POi{man_i}(:,6)), '.', 'markersize', 12, 'color', colors.red2)
            plot(rowNorm(man_stable_p_zEquals0_POi{man_i}(:,1:3) - [1-prms.u,0,0]), abs(man_stable_p_zEquals0_POi{man_i}(:,6)), '.', 'markersize', 12, 'color', colors.red2)
        end

    end
end

%%%%%% Process for overlaying this data on grid-search data
% old_legend=findobj(gcf, 'Type', 'Legend');
% % (plot)
[legh,objh] = legend([old_legend.String,'Hg2\_2T Manifolds'], 'FontSize', 16);
lineh = findobj(objh,'type','line');
lineh = lineh([1, 3, 5]);
set(lineh,'linestyle','-','linewidth',3);
% plot([secondary.R_n secondary.R_n],[0 3],'k', 'linewidth',1.5); % zd
plot([secondary.R_n secondary.R_n],[0 60],'k', 'linewidth',1.5); % zd
%%%%%



if 1+1==1
    if plot_unstable_manifolds
        plot(R_i_unstable_m_stack(:,1), R_i_unstable_m_stack(:,2), '.', 'markersize', 12, 'color', colors.red2)
        plot(R_i_unstable_p_stack(:,1), R_i_unstable_p_stack(:,2), '.', 'markersize', 12, 'color', colors.red2)
    end

    if plot_stable_manifolds
        plot(R_i_stable_m_stack(:,1), R_i_stable_m_stack(:,2), '.', 'markersize', 12, 'color', colors.red2)
        plot(R_i_stable_p_stack(:,1), R_i_stable_p_stack(:,2), '.', 'markersize', 12, 'color', colors.red2)
    end

end


%%%%% Process for overlaying this data on grid-search data
% old_legend=findobj(gcf, 'Type', 'Legend');
% % % (plot)
% [legh,objh] = legend([old_legend.String,'Se7\_3T Manifolds'], 'FontSize', 16);
% lineh = findobj(objh,'type','line');
% lineh = lineh([1, 3, 5]);
% set(lineh,'linestyle','-','linewidth',3);
% plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5); % inclination
%%%%


if 1+1 == 1

    for man_i = 1:n_nodes
        latlons_un_m = zeros(size(man_unstable_m_POi{man_i},1),2);
        for kk = 1:size(man_unstable_m_POi{man_i},1)
            r = man_unstable_m_POi{man_i}(kk,1:3);
            [lat_deg, lon_deg] = BCR2latlon(r, 'secondary', prms.u);
            latlons_un_m(kk,:) = [lat_deg, lon_deg];
        end

        latlons_un_p = zeros( size(man_unstable_p_POi{man_i},1),2);
        for kk = 1:size(man_unstable_p_POi{man_i},1)
            r = man_unstable_p_POi{man_i}(kk,1:3);
            [lat_deg, lon_deg] = BCR2latlon(r, 'secondary', prms.u);
            latlons_un_p(kk,:) = [lat_deg, lon_deg];
        end

        latlons_st_m = zeros( size(man_stable_m_POi{man_i},1),2);
        for kk = 1:size(man_stable_m_POi{man_i},1)
            r = man_stable_m_POi{man_i}(kk,1:3);
            [lat_deg, lon_deg] = BCR2latlon(r, 'secondary', prms.u);
            latlons_st_m(kk,:) = [lat_deg, lon_deg];
        end

        latlons_st_p = zeros(size(man_stable_p_POi{man_i},1) ,2);
        for kk = 1:size(man_stable_p_POi{man_i},1)
            r = man_stable_p_POi{man_i}(kk,1:3);
            [lat_deg, lon_deg] = BCR2latlon(r, 'secondary', prms.u);
            latlons_st_p(kk,:) = [lat_deg, lon_deg];
        end
        
%         plot3(rowNorm(man_unstable_m_POi{man_i}(:,1:3) - [1-prms.u, 0, 0]),latlons_un_m(:,1), ones(size(latlons_un_m(:,1))).*5000, '.', 'markersize', 12, 'color', colors.red2)
%         plot3(rowNorm(man_unstable_p_POi{man_i}(:,1:3) - [1-prms.u, 0, 0]),latlons_un_p(:,1), ones(size(latlons_un_p(:,1))).*5000, '.', 'markersize', 12, 'color', colors.red2)
%         plot3(rowNorm(man_stable_m_POi{man_i}(:,1:3) - [1-prms.u, 0, 0]),latlons_st_m(:,1), ones(size(latlons_st_m(:,1))).*5000, '.', 'markersize', 12, 'color', colors.red2)
%         plot3(rowNorm(man_stable_p_POi{man_i}(:,1:3) - [1-prms.u, 0, 0]),latlons_st_p(:,1), ones(size(latlons_st_p(:,1))).*5000, '.', 'markersize', 12, 'color', colors.red2)

        plot(latlons_un_m(:,2),latlons_un_m(:,1), '.', 'markersize', 12, 'color', colors.red2)
        plot(latlons_un_p(:,2),latlons_un_p(:,1), '.', 'markersize', 12, 'color', colors.red2)
        plot(latlons_st_m(:,2),latlons_st_m(:,1), '.', 'markersize', 12, 'color', colors.red2)
        plot(latlons_st_p(:,2),latlons_st_p(:,1), '.', 'markersize', 12, 'color', colors.red2)
        
        plot(latlons_un_m(:,2),-latlons_un_m(:,1), '.', 'markersize', 12, 'color', colors.red2)
        plot(latlons_un_p(:,2),-latlons_un_p(:,1), '.', 'markersize', 12, 'color', colors.red2)
        plot(latlons_st_m(:,2),-latlons_st_m(:,1), '.', 'markersize', 12, 'color', colors.red2)
        plot(latlons_st_p(:,2),-latlons_st_p(:,1), '.', 'markersize', 12, 'color', colors.red2)
    end
end
% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
%%% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\n(Elapsed time: %1.4f seconds)\n',tocWhole)
















