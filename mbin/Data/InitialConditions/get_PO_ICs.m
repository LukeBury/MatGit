function [PO_ICs] = get_PO_ICs()
%%% Description
%       Returns initial conditions (X0 and Tp) for various types of
%       periodic orbits at various systems.
%
%       **NOTE: NO TOLERANCE GUARANTEED - MAY NEED TO SHOOT AND CORRECT
%       FROM THESE CONDITIONS DEPENDING ON NEEDS
%       
% ------------------------------------------------------------------------
%%% Inputs
%       
% ------------------------------------------------------------------------
%%% Outputs
%       PO_ICs -  {struct} struct with fields that correspend to various
%       systems
% ------------------------------------------------------------------------
% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% ========================================================================
%%% Earth-Moon
% ========================================================================
%%% L2 Lyapunov
PO_ICs.EarthMoon.CR3BP.L2_Lyapunov = [1.155656895105045;...
                                0;
                                0;
                                0;
                                0.000093165346366;
                                0;
                                3.373237818307939];
%%% L2 Northern Halo
PO_ICs.EarthMoon.CR3BP.L2_NHalo    = [1.180888549172871;...
                                -0.000000000000002;
                                -0.000137978545959;
                                0.000000000000128;
                                -0.155846477971307;
                                0;
                                3.415510771287679];                                    
%%% L2 Southern Halo
PO_ICs.EarthMoon.CR3BP.L2_SHalo    = [1.180888549172871;...
                                -0.000000000000002;
                                0.000137978545959;
                                0.000000000000128;
                                -0.155846477971307;
                                0;
                                3.415510771287679];   

% ========================================================================
%%% Jupiter-Europa
% ========================================================================
%%% L2 Lyapunov
PO_ICs.JupiterEuropa.CR3BP.L2_Lyapunov = [1.020387259845193;...
                                    0.000000000000080;
                                    0;
                                    0.000000000000051;
                                    0.000479945028703;
                                    0;
                                    3.076504388783907];
                                
%%% L2 Vertical
PO_ICs.JupiterEuropa.CR3BP.L2_Vertical = [1.020461555840290;...
                                    0.000000000000006;
                                    0;
                                    0.000000000000002;
                                    -0.000000295384993;
                                    -0.000288955625123;
                                    3.189373028060691];
                                
%%% L2 Northern Halo
PO_ICs.JupiterEuropa.CR3BP.L2_NHalo    = [1.017028544429315;...
                                    0;
                                    0.000023895588445;
                                    0;
                                    0.020011714717120;
                                    0;
                                    3.124212848627294];                                    
%%% L2 Southern Halo
PO_ICs.JupiterEuropa.CR3BP.L2_SHalo    = [1.017028544062047;...
                                    0;
                                    -0.000024035651262;
                                    0;
                                    0.020011715941382;
                                    0;
                                    3.124212845842041];   

%%% L2 Eastern Axial
PO_ICs.JupiterEuropa.CR3BP.L2_EasternAxial = [1.016217903229837;...
                                        0;
                                        -0.000000000000001;
                                        -0.000000000000001;
                                        -0.019672042942108;
                                        -0.055607116759711;
                                        4.266381393079572];   

%%% L2 Western Axial
PO_ICs.JupiterEuropa.CR3BP.L2_WesternAxial = [1.003798418598181;...
                                        0;
                                        -0.031825036366238;
                                        -0.002599833461840;
                                        0.013105232571822;
                                        -0.005032273682898;
                                        4.266374845956528];  

% ========================================================================
%%% Saturn-Enceladus
% ========================================================================
%%% L2 Lyapunov
PO_ICs.SaturnEnceladus.CR3BP.L2_Lyapunov = [1.003963326025942;...
                                      -0.000000000000001;
                                      0;
                                      -0.000000000000001;
                                      0.000184760844793;
                                      0;
                                      3.041677197603882];

%%% L2 Lyapunov
PO_ICs.SaturnEnceladus.CR3BP.L1_Lyapunov = [0.996033085409493;
                                            -0.000000000000000;
                                            0.000000000000000;
                                            0.000000000004301;
                                            -0.000095367933369;
                                            0.000000000000000;
                                            3.024489705865199]; 
                                
%%% L2 Vertical
PO_ICs.SaturnEnceladus.CR3BP.L2_Vertical = [1.003991359461843;...
                                      0.000000000000001;
                                      0;
                                      -0.000000000000001;
                                      -0.000000172368141;
                                      -0.000099389248238;
                                      3.151015108553384];

%%% L1 Vertical
PO_ICs.SaturnEnceladus.CR3BP.L1_Vertical = [0.996019267540287;...
                                      0.000000000000002;
                                      0;
                                      0;
                                      0.000000995749908;
                                      0.000240920423191;
                                      3.132417557725700];
%%% L2 Northern Halo
PO_ICs.SaturnEnceladus.CR3BP.L2_NHalo  = [1.004446422097840;...
                                    0;
                                    -0.000777541208133;
                                    0.000000000000305;
                                    -0.003759817746713;
                                    -0.000000000000041;
                                    3.086351054759526];                                

%%% L2 Southern Halo
PO_ICs.SaturnEnceladus.CR3BP.L2_SHalo  = [1.004446422097840;...
                                    0;
                                    0.000777541208133;
                                    0.000000000000305;
                                    -0.003759817746713;
                                    0.000000000000041;
                                    3.086351054759526];    
                                
%%% L2 Northern Halo
PO_ICs.SaturnEnceladus.CR3BP.L1_NHalo  = [0.995560588991134;
                                          0.000000000000000;
                                          -0.000482688805533;
                                          0.000000000000728;
                                          0.003610891192393;
                                          0.000000000000126;
                                          3.071808037967654];                           

%%% L2 Southern Halo
PO_ICs.SaturnEnceladus.CR3BP.L1_SHalo  = [0.995560588991134;
                                          0.000000000000000;
                                          0.000482688805533;
                                          0.000000000000728;
                                          0.003610891192393;
                                          -0.000000000000126;
                                          3.071808037967654];
                                                                
%%% DRO
PO_ICs.SaturnEnceladus.CR3BP.DRO  = [0.998312974740448;...
                                    0;
                                    0;
                                    -0.000000000000460;
                                    0.012505367434921;
                                    0;
                                    0.885611959509138];     
                                       
%%% L2 Eastern Axial
PO_ICs.SaturnEnceladus.CR3BP.L2_EasternAxial = [1.000817566039729;...
                                        0;
                                        0.006260935119950;
                                        -0.000135997227746;
                                        0.002501263457997;
                                        0.000273281110328;
                                        4.248655828819825];  
                                    
%%% L2 Western Axial
PO_ICs.SaturnEnceladus.CR3BP.L2_WesternAxial = [1.003125755237350;...
                                        0;
                                        0;
                                        0.000000000000004;
                                        -0.003490623466145;
                                        0.011041263805649;
                                        4.248656068496835];  
                                    
%%% Unlabeled 1 Trailing (south - kind of like a Halo, but rotated 90 
%%% degrees and centered about the secondary
PO_ICs.SaturnEnceladus.CR3BP.unlabeled1TrailingS = [1.007302090959311;...
                                                    0;
                                                    -0.004886898788232;
                                                    0.009376480285938;
                                                    -0.015363209502224;
                                                    0.004755187845638;
                                                    5.881033454492806];    
                                
%%% Unlabeled 1 Leading (south) - kind of like a Halo, but rotated 90 
%%% degrees and centered about the secondary
PO_ICs.SaturnEnceladus.CR3BP.unlabeled1LeadingS = [0.993453790361993;...
                                                   0;
                                                   -0.004042085896919;
                                                   -0.008637622615002;
                                                   0.013981142614446;
                                                   0.004616467492074;
                                                   5.722382510772366];   
                                       
% From unlabeled1TrailingS Biffurcation ... this is symmetric with that
% about z = 0 plane
PO_ICs.SaturnEnceladus.CR3BP.unlabeled1LeadingN = [0.996374878694517;
                                                   0.000000000000000;
                                                   0.001349243394813;
                                                   0.007526745713117;
                                                   0.009848794483643;
                                                   0.004560033495253;
                                                   4.257566451624301];

%%% Unlabeled 1 Trailing (North - kind of like a Halo, but rotated 90 
%%% degrees and centered about the secondary
PO_ICs.SaturnEnceladus.CR3BP.unlabeled1TrailingN = [0.996374878694517;
                                                    0.000000000000000;
                                                    0.001349243394813;
                                                    -0.007526745713117;
                                                    0.009848794483643;
                                                    -0.004560033495253;
                                                    4.257566451624301];    

%%% Some kind of double-period L1 Halo at Enceladus
PO_ICs.SaturnEnceladus.CR3BP.L1_NHalo_DoublePeriod = [0.999996805811347;
                                                      -0.000000000000000;
                                                      0.000971423431191;
                                                      0.000000000000002;
                                                      -0.018821623815005;
                                                      0.000000000000001;
                                                      4.259671208643186];   
                                                  
                                                  
%%% Some kind of L1 Northern Pretzel orbit at Enceladus
%* May have useful manifolds for landing on Enceladus south pole
PO_ICs.SaturnEnceladus.CR3BP.L1_NHalo_DoublePeriod = [0.999949172678353;
                                                      -0.000011913400115;
                                                      0.001050500079262;
                                                      -0.000018452263800;
                                                      -0.017998820732015;
                                                      -0.000114268319427;
                                                       4.257875322097243];                                                    
                                                  
%%% Double period solution from bifurcation of unlabeled1LeadingS near 
%%% enceladus             
PO_ICs.SaturnEnceladus.CR3BP.SLeadingSneakerToe = [1.003278945892397;
                                                   -0.000000000000000;
                                                   -0.001160178622759;
                                                   -0.008518652394546;
                                                   -0.009390206737848;
                                                   -0.004730702228170;
                                                   8.797668183419933];      
%%% Double period solution from bifurcation of unlabeled1LeadingS near 
%%% enceladus             
PO_ICs.SaturnEnceladus.CR3BP.SLeadingSneakerHeel = [1.003216735533549;
                                                    -0.000000000000000;
                                                    -0.001405927038742;
                                                    -0.008312954460291;
                                                    -0.009308705525712;
                                                    -0.003780793737525;
                                                    7.953030660306376];                                           

%%% Double period solution from bifurcation of unlabeled1LeadingS near 
%%% enceladus             
PO_ICs.SaturnEnceladus.CR3BP.unknown1 = [0.939388822542880;
                                       -0.409996929538443;
                                       0.077223984876444;
                                       -0.099807452215460;
                                       -0.020434635060696;
                                       -0.071232233570179;
                                       6.697427552781499];                                           
                                           

% ========================================================================
%%% Saturn-Titan
% ========================================================================
%%% L2 Lyapunov
PO_ICs.SaturnTitan.CR3BP.L2_Lyapunov = [1.043225594807153;...
                                      0;
                                      0;
                                      0;
                                      0.000184026217029;
                                      0;
                                      3.124245125264076];

%%% L2 Vertical
PO_ICs.SaturnTitan.CR3BP.L2_Vertical = [1.043254890948189;...
                                      -0.000000000000002;
                                      0;
                                      -0.000000000000001;
                                      -0.000000062689698;
                                      -0.000188867650859;
                                      3.241966380524457];   



% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------
% --------------------------
% 
% --------------------------


end % function