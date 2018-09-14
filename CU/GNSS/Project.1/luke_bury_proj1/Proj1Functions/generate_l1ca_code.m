function code = generate_l1ca_code(prn)
    % Generates GPS L1 C/A code for given PRN.
    
    % Parameters
    % ----------------------------------------------------------------------------
    % prn : int 
    %     the signal PRN

    % Returns
    % ----------------------------------------------------------------------------
    % output : ndarray of shape(1023,)
    %     the complete code sequence
    %  The following obtains the phase select indices particular to the G1 MLS
    %  used to generate this PRN's Gold code.
    load('L1_CODE_PHASE_ASSIGNMENTS.mat');
    phaseSelector = L1_CODE_PHASE_ASSIGNMENTS{prn}{3};
    
    % Converting from cell array to regular array
    phaseSelector = cell2mat(phaseSelector);
    
    % Creating G1 and G2 MLS
    G1_MLS = generate_mls(10,[3 10],[10]);
    G2_MLS = generate_mls(10,[2,3,6,8,9,10],phaseSelector);
    
    % Creating Gold Code from G1 and G2
    code = mod(G1_MLS + G2_MLS,2);

end
