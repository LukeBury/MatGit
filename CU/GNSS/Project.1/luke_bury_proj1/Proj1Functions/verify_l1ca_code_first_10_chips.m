%X  This function will be used to check the correctness of `generate_l1ca_code`.
%X  Don't worry if this code seems confusing -- all it does is check that the
%X  first 10 bits of the generated Gold code match their supposed octal 
%X  representation according to the ICD specifications.
function res = verify_l1ca_code_first_10_chips(prn)
    code = generate_l1ca_code(prn);
    first_10_bits = code(1:10);
    
    % Converting binary vector to octal
    first_10_bits_oct = bin2oct(first_10_bits);

    load('L1_CODE_PHASE_ASSIGNMENTS.mat');
    octal = L1_CODE_PHASE_ASSIGNMENTS{prn}{7};

    res = first_10_bits_oct == octal;
end
