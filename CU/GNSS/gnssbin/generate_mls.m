function values = generate_mls(N, feedback_taps, output_taps)
    % Generates maximum-length sequence (MLS) for the given linear feedback
    % shift register (LFSR) length, feedback taps, and output taps.  The initial
    % state of the LFSR is taken to be all ones.
    % 
    % Parameters
    % ----------------------------------------------------------------------------
    % N : int
    %     length of LFSR
    % feedback_taps : array or ndarray of shape (L,)
    %     the L taps to use for feedback to the shift register's first value
    % output_taps : array or ndarray of shape (M,)
    %     the M taps to use for choosing the code output
    % 
    % Returns
    % ----------------------------------------------------------------------------
    % output : ndarray of shape (2**N - 1,)
    %     the binary MLS values
    %
    shift_register = ones(1, N);
    values = [];
    for i = 1:2^N-1
            %%% Grab last index of register and amend to output sequence
            newValue = mod(sum(shift_register(output_taps)),2);
            values = [values, newValue];
            
            %%% Calculate a new bit for the register
            shiftval = mod(sum(shift_register(feedback_taps)),2);
            
            %%% Shift the register and insert new bit
            shift_register = [shiftval, shift_register(1:end-1)];
    end

end
