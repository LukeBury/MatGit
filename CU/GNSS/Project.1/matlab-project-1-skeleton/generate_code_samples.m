function [code_samples] = generate_code_samples(fs, L, code, fc, varargin)
    % Generates samples of code sequence given sampling duration / rate and
    % code sequence, code chipping rate, and initial code phase (optional).
    %
    % Parameters
    % ----------------------------------------------------------------------------
    % fs : float 
    %     sampling rate
    % L : float
    %     duration of sampled signal (seconds)
    % code : array of shape (M,)
    %     code sequence to sample
    % fc : float
    %     code chipping rate
    % c0 (varargin) : float
    %     (optional) defaults to zero -- the initial code phase in the sampled
    %     time-series
    %  
    % Returns
    % ----------------------------------------------------------------------------
    % output : ndarray of shape (N,)
    %     the code samples
    c0 = 0;
    if nargin > 4
      c0 = varargin(1);
    end
    
    % Total number of sampled instances -> length of code_samples
    nSamples = L*fs;
    
    % Creating time vector
    time = linspace(0,L,nSamples); % sec
    
    % Initializing code_samples array
    code_samples = zeros(1,length(time));

    for kk = 1:length(time)
        % Determining chip value at time(kk) .. not bounded to index 1024
        nChip = floor(time(kk)*fc) + c0;
        
        % Resetting chip value range from 1-1024
        multiple = floor(nChip/1023);
        nChip = nChip - multiple*1023 + 1;
        
        % Storing code sample
        code_samples(kk) = code(nChip);
            
    end    

end
