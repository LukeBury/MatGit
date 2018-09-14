clear
clc
close all
% ASEN 5090 Assignment 1
 
% 1. Complete the function `generate_mls` that generates maximum-length binary sequences.
% 2. Use the function `generate_mls` from part 1 to complete the function `generate_l1ca_code`, which generates and returns the binary the GPS L1 C/A code sequence given the GPS PRN number.
% 3. Complete the function `generate_code_samples` that generates a time-series of code samples given the sampling rate, sampling duration, code sequence, code chipping rate, and initial code phase.
% 4. Choose 2 PRNs, sample their code sequences, and convert the samples to +-1 (i.e. a modulated signal of constant magnitude).  Plot the auto-correlation for each and the cross-correlation of the two signals.

%  Verify that the C/A code sequences match their octal representation for the
%  first 10 bits.
passed = [];
failed = [];
for prn = 1:32
  res = verify_l1ca_code_first_10_chips(prn);
  if res
    passed = [passed, prn];
  else
    failed = [failed, prn];
  end
end
if length(passed) == 32
    disp('All PRNs passed');
elseif ~isempty(failures)
    fprintf('PRNs %s did NOT pass', mat2str(failed));
end


%%
prn_1 = 20;
prn_2 = 21;

% Setting Parameters
fs = 5e6;     % 5 MHz sampling rate
L = 3e-3;     % 5 milliseconds sampling duration
fc = 1.023e6; % GPS L1 C/A code chipping rate (Hz)

%%% Generating code_1 and code_1 sample
code_1 = generate_l1ca_code(prn_1);
samples_1 = generate_code_samples(fs, L, code_1, fc);

% Converting 1 -> -1 and 0 -> 1
samples_1(samples_1==1) = -1;
samples_1(samples_1==0) = 1;

%%% Generating code_2 and code_2 sample
code_2 = generate_l1ca_code(prn_2);
samples_2 = generate_code_samples(fs, L, code_2, fc);

% Converting 1 -> -1 and 0 -> 1
samples_2(samples_2==1) = -1;
samples_2(samples_2==0) = 1;

%%% Computing autocorrelations and cross-correlations between code samples
auto_1 = circular_correlation(samples_1, samples_1);
auto_2 = circular_correlation(samples_2, samples_2);
cross = circular_correlation(samples_1, samples_2);

% Plotting autocorrelations and cross-correlations
figure
subplot(3,1,1)
plot(auto_1./max(auto_1))
PlotBoi2('','PRN20 Code',18)
subplot(3,1,2)
plot(auto_2./max(auto_2))
PlotBoi2('','PRN21 Code',18)
subplot(3,1,3)
plot(cross)
PlotBoi2('Delay Sample Index','Cross-Correlation',18)

