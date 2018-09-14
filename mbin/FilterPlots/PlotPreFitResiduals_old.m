%%% Inputs
% 1) Vector of Times (nT x 1)
% 2) Matrix of Pre-Fit Residuals (nM x nT)
% 4) Error associated w/ Measurement Type 1 (sigma) (1 x 1)
% 5) Error associated w/ Measurement Type 2 (sigma) (1 x 1)
% 6) Desired Font Size

function PlotPreFitResiduals_old(time, yks, sig1, sig2, fs)
 
if size(yks,1) ~= 2
    warning('Function assumes 2 measurement types')
end
yks(yks == 0) = NaN;
time = time./3600;
figure

subplot(2,1,1); hold all; xlim([time(1) time(end)]);
plot(time, yks(1,:),'mx')
plot([time(1) time(end)],[3*sig1 3*sig1],'--r')
plot([time(1) time(end)],[-3*sig1 -3*sig1],'--r')
title('Prefit Residuals')
PlotBoi2('','Range Residuals, km',fs)

subplot(2,1,2); hold all; xlim([time(1) time(end)]);
plot(time, yks(2,:),'mx')
plot([time(1) time(end)],[3*sig2 3*sig2],'--r')
plot([time(1) time(end)],[-3*sig2 -3*sig2],'--r')
PlotBoi2('Times, hours','Range Rate Residuals, km/s',fs)
    
end