%%% Inputs
% 1) Vector of Times (nT x 1)
% 2) Matrix of State Errors (nT x nS)
% 3) Matrix of Covariances (nS x nS x nT)
% 4) Desired Font Size

function PlotEstimateError(time, error, covariance, fs)
if rem(size(error,2),2) ~= 0
    warning('Function assumes even number of states')
    return
end
nS = size(error,2);
figure
s = 1;
for sP = 1:2:nS-1 % 1, 3, 5
    subplot(nS/2,2,sP); hold all; xlim([time(1) time(end)]);
    plot(time, error(:,s))
    plot(time,3*sqrt(squeeze(covariance(s,s,:))),'--r')
    plot(time,-3*sqrt(squeeze(covariance(s,s,:))),'--r')
    if s == 1
        PlotBoi2('','X Position Error, km',fs)
    elseif s == 2
        PlotBoi2('','Y Position Error, km',fs)
    elseif s == 3
        PlotBoi2('Time, sec','Z Position Error, km',fs)
    end
    s = s+1;
end

for sP = 2:2:nS % 2, 4, 6
    subplot(nS/2,2,sP); hold all; xlim([time(1) time(end)]);
    plot(time, error(:,s))
    plot(time,3*sqrt(squeeze(covariance(s,s,:))),'--r')
    plot(time,-3*sqrt(squeeze(covariance(s,s,:))),'--r')
    if s == 4
        PlotBoi2('','X Velocity Error, km/s',fs)
        legend('Truth','3\sigma Covariance')
    elseif s == 5
        PlotBoi2('','Y Velocity Error, km/s',fs)
    elseif s == 6
        PlotBoi2('Time, sec','Z Velocity Error, km/s',fs)
    end
    s = s+1;
end

end