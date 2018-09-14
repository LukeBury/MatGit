%%% Inputs
% 1) Vector of Times (nT x 1)
% 2) Matrix of Pre-Fit Residuals (nM x nT)
% 3) Vector of station number at each time [nData x 2]
% 4) Error associated w/ Measurement Type 1 (sigma) (1 x 1)
% 5) Error associated w/ Measurement Type 2 (sigma) (1 x 1)
% 6) Desired Font Size

function PlotPreFitResiduals(time, yks, statNums, sig1, sig2, fs)
 
if size(yks,1) ~= 2
    warning('Function assumes 2 measurement types')
end
%%% Parse yks by stations
statPrs = [statNums(1) 1;statNums(1) 0];
row = 3;
for k = 1:length(statNums)
    if statNums(k) ~= statPrs(row-1,1)
        kt = find(time == statNums(k-1,2));
        statPrs(row-1, 2) = kt;
        statPrs(row, 1) = statNums(k);
        row = row+1;
    end
end
statPrs(end,2) = find(time == statNums(end,2));

yks(yks == 0) = NaN;
time = time./3600;
s1_On = 0; s2_On = 0; s3_On = 0;
figure
subplot(2,1,1); hold all; xlim([time(1) time(end)]);
for k = 1:size(statPrs,1)-1
    t1 = statPrs(k,2)+1;
    t2 = statPrs(k+1,2);
    if statPrs(k+1,1) == 1
        p1 = plot(time(t1:t2),yks(1,t1:t2),'mx');
        s1_On = 1;
    elseif statPrs(k+1,1) == 2
        p2 = plot(time(t1:t2),yks(1,t1:t2),'rx');
        s2_On = 1;
    elseif statPrs(k+1,1) == 3
        p3 = plot(time(t1:t2),yks(1,t1:t2),'bx');
        s3_On = 1;
    end
end

plot([time(1) time(end)],[3*sig1 3*sig1],'--r')
plot([time(1) time(end)],[-3*sig1 -3*sig1],'--r')
%%% Legend switch
if     s1_On == 1 && s2_On == 0 && s3_On == 0
    legend([p1],'Canberra')
elseif s1_On == 0 && s2_On == 1 && s3_On == 0
    legend([p2], 'Madrid')
elseif s1_On == 0 && s2_On == 0 && s3_On == 1
    legend([p3], 'Goldstone')
elseif s1_On == 1 && s2_On == 1 && s3_On == 0
    legend([p1, p2], 'Canberra','Madrid')
elseif s1_On == 1 && s2_On == 0 && s3_On == 1
    legend([p1, p3], 'Canberra','Goldstone')
elseif s1_On == 0 && s2_On == 1 && s3_On == 1
    legend([p2, p3],'Madrid','Goldstone')
elseif s1_On == 1 && s2_On == 1 && s3_On == 1
    legend([p1, p2, p3],'Canberra','Madrid','Goldstone')
end
title('Prefit Residuals')
PlotBoi2('','Range Residuals, km',fs)

subplot(2,1,2); hold all; xlim([time(1) time(end)]);
for k = 1:size(statPrs,1)-1
    t1 = statPrs(k,2)+1;
    t2 = statPrs(k+1,2);
    if statPrs(k+1,1) == 1
        p1 = plot(time(t1:t2),yks(2,t1:t2),'mx');
    elseif statPrs(k+1,1) == 2
        p2 = plot(time(t1:t2),yks(2,t1:t2),'rx');
    elseif statPrs(k+1,1) == 3
        p3 = plot(time(t1:t2),yks(2,t1:t2),'bx');
    end
end
plot([time(1) time(end)],[3*sig2 3*sig2],'--r')
plot([time(1) time(end)],[-3*sig2 -3*sig2],'--r')
PlotBoi2('Times, hours','Range Rate Residuals, km/s',fs)
    
end