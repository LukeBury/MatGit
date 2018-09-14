% LABORATORY 6
% code to process the measured acceleration signal
% plots time history 
% integrates to find displacement
% calculates power spectrum

% sirohi 102214

%%%clear all
close all

fs = 1000;              % sampling frequency, Hz
T = 3;                  % duration of signal, s

dt = 1/fs;              % calculate time interval, s
tvec = [0:dt:T-dt]';    % generate time vector

%%%data = load('A1_data.txt','-ascii'); 
a1 = [A1(145:3000);A1(2857:3000)];
a2 = [A2(320:3000);A2(2682:3000)];
a3 = [A3(380:3000);A3(2622:3000)];
aavg = (a1+a2+a3)./3; %1
%data = a1; %1
%data = a2; %1
%data = a3; %1
%data = aavg; %1


b1_1 = [B1_1(200:3000);B1_1(2802:3000)];
b2_1 = [B2_1(140:3000);B2_1(2862:3000)];
b3_1 = [B3_1(150:3000);B3_1(2852:3000)];
bavg_1 = (b1_1+b2_1+b3_1)./3; 

b1_2 = [B1_2(195:3000);B1_2(2807:3000)];
b2_2 = [B2_2(140:3000);B2_2(2862:3000)];
b3_2 = [B3_2(145:3000);B3_2(2857:3000)];
bavg_2 = (b1_2+b2_2+b3_2)./3; 

%data = b1_1; %2
%data = b1_2; %3
%data = b2_1; %2
%data = b2_2; %3
%data = b3_1; %2
%data = b3_2; %3

%data = aavg;   %1
data = bavg_1; %2
%data = bavg_2; %3

dataset = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dataset == 1
tipAccel = data(:,1)-mean(data(:,1));
tipVel = zeros(length(tvec),1);
tipDisp = zeros(length(tvec),1);

% % integrate to find velocity and displacement
% % for ii = 2:length(tvec),
% %     tipVel(ii) = tipVel(ii-1) + dt/2*(tipAccel(ii-1)+tipAccel(ii));
% %     tipDisp(ii) = tipDisp(ii-1) + dt/2*(tipVel(ii-1)+tipVel(ii));
% % end
tipVel = cumtrapz(tipAccel)*dt;
tipDisp = cumtrapz(tipVel)*dt;

figure(1)
subplot(311);
plot(tvec,tipAccel,'k-');
ylabel('Tip Acceleration, m^2/s','FontName','Cambria','Fontsize',18);
grid on
subplot(312);
plot(tvec,tipVel,'k-');
ylabel('Tip Velocity, m/s','FontName','Cambria','Fontsize',18);
grid on
subplot(313);
plot(tvec,tipDisp,'k-');
ylabel('Cumulative Tip Displacement, m','FontName','Cambria','Fontsize',18);
xlabel('Time, s','FontName','Cambria','Fontsize',18);
grid on

dftAccel = fft(tipAccel);                               % calculate FFT

dftAccel = dftAccel(1:length(tvec)/2+1);                % throw away negative frequencies
%magAccel = 1/length(tvec)*abs(dftAccel);               % only magnitude at each frequency
psdAccel = (1/(fs*length(tvec))) * abs(dftAccel).^2;    % calculate power and scale it
psdAccel(2:end-1) = 2*psdAccel(2:end-1);                % multiply by two to account for negative frequencies
fvec = 0:fs/length(tvec):fs/2;                          % generate properly scaled frequency vector

% figure(4)
% plot(fvec,real(dftAccel))
% freal = real(dftAccel);
% for i=0:size(freal)
%     if (freal(i) >= 5)
%         i
%         dftAccel(i)
%         freal(i)
%         
%     end
% end


figure(2)
subplot(211)
plot(tvec,tipAccel,'k-');
% ylim([-.04,.04])
ylabel('Tip Acceleration, m^2/s','FontName','Cambria','Fontsize',18);
xlabel('Time, s','FontName','Cambria','Fontsize',18);
grid on
subplot(212);
plot(fvec,10*log10(psdAccel))
ylim([-120,-20])
%plot(fvec,magAccel, 'b-');
grid on
xlabel('Frequency, Hz','FontName','Cambria','Fontsize',18)
ylabel('Power/Frequency (dB/Hz)','FontName','Cambria','Fontsize',18)



%%% One Second, Bro
figure(3)
subplot(211)
tvec2 = tvec(2000:3000);
plot(tvec(2000:3000),tipAccel(2000:3000),'k-');

xlim([2,3])
ylim([-.004,.004])
ylabel('Tip Acceleration, m^2/s','FontName','Cambria','Fontsize',18);
xlabel('Last Second of Time, s','FontName','Cambria','Fontsize',18);
grid on
subplot(212);
fvec2 = 0:fs/length(tvec2):fs/2;
%plot(fvec(1000:1501),10*log10(psdAccel(1000:1501))) %old
plot(fvec2,10*log10(psdAccel(1001:1501)))

%xlim([330,500]) %old
xlim([0,500])
ylim([-110,-70])
%plot(fvec,magAccel, 'b-');
grid on
xlabel('Last Second Frequency, Hz','FontName','Cambria','Fontsize',18)
ylabel('Power/Frequency (dB/Hz)','FontName','Cambria','Fontsize',18)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else if dataset == 2 
        
tipAccel = data(:,1)-mean(data(:,1));
tipVel = zeros(length(tvec),1);
tipDisp = zeros(length(tvec),1);

% % integrate to find velocity and displacement
% % for ii = 2:length(tvec),
% %     tipVel(ii) = tipVel(ii-1) + dt/2*(tipAccel(ii-1)+tipAccel(ii));
% %     tipDisp(ii) = tipDisp(ii-1) + dt/2*(tipVel(ii-1)+tipVel(ii));
% % end
tipVel = cumtrapz(tipAccel)*dt;
tipDisp = cumtrapz(tipVel)*dt;

figure(1)
subplot(311);
plot(tvec,tipAccel,'k-');
ylabel('Tip Acceleration, m^2/s','FontName','Cambria','Fontsize',18);
grid on
subplot(312);
plot(tvec,tipVel,'k-');
ylabel('Tip Velocity, m/s','FontName','Cambria','Fontsize',18);
grid on
subplot(313);
plot(tvec,tipDisp,'k-');
ylabel('Cumulative Tip Displacement, m','FontName','Cambria','Fontsize',18);
xlabel('Time, s','FontName','Cambria','Fontsize',18);
grid on

dftAccel = fft(tipAccel);                               % calculate FFT
dftAccel = dftAccel(1:length(tvec)/2+1);                % throw away negative frequencies
%magAccel = 1/length(tvec)*abs(dftAccel);               % only magnitude at each frequency
psdAccel = (1/(fs*length(tvec))) * abs(dftAccel).^2;    % calculate power and scale it
psdAccel(2:end-1) = 2*psdAccel(2:end-1);                % multiply by two to account for negative frequencies
fvec = 0:fs/length(tvec):fs/2;                          % generate properly scaled frequency vector

figure(2)
subplot(211)
plot(tvec,tipAccel,'k-');
xlim([0,3])
ylim([-2.5,1.5])
ylabel('Tip Acceleration, m^2/s','FontName','Cambria','Fontsize',18);
xlabel('Time, s','FontName','Cambria','Fontsize',18);
grid on
subplot(212);
plot(fvec,10*log10(psdAccel))
ylim([-100,0])
%plot(fvec,magAccel, 'b-');
grid on
xlabel('Frequency, Hz','FontName','Cambria','Fontsize',18)
ylabel('Power/Frequency (dB/Hz)','FontName','Cambria','Fontsize',18)




%%%%%%%%%%%%% One Second, Bro
figure(3)
subplot(211)
plot(tvec(2000:3000),tipAccel(2000:3000),'k-');

xlim([2,3])
ylim([-.15,.15])
ylabel('Tip Acceleration, m^2/s','FontName','Cambria','Fontsize',18);
xlabel('Last Second of Time, s','FontName','Cambria','Fontsize',18);
grid on
subplot(212);

tvec2 = tvec(2000:3000);
fvec2 = 0:fs/length(tvec2):fs/2;
plot(fvec2,10*log10(psdAccel(1001:1501)))
%plot(fvec(1000:1501),10*log10(psdAccel(1000:1501)))

xlim([0,500])
ylim([-70,-30])
%plot(fvec,magAccel, 'b-');
grid on
xlabel('Last Second Frequency, Hz','FontName','Cambria','Fontsize',18)
ylabel('Power/Frequency (dB/Hz)','FontName','Cambria','Fontsize',18)

    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else if dataset == 3 
            
tipAccel = data(:,1)-mean(data(:,1));
tipVel = zeros(length(tvec),1);
tipDisp = zeros(length(tvec),1);

% % integrate to find velocity and displacement
% % for ii = 2:length(tvec),
% %     tipVel(ii) = tipVel(ii-1) + dt/2*(tipAccel(ii-1)+tipAccel(ii));
% %     tipDisp(ii) = tipDisp(ii-1) + dt/2*(tipVel(ii-1)+tipVel(ii));
% % end
tipVel = cumtrapz(tipAccel)*dt;
tipDisp = cumtrapz(tipVel)*dt;

figure(1)
subplot(311);
plot(tvec,tipAccel,'k-');
ylabel('Tip Acceleration, m^2/s','FontName','Cambria','Fontsize',18);
grid on
subplot(312);
plot(tvec,tipVel,'k-');
ylabel('Tip Velocity, m/s','FontName','Cambria','Fontsize',18);
grid on
subplot(313);
plot(tvec,tipDisp,'k-');
ylabel('Cumulative Tip Displacement, m','FontName','Cambria','Fontsize',18);
xlabel('Time, s','FontName','Cambria','Fontsize',18);
grid on

dftAccel = fft(tipAccel);                               % calculate FFT
dftAccel = dftAccel(1:length(tvec)/2+1);                % throw away negative frequencies
%magAccel = 1/length(tvec)*abs(dftAccel);               % only magnitude at each frequency
psdAccel = (1/(fs*length(tvec))) * abs(dftAccel).^2;    % calculate power and scale it
psdAccel(2:end-1) = 2*psdAccel(2:end-1);                % multiply by two to account for negative frequencies
fvec = 0:fs/length(tvec):fs/2;                          % generate properly scaled frequency vector

figure(2)
subplot(211)
plot(tvec,tipAccel,'k-');
xlim([0,3])
ylim([-2,2])
ylabel('Tip Acceleration, m^2/s','FontName','Cambria','Fontsize',18);
xlabel('Time, s','FontName','Cambria','Fontsize',18);
grid on
subplot(212);
plot(fvec,10*log10(psdAccel))
ylim([-100,0])
%plot(fvec,magAccel, 'b-');
grid on
xlabel('Frequency, Hz','FontName','Cambria','Fontsize',18)
ylabel('Power/Frequency (dB/Hz)','FontName','Cambria','Fontsize',18)



%%% One Second, Bro
figure(3)
subplot(211)
plot(tvec(2000:3000),tipAccel(2000:3000),'k-');

xlim([2,3])
ylim([-.05,.05])
ylabel('Tip Acceleration, m^2/s','FontName','Cambria','Fontsize',18);
xlabel('Last Second of Time, s','FontName','Cambria','Fontsize',18);
grid on
subplot(212);

tvec2 = tvec(2000:3000);
fvec2 = 0:fs/length(tvec2):fs/2;
plot(fvec2,10*log10(psdAccel(1001:1501)))
%plot(fvec(1000:1501),10*log10(psdAccel(1000:1501)))

xlim([0,500])
ylim([-70,-30])
%plot(fvec,magAccel, 'b-');
grid on
xlabel('Last Second Frequency, Hz','FontName','Cambria','Fontsize',18)
ylabel('Power/Frequency (dB/Hz)','FontName','Cambria','Fontsize',18)
    
    end
end

        
end
        
        
        
        