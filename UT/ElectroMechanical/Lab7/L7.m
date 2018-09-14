clc
close all

w = [5:100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Workspace Key
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cal0, cal3, cal2, cal1 (=) 0g, 70g, 150g, 250g
%
% acc1 = acceleration data set 1
% acc2 = acceleration data set 2
% 
% force1 = force data set 1
% force2 = force data set 2
%
% w = frequencies from 5 to 100 Hz



%%% Old Natural Frequencies in vertical direction
% Piezo: w = 9 Hz
% MEMS: w = 9.33, 61.67 Hz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weight =[0,.07,.15,.25];
cal = [cal0,cal3,cal2,cal1];
p = polyfit(cal,weight,1)

fit = @(x) (p(1)*x + p(2)) * 9.8; % Enter voltage (V), returns corresponding force (N)

%%% Calibration Constants
% m = p(1)
% b = p(2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(w,acc1)
ylabel('Tip Acceleration')
figure(2)
plot(w,acc2)
ylabel('Tip Acceleration')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Tip displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp1=zeros(96,1);
disp2=zeros(96,1);
for i=1:96; 
    disp1(i) = acc1(i)/((w(i)^2)* (cos(w(i)*5) + sqrt(-1)*sin(w(i)*5)));
    disp2(i) = acc2(i)/((w(i)^2)* (cos(w(i)*5) + sqrt(-1)*sin(w(i)*5)));
end

figure(3)
plot(w,disp1)
ylabel('Tip Displacement (m)');
figure(4)
plot(w,disp2)
ylabel('Tip Displacement (m)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F1=zeros(96,1);
F2=zeros(96,1);
for i=1:96; 
    F1(i) = fit(force1(i));
    F2(i) = fit(force2(i));
end
F1

figure(5)
plot(w,F1)
ylabel('Force (N)')

figure(6)
plot(w,F2)
ylabel('Force (N)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Tip Displacement / Force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DF1 = zeros(96,1);
DF2 = zeros(96,1);

for i=1:96;
    DF1(i) = abs(real(disp1(i)))/abs(F1(i));
    DF2(i) = abs(real(disp2(i)))/abs(F2(i));
end

figure(7)
plot(w,DF1)
ylabel('Displacement/Force (m/N)')

figure(8)
plot(w,DF2)
ylabel('Displacement/Force (m/N)')




% % % % fs = 1000;              % sampling frequency, Hz
% % % % T = 3;                  % duration of signal, s
% % % % 
% % % % dt = 1/fs;              % calculate time interval, s
% % % % tvec = [0:dt:T-dt]';    % generate time vector
% % % % 
% % % % %%%data = load('A1_data.txt','-ascii'); 
% % % % data = A1;
% % % % tipAccel = data(:,1)-mean(data(:,1));
% % % % tipVel = zeros(length(tvec),1);
% % % % tipDisp = zeros(length(tvec),1);
% % % % 
% % % % % % integrate to find velocity and displacement
% % % % % % for ii = 2:length(tvec),
% % % % % %     tipVel(ii) = tipVel(ii-1) + dt/2*(tipAccel(ii-1)+tipAccel(ii));
% % % % % %     tipDisp(ii) = tipDisp(ii-1) + dt/2*(tipVel(ii-1)+tipVel(ii));
% % % % % % end
% % % % tipVel = cumtrapz(tipAccel)*dt;
% % % % tipDisp = cumtrapz(tipVel)*dt;
% % % % 
% % % % figure(1)
% % % % subplot(311);
% % % % plot(tvec,tipAccel,'k-');
% % % % ylabel('Tip Acceleration, m^2/s');
% % % % grid on
% % % % subplot(312);
% % % % plot(tvec,tipVel,'k-');
% % % % ylabel('Tip Velocity, m/s');
% % % % grid on
% % % % subplot(313);
% % % % plot(tvec,tipDisp,'k-');
% % % % ylabel('Tip Displacement, m');
% % % % xlabel('Time, s');
% % % % grid on
% % % % 
% % % % dftAccel = fft(tipAccel);                               % calculate FFT
% % % % dftAccel = dftAccel(1:length(tvec)/2+1);                % throw away negative frequencies
% % % % %magAccel = 1/length(tvec)*abs(dftAccel);               % only magnitude at each frequency
% % % % psdAccel = (1/(fs*length(tvec))) * abs(dftAccel).^2;    % calculate power and scale it
% % % % psdAccel(2:end-1) = 2*psdAccel(2:end-1);                % multiply by two to account for negative frequencies
% % % % fvec = 0:fs/length(tvec):fs/2;                          % generate properly scaled frequency vector
% % % % 
% % % % figure(2)
% % % % subplot(211)
% % % % plot(tvec,tipAccel,'k-');
% % % % ylabel('Tip Acceleration, m^2/s');
% % % % xlabel('Time, s');
% % % % grid on
% % % % subplot(212);
% % % % plot(fvec,10*log10(psdAccel))
% % % % %plot(fvec,magAccel, 'b-');
% % % % grid on
% % % % xlabel('Frequency, Hz')
% % % % ylabel('Power/Frequency (dB/Hz)')






