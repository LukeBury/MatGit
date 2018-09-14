close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Amplifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Averaging amp data to shorten vertical spectrum. n = how far out the
%%% avg will go
AmpTest = zeros(8400,1);
n=25;
for i=1:n:8400
    av=0;
    for r = 0:(n-1)
        av = av + Amp1(i+r);
    end
    avg = av/n;
    
    for w = 0:(n-1)
    AmpTest(i+w) = avg;
    end
    
end


%%% Other method to average amp data, not as good (probably) 
% for i=1:8400
%     sum1=0;
%     sum2=0;
%     %av = (Amp1(i)+Amp1(i+1)+Amp1(i+2)+Amp1(i+3)+Amp1(i+4)+Amp1(i+5)+Amp1(i+6)+Amp1(i+7)+Amp1(i+8)+Amp1(i+9))/10;
%     for k = 0:3
%         if (i+k) <= 8400
%         sum1 = sum1 + Amp1(i+k);
%         end
%         
%         if (i-k) >= 1
%         sum2 = sum2 + Amp1(i-k);
%         end
%     end
%     av = (sum1+sum2+Amp1(i))/9;
%     AmpTest(i)   = av;
%   
% end

%%% Amp1 = raw amp data
figure(1)
plot(Amp1)
xlabel('Measurements','FontName','Cambria','Fontsize',18)
ylabel('Temperature','FontName','Cambria','Fontsize',18)
set(gcf,'color','white')

hx = graph2d.constantline(725, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(2395, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(3475, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(4210, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(5010, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(6010, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(7295, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(8330, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

%%% Shortened amp data with verical lines representing breaks in songs
figure(2)
plot(AmpTest)
hx = graph2d.constantline(725, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(2395, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(3475, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(4210, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(5010, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(6010, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(7295, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(8330, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

xlabel('Measurements','FontName','Cambria','Fontsize',18)
ylabel('Temperature','FontName','Cambria','Fontsize',18)
set(gcf,'color','white')
%%% Best-fit to amp data, probably not useful
%  plot(Amp1)
%  
%  figure
% 
%  plot(AmpTest)
%  hold all
%  
%  x=[1:8400]';
%  poly = polyfit(x,AmpTest,7);
%  y = poly(1).*(x.^7) + poly(2).*(x.^6) + poly(3).*(x.^5) + poly(4).*(x.^4) + poly(5).*(x.^3) + poly(6).*(x.^2) + poly(7).*(x.^1) + poly(8);
%  
%  plot(x,y,'linewidth',2)

 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Longnew = Long shifted over to account for extra songs
%%% LongTest = Long w/ averages
LongTest = zeros(8400,1);
for k=1:30:8400
    
    sum=0;
    for a = 0:29
        sum = sum + Longnew(k+a);
    end
    avg = sum/30;
    
    for b = 0:29
        LongTest(k+b) = avg;
    end
    
   
end
figure(3)
plot(Longnew)
hx = graph2d.constantline(725, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(2395, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(3475, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(4210, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(5010, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(6010, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(7295, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(8330, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

xlabel('Measurements','FontName','Cambria','Fontsize',18)
ylabel('Temperature','FontName','Cambria','Fontsize',18)
set(gcf,'color','white')

figure(4)
hold all
plot(LongTest)


%%% Fit to Longnew
q=[1:8400]';
pol = polyfit(x,Longnew,9);
w = pol(1).*(x.^9)+pol(2).*(x.^8)+pol(3).*(x.^7)+pol(4).*(x.^6)+pol(5).*(x.^5)+pol(6).*(x.^4)+pol(7).*(x.^3)+pol(8).*(x.^2)+pol(9).*x+pol(10);
  
plot(q,w,'linewidth',2)

xlabel('Measurements','FontName','Cambria','Fontsize',18)
ylabel('Temperature','FontName','Cambria','Fontsize',18)
set(gcf,'color','white')

hx = graph2d.constantline(725, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(2395, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(3475, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(4210, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(5010, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(6010, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(7295, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

hx = graph2d.constantline(8330, 'LineStyle','-', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Use Me - Bill Withers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;              % sampling frequency, Hz
T = 10;                  % duration of signal, s
dt = 1/fs;              % calculate time interval, s
tvec = [0:dt:T-dt]';    % generate time vector
L=441000;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y

UseTT = UseT - mean(UseT);

figure(5)
subplot(211)
Y = fft(UseP,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2+1));
YYY = YY./max(YY);
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f,YYY) 
title('Bill Withers - Use Me','FontName','Cambria','Fontsize',14)
xlabel('Frequency (Hz)','FontName','Cambria','Fontsize',16)
ylabel('Normalized Pressure','FontName','Cambria','Fontsize',16)
axis([0 800 0 1])
subplot(212)
Y = fft(UseTT,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2+1));
YYY = YY./max(YY);
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f,YYY) 
xlabel('Frequency (Hz)','FontName','Cambria','Fontsize',16)
ylabel('Normalized Temperature','FontName','Cambria','Fontsize',16)
set(gcf,'color','white')

A = 2*abs(Y(1:NFFT/2+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Is There Any Love - Kid Cudi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LoveTT = LoveT - mean(LoveT);

figure(6)
subplot(211)
Y = fft(LoveP,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2+1));
YYY = YY./max(YY);
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f,YYY) 
title('Is There Any Love - Kid Cudi','FontName','Cambria','Fontsize',14)
xlabel('Frequency (Hz)','FontName','Cambria','Fontsize',16)
ylabel('Normalized Pressure','FontName','Cambria','Fontsize',16)
axis([0 800 0 1])
subplot(212)
Y = fft(LoveTT,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2+1));
YYY = YY./max(YY);
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f,YYY) 
xlabel('Frequency (Hz)','FontName','Cambria','Fontsize',16)
ylabel('Normalized Temperature','FontName','Cambria','Fontsize',16)
set(gcf,'color','white')

B = 2*abs(Y(1:NFFT/2+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% I'm a Wheel - Wilco
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WheelTT = WheelT - mean(WheelT);

figure(7)
subplot(211)
Y = fft(WheelP,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2+1));
YYY = YY./max(YY);
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f,YYY) 
title('I''m a Wheel - Wilco','FontName','Cambria','Fontsize',14)
xlabel('Frequency (Hz)','FontName','Cambria','Fontsize',16)
ylabel('Normalized Pressure','FontName','Cambria','Fontsize',16)
axis([0 800 0 1])
subplot(212)
Y = fft(WheelTT,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2+1));
YYY = YY./max(YY);
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f,YYY) 
xlabel('Frequency (Hz)','FontName','Cambria','Fontsize',16)
ylabel('Normalized Temperature','FontName','Cambria','Fontsize',16)
set(gcf,'color','white')

C = 2*abs(Y(1:NFFT/2+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% The World is Yours - Nas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nas1TT = Nas1T - mean(Nas1T);

figure(8)
subplot(211)
Y = fft(Nas1P,NFFT)/L;
YY = 2*abs(Y(1:NFFT/2+1));
YYY = YY./max(YY);
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f,YYY) 
title('The World is Yours - Nas','FontName','Cambria','Fontsize',14)
xlabel('Frequency (Hz)','FontName','Cambria','Fontsize',16)
ylabel('Normalized Pressure','FontName','Cambria','Fontsize',16)
%axis([0 800 0 4e-3])
axis([0 800 0 1])
subplot(212)
Y = fft(Nas1TT,NFFT)/L;
YY= 2*abs(Y(1:NFFT/2+1));
YYY = YY./max(YY);
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f,YYY) 
xlabel('Frequency (Hz)','FontName','Cambria','Fontsize',16)
ylabel('Normalized Temperature','FontName','Cambria','Fontsize',16)
set(gcf,'color','white')

D = 2*abs(Y(1:NFFT/2+1));







