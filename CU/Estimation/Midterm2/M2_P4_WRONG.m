clear
clc
close all
load exam2_jerkcar.mat
% dt, P0, uvec, xhat0, xhistgroundtrack, yacchist, yposhist
% ------------------------------------------------------------------------
%%% 4b
% ------------------------------------------------------------------------
A1 = [0 1 0; 0 0 1; 0 0 0];
B1 = [0; 0; 1];
Gamma1 = [0; 0; 1];
W1 = 5e-4; % (m/s^3)^2

M1 = [-A1, Gamma1*W1*Gamma1'; zeros(3,3), A1']*dt

MM1 = expm(M1)

F1 = MM1(4:6,4:6)'

Q1 = F1 * MM1(1:3,4:6)

% ------------------------------------------------------------------------
%%% 4c
% ------------------------------------------------------------------------
A = [0 1 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
C = [1 0 0 0; 0 0 1 0];
D = zeros(2,2);
Gamma = [0 0; 0 0; 1 0; 0 1];
W = [5e-4 0; 0 0.053];

sys = ss(A,B,C,D);
sysd = c2d(sys,dt,'ZOH');

% F3 = sysd.A
G = sysd.B
% H3 = sysd.C;
% M3 = sysd.D


M = [-A, Gamma*W*Gamma'; zeros(4,4), A']*dt

MM = expm(M)

F = MM(5:8,5:8)'

Q = F * MM(1:4,5:8)

% G = [0 0; 0 0; 1 0; 0 1];

% ------------------------------------------------------------------------
%%% 4f
% ------------------------------------------------------------------------


tvec = dt:dt:20;
mu0a = xhat0;
P0a = P0;
Ra = 0.05; % (m/s^2)^2
Rp = 0.5; % m^2
Tau = 0.05; % Bias time constant

muaHistfilt = zeros(4,length(tvec));
    muaHistfilt(:,1) = mu0a;
PaHistfilt = zeros(4,4,length(tvec));
    PaHistfilt(:,:,1) = P0a;

for kk = 1:length(tvec)
    if isnan(yposhist(kk))
        H = [0 0 1 0];
        R = Ra;
    else
        H = [1 0 0 0; 0 0 1 0];
        R = [Rp 0; 0 Ra];
    end
    u = [uvec(kk); Tau]; % Create the new control vector

    %%% Prediction Step
    mua_kkp1_minus = F*muaHistfilt(:,kk) + G*u; % 28x1
    Pa_kkp1_minus = F*PaHistfilt(:,:,kk)*F' + Q; % 28x28
    
    %%% Kalman Gain
    Kkkp1 = Pa_kkp1_minus * H' / (H*Pa_kkp1_minus*H' + R);

    %%% Measurement Update
    if isnan(yposhist(kk))
        ykkp1 = yacchist(kk);
    else
        ykkp1 = [yposhist(kk); yacchist(kk)];
    end
    mua_kkp1_plus = mua_kkp1_minus + Kkkp1*(ykkp1 - H*mua_kkp1_minus);
    Pa_kkp1_plus = (eye(4) - Kkkp1*H)*Pa_kkp1_minus;
    muaHistfilt(:,kk+1) = mua_kkp1_plus;
    PaHistfilt(:,:,kk+1) = Pa_kkp1_plus;
%     
%     
%     if mod(kk,10)==0 || kk == 1
%         [Xell,Yell] = calc_gsigma_ellipse_plotpoints...
%         ([muaHistfilt(1,kk),muaHistfilt(3,kk)]',...
%          [PaHistfilt(1,1,kk),PaHistfilt(1,3,kk);PaHistfilt(3,1,kk),PaHistfilt(3,3,kk)],...
%           2,100);
%         plot(Xell,Yell,'b.-') %%plot new ellipse
%         
%     end
end
hold all
plot(muaHistfilt(1,:),'r.')
plot(xhistgroundtruth(1,:),'.')
plot(yposhist,'c.')
title('Position')
legend('Did','Meant to do')

figure
hold all
plot(muaHistfilt(2,:),'r.')
plot(xhistgroundtruth(2,:),'.')
title('Velocity')
legend('Did','Meant to do')

figure
hold all
plot(muaHistfilt(3,:),'r.')
plot(xhistgroundtruth(3,:),'.')
plot(yacchist,'c.')
title('Acceleration')
legend('Did','Meant to do')

figure
hold all
plot(muaHistfilt(4,:),'r.')
plot(xhistgroundtruth(4,:),'.')
title('Bias')
legend('Did','Meant to do')
% %Plotting Measurements
% for i = 1:12
%     p3 = plot(Yfulllin(2*i-1,:)+muaHistfilt(1,:), Yfulllin(2*i,:)+muaHistfilt(3,:),'c.');
% end
% p1 = plot(muaHistfilt(1,:),muaHistfilt(3,:),'linewidth',2);
% p2 = plot(xRTrueHist(1,:),xRTrueHist(3,:),'linewidth',2);
% p4 = plot(xLtrue(1,:),xLtrue(2,:),'kx','linewidth',2,'markersize',8);
% legend([p1 p2 p4 p3],'KF','True','Landmarks','Measurements')
% PlotBoi2('East Position, [m]','North Position, [m]',18)
% 
%plot state errors
error_Hist = muaHistfilt(:,2:2000)-xhistgroundtruth(:,1:1999);
figure

for ss=1:4
    
    eval(['subplot(41',num2str(ss),')'])
    plot(error_Hist(ss,:),'b'),hold on
    plot(2*sqrt(squeeze(PaHistfilt(ss,ss,:)))','b--'),
    plot(-2*sqrt(squeeze(PaHistfilt(ss,ss,:)))','b--')
    legend('Estimated State Error','2-Sigma Error Bounds')
    
    if ss == 1
        PlotBoi2('','Position Error',9)
    elseif ss == 2
        PlotBoi2('','Velocity Error',9)
    elseif ss == 3
        PlotBoi2('','Acceleration Error',9)
    elseif ss == 4
        PlotBoi2('Time Steps','Accelerometer Bias Error',9)
    end
end
