clear
clc
close all
load exam2_jerkcar.mat
% dt, P0, uvec, xhat0, xhistgroundtrack, yacchist, yposhist, Xdes
% ------------------------------------------------------------------------
%%% 4b
% ------------------------------------------------------------------------
A1 = [0 1 0; 0 0 1; 0 0 0];
B1 = [0; 0; 1];
C1 = [0 0 0; 0 0 0];
D1 = [0; 0];
Gamma1 = [0; 0; 1];
Wj = 5e-4; % (m/s^3)^2

sys1 = ss(A1,B1,C1,D1);
sysd = c2d(sys1,dt,'ZOH');
G1 = sysd.B

M1 = [-A1, Gamma1*Wj*Gamma1'; zeros(3,3), A1']*dt

MM1 = expm(M1)

F1 = MM1(4:6,4:6)'

Q1 = F1 * MM1(1:3,4:6)

% ------------------------------------------------------------------------
%%% 4c
% ------------------------------------------------------------------------

Tau = 0.05; % Bias time constant
Wb = 0.053;
A = [0 1 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 Tau];
B = [0; 0; 1; 0];
C = [1 0 0 0; 0 0 1 0];
D = zeros(2,1);

Gamma = [0 0; 0 0; 1 0; 0 1];
W = [Wj 0; 0 Wb];

sys = ss(A,B,C,D);
sysd = c2d(sys,dt,'ZOH');
G = sysd.B

M = [-A, Gamma*W*Gamma'; zeros(4,4), A']*dt

MM = expm(M)

F = MM(5:8,5:8)'

Q = F * MM(1:4,5:8)


% ------------------------------------------------------------------------
%%% 4f
% ------------------------------------------------------------------------
tvec = dt:dt:20;
Ra = 0.05; % (m/s^2)^2
Rp = 0.5; % m^2

muaHistfilt = zeros(4,length(tvec));
    muaHistfilt(:,1) = xhat0;
PaHistfilt = zeros(4,4,length(tvec));
    PaHistfilt(:,:,1) = P0;

for kk = 1:length(tvec)
    if isnan(yposhist(kk))
        H = [0 0 1 1];
        R = Ra;
    else
        H = [1 0 0 0; 0 0 1 1];
        R = [Rp 0; 0 Ra];
    end

    %%% Prediction Step
    mua_kkp1_minus = F*muaHistfilt(:,kk) + G*uvec(kk); % 28x1
    Pa_kkp1_minus = F*PaHistfilt(:,:,kk)*F' + Q; % 28x28
    
    %%% Kalman Gain
    Kkkp1 = Pa_kkp1_minus * H' / (H*Pa_kkp1_minus*H' + R);

    %%% Measurement Update
    if isnan(yposhist(kk)) % If no GPS data
        ykkp1 = yacchist(kk);
    else
        ykkp1 = [yposhist(kk); yacchist(kk)]; % If GPS data
    end
    
    mua_kkp1_plus = mua_kkp1_minus + Kkkp1*(ykkp1 - H*mua_kkp1_minus);
    Pa_kkp1_plus = (eye(4) - Kkkp1*H)*Pa_kkp1_minus;
    muaHistfilt(:,kk+1) = mua_kkp1_plus;
    PaHistfilt(:,:,kk+1) = Pa_kkp1_plus;

end

%%% Plotting
% figure
% hold all
% plot(muaHistfilt(1,:),'r.')
% plot(xhistgroundtruth(1,:),'.')
% plot(yposhist,'c.')
% title('Position')
% legend('Did','Meant to do')
% 
% figure
% hold all
% plot(muaHistfilt(2,:),'r.')
% plot(xhistgroundtruth(2,:),'.')
% title('Velocity')
% legend('Did','Meant to do')
% 
% figure
% hold all
% plot(muaHistfilt(3,:),'r.')
% plot(xhistgroundtruth(3,:),'.')
% plot(yacchist-xhistgroundtruth(4,:),'c.')
% title('Acceleration')
% legend('Did','Meant to do')
% 
% figure
% hold all
% plot(muaHistfilt(4,:),'r.')
% plot(xhistgroundtruth(4,:),'.')
% title('Bias')
% legend('Did','Meant to do')

error_Hist = muaHistfilt(:,2:2001)-xhistgroundtruth(:,1:2000);
figure

for ss=1:4
    
    eval(['subplot(41',num2str(ss),')'])
    plot(error_Hist(ss,:),'b'),hold on
    plot(2*sqrt(squeeze(PaHistfilt(ss,ss,:)))','b--'),
    plot(-2*sqrt(squeeze(PaHistfilt(ss,ss,:)))','b--')
    legend('Estimated State Error','2-Sigma Error Bounds')
    
    if ss == 1
        PlotBoi2('','Position Error [m]',16)
    elseif ss == 2
        PlotBoi2('','Velocity Error [m/s]',16)
    elseif ss == 3
        PlotBoi2('','Acceleration Error [m/s^2]',16)
    elseif ss == 4
        PlotBoi2('Time Steps','Accelerometer Bias Error',16)
    end
    ylim([min(-2*sqrt(squeeze(PaHistfilt(ss,ss,:)))) max(2*sqrt(squeeze(PaHistfilt(ss,ss,:))))])
    xlim([0 2000]);
end







