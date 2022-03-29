%hw6_problem4template.m

clc,clear
close all

%% 0. Nominal scenario for A/C a and A/C b executing coord turn
dT = 0.5; %sec
tvec = 0:dT:100;

load hw6problem4_data.mat 

Omega_a = 0.045; %rad/s
Aa = [0 1 0 0; 0 0 0 -Omega_a; 0 0 0 1;0 Omega_a 0 0];

Omega_b = -0.045; %rad/s
Ab = [0 1 0 0; 0 0 0 -Omega_b; 0 0 0 1;0 Omega_b 0 0];

%% a. find the Process noise matrix
Gamma_a = [0 0 
           1 0;
           0 0;
           0 1];
Gamma_b = [0 0 
           1 0;
           0 0;
           0 1];
Wa = 10*[2.0 0.05;   %%model directional wind disturbances 
         0.05 0.5];
Wb = 10*[2.0 0.05;
         0.05 0.5];
     
 Ma = [-Aa Gamma_a*Wa*Gamma_a'; zeros(4,4) Aa']*dT;
 Mb = [-Ab Gamma_b*Wb*Gamma_b'; zeros(4,4) Ab']*dT;
 Ba = expm(Ma);
 Bb = expm(Mb);
 Fa = Ba(5:8,5:8)'
 Fb = Bb(5:8,5:8)'
 
%%Find Qa and Qb 

Qa = Fa*Ba(1:4,5:8)
Qb = Fb*Bb(1:4,5:8)

%% b. Pure prediction for a/c A
mu0a = [0; 85*cos(pi/4); 0; -85*sin(pi/4)]; %%initial mean of A at k=0
P0a = 900 * diag([10, 2, 10, 2],0); %%initial covar of A at k=0

muaHistpred = zeros(4,length(tvec));
    muaHistpred(:,1) = mu0a;
PaHistpred = zeros(4,4,length(tvec));
    PaHistpred(:,:,1) = P0a;

%%Draw samples to predict through state dynamics (for both A and B)
nsamps = 50;
xa_samps = mvnrnd(mu0a',P0a,nsamps)';
xa_sampHist = zeros(4,nsamps,length(tvec));
    xa_sampHist(:,:,1) = xa_samps;
    
%%Intialize plot for A    
figure(), 
[Xell,Yell] = calc_gsigma_ellipse_plotpoints...
        ([mu0a(1),mu0a(3)]',[P0a(1,1),P0a(1,3);P0a(3,1),P0a(3,3);],2,100);
plot(Xell,Yell,'b.-'),hold on
plot(xa_samps(1,:),xa_samps(3,:),'r.')
PlotBoi2('East Position, [m]','North Position, [m]',18)
for kk=1:length(tvec)-1
    %%Pure prediction update for A
    muaHistpred(:,kk+1) = Fa^(kk)*mu0a;
    PaHistpred(:,:,kk+1) = (Fa^(kk)) * P0a * (Fa^(kk))' + Qa;
    
    %%Sample update
    %%simulate process noise vector samples
%     wkka_samp = mvnrnd(muaHistpred(:,kk)',PaHistpred(:,:,kk),50)';
%     wkka_samp = mvnrnd(zeros(1,4),PaHistpred(:,:,kk),50)';
    wkka_samp = mvnrnd(zeros(4,1),Qa,50)';
    xa_sampHist(:,:,kk+1) = Fa*(xa_sampHist(:,:,kk)) + wkka_samp;
    
    %%update plot for A
    if mod(kk,25)==0
        [Xell,Yell] = calc_gsigma_ellipse_plotpoints...
        ([muaHistpred(1,kk+1),muaHistpred(3,kk+1)]',...
         [PaHistpred(1,1,kk+1),PaHistpred(1,3,kk+1);PaHistpred(3,1,kk+1),PaHistpred(3,3,kk+1)],...
          2,100);
        plot(Xell,Yell,'b.-') %%plot new ellipse   
        %%plot new Monte Carlo sample points
        plot(xa_sampHist(1,:,kk+1),xa_sampHist(3,:,kk+1),'r.')
    end
    
end
axis square
fprintf('bii:\n')
fprintf('The shape (circular) indicates we have essentially\n equal uncertainty in radial and angular directions. \nThe circles grow because uncertainty compounds over \ntime. The points do no stay inside the ellipses, which\n is problematic\n')
fprintf('you''d expect 5/100 points outside ellipse, and\n we''re seeing much more\n\n')
fprintf('nvm dots are ok!\n')
%% c. Measurement and KFs
rng(100)
%%% xasingle_truth; %%--> get from data log

%%%i. simulate measurements of A from ground tracking station
Ra = [20 0.05; 0.05 20]; % m^2
Ha = [1 0 0 0; 0 0 1 0];
% xasingle_truth = ground truth trajectory A
yaHist = zeros(2,length(tvec)-1);
for kk=2:length(tvec)
    vakk_samp = mvnrnd(zeros(2,1),Ra,1)'; %%simulate random ground station measurement noise
    yaHist(:,kk-1) = Ha*(xasingle_truth(:,kk)) + vakk_samp; %%simulate noisy ground station measurement
end

%%plot raw measurements
figure(),
plot(tvec(1:41),yaHist(1,1:41),'rx'), hold on
plot(tvec(1:41),yaHist(2,1:41),'bx'),
legend('East Position Measurement','North Position Measurement')
PlotBoi2('Time [sec]','Distance [m]',18)
%%
%%% ii. Run the KF
mu0a = [0; 85*cos(pi/4); 0; -85*sin(pi/4)]; %%initial mean of A at k=0
P0a =  900 * diag([10, 2, 10, 2],0); %%initial covar for A at k=0
muaHistfilt = zeros(4,length(tvec));
    muaHistfilt(:,1) = mu0a;
PaHistfilt = zeros(4,4,length(tvec));
    PaHistfilt(:,:,1) = P0a;
    
for kk=1:length(tvec)-1
    %%Prediction step
    mua_kkp1_minus = Fa * muaHistfilt(:,kk);
    Pa_kkp1_minus = Fa * PaHistfilt(:,:,kk) * Fa' + Qa;
    
    %%Kalman gain
    Kkkp1 = Pa_kkp1_minus * Ha' * inv(Ha*Pa_kkp1_minus*Ha' + Ra);
    
    %%Measurement update
    ykkp1 = yaHist(:,kk);
    mua_kkp1_plus = mua_kkp1_minus + Kkkp1*(ykkp1 - Ha*mua_kkp1_minus);
    Pa_kkp1_plus = (eye(4) - Kkkp1*Ha)*Pa_kkp1_minus;
    muaHistfilt(:,kk+1) = mua_kkp1_plus;
    PaHistfilt(:,:,kk+1) = Pa_kkp1_plus;
end

figure(),
%plot state errors
errora_Hist = muaHistfilt-xasingle_truth;
for ss=1:4
    eval(['subplot(41',num2str(ss),')'])
    plot(tvec,errora_Hist(ss,:),'b'),hold on
    plot(tvec,2*sqrt(squeeze(PaHistfilt(ss,ss,:)))','b--'),
    plot(tvec,-2*sqrt(squeeze(PaHistfilt(ss,ss,:)))','b--')
    legend('2-Sigma Error Bounds','Estimated State Error')
    xlabel('Time [sec]')
    ylabel('Error')    
end
fprintf('cii:\n')
fprintf('The position error is appears more constant over time\n')
fprintf('But the position and velocity states seem nearly as certain\n')

%% d. combined A and B tracking
%%% xadouble_truth; %---> get from data log
%%% xbdouble_truth; %---> get from data log
mu0a = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
P0a = 900 * diag([10, 2, 10, 2],0);
mu0b = [3200; 85*cos(pi/4); 3200; -85*sin(pi/4)];
P0b = 900 * diag([11, 4, 11, 4],0);
%%%i. simulate new ground-to-A and A-to-B transponder measurements 
Rd = [10 0.15; 0.15 10]; % d = difference in positions
ya2Hist = zeros(2,length(tvec)-1); % m^2
ydHist = zeros(2,length(tvec)-1);
for kk=2:length(tvec)
    vakk_samp = mvnrnd(zeros(2,1),Ra,1)'; %%simulate random ground station measurement noise
    ya2Hist(:,kk-1) = Ha*(xadouble_truth(:,kk-1)) + vakk_samp; %%simulate the noisy ground measurement
    vkkd_samp = mvnrnd(zeros(2,1),Rd,1)'; %%simulate random transponder noise
    ydHist(:,kk-1) = xadouble_truth([1 3],kk-1)-xbdouble_truth([1 3],kk-1) + vkkd_samp; %%simulate the noisy transponder measurement
end

Hdiff = [1 0 0 0 -1 0 0 0;0 0 1 0 0 0 -1 0];% 2x8 (performing a-b)

%%New augmented DT matrices
Faug = [Fa zeros(4,4); zeros(4,4) Fb];
Qaug = [Qa zeros(4,4); zeros(4,4) Qb];
Haug = [Ha zeros(2,4); Hdiff];
Raug = [Ra zeros(2,2); zeros(2,2) Rd];
mu0aug = [mu0a; mu0b]; %%initial augmented state mean at k=0
P0aug = [P0a zeros(4,4); zeros(4,4) P0b];  %%initial augmented state covar at k=0

%%Augmented measurements
yaugHist = [ya2Hist; ydHist];

%%Initialize and run the KF on augmented states
%%%% Filter with both sets of measurements
muaugHistfilt = zeros(8,length(tvec));
    muaugHistfilt(:,1) = mu0aug;
PaugHistfilt = zeros(8,8,length(tvec));
    PaugHistfilt(:,:,1) = P0aug;
    
%%%Filter with only yd set of measurements
muaugHistfilt_ydo = zeros(8,length(tvec));
    muaugHistfilt_ydo(:,1) = mu0aug;
PaugHistfilt_ydo = zeros(8,8,length(tvec));
    PaugHistfilt_ydo(:,:,1) = P0aug;

for kk=1:length(tvec)-1
    %%Prediction steps
    %%%%for filter with both sets of measurements
    muaug_kkp1_minus = Faug * muaugHistfilt(:,kk);
    Paug_kkp1_minus = Faug * PaugHistfilt(:,:,kk) * Faug' + Qaug;%Fa * PaHistfilt(:,:,kk) * Fa' + Qa;
    %%%%for filter with only transponder measurements
    muaug_ydo_kkp1_minus = Faug * muaugHistfilt_ydo(:,kk);
    Paug_ydo_kkp1_minus = Faug * PaugHistfilt_ydo(:,:,kk) * Faug' + Qaug;
    
    
    %%Kalman gain
    %%%%for filter with both sets of measurements
    Kkkp1_aug = Paug_kkp1_minus * Haug' * inv(Haug*Paug_kkp1_minus*Haug' + Raug);
    %%%%for filter with only transponder measurements
    Kkkp1_ydo_aug = Paug_ydo_kkp1_minus * Hdiff' * inv(Hdiff*Paug_ydo_kkp1_minus*Hdiff' + Rd);
    
    %%Measurement update
    %%%%for filter with both sets of measurements, y = 4x1, East & North
    %%%%for a & d
    ykkp1 = yaugHist(:,kk);
    muaug_kkp1_plus = muaug_kkp1_minus + Kkkp1_aug*(ykkp1 - Haug*muaug_kkp1_minus);
    Paug_kkp1_plus = (eye(8) - Kkkp1_aug*Haug)*Paug_kkp1_minus;
    muaugHistfilt(:,kk+1) = muaug_kkp1_plus;
    PaugHistfilt(:,:,kk+1) = Paug_kkp1_plus;
   
    %%%%for filter with only transponder measurements, y = 2x1, East &
    %%%%North for d
    ykkp1_ydo = ydHist(:,kk);
    muaug_ydo_kkp1_plus = muaug_ydo_kkp1_minus + Kkkp1_ydo_aug*(ykkp1_ydo - Hdiff*muaug_ydo_kkp1_minus);
    Paug_ydo_kkp1_plus = (eye(8) - Kkkp1_ydo_aug*Hdiff)*Paug_ydo_kkp1_minus;
    muaugHistfilt_ydo(:,kk+1) = muaug_ydo_kkp1_plus;
    PaugHistfilt_ydo(:,:,kk+1) = Paug_ydo_kkp1_plus;
end

%plot state errors for KF using both sets of measurements (di)
figure(),
erroraug_Hist = muaugHistfilt - [xadouble_truth; xbdouble_truth];
for ss=1:8
    eval(['subplot(42',num2str(ss),')'])
    plot(tvec,erroraug_Hist(ss,:),'b'),hold on
    plot(tvec,2*sqrt(squeeze(PaugHistfilt(ss,ss,:)))','b--'),
    plot(tvec,-2*sqrt(squeeze(PaugHistfilt(ss,ss,:)))','b--')
    ylabel(['ss=',num2str(ss)])
end

%plot state errors for KF using only transponder data (dii)
figure(),
erroraug_ydo_Hist = muaugHistfilt_ydo - [xadouble_truth; xbdouble_truth];
for ss=1:8
    eval(['subplot(42',num2str(ss),')'])
    plot(tvec,erroraug_ydo_Hist(ss,:),'r'),hold on
    plot(tvec,2*sqrt(squeeze(PaugHistfilt_ydo(ss,ss,:)))','r--'),
    plot(tvec,-2*sqrt(squeeze(PaugHistfilt_ydo(ss,ss,:)))','r--')
    ylabel(['ss=',num2str(ss)])
end
fprintf('di:\n')
fprintf('BE SURE TO EXPLAIN HOW YOU GOT THESE MATRICES\n')
fprintf('P0aug, Faug, Qaug are blocks b/c x is 8x1\n')
fprintf('Got Hdiff by relating transponder measurements to A and B states\n\n')
fprintf('dii:\n')
fprintf('With just transponder it''s significantly worse\n')
fprintf('b/c it''s only one measurement set instead of two...less information to update estimate\n')
fprintf('the Kalman gain has more data to manipulate in di\n')

%% d.iii. Compute probability of collision in case of augmented and unaugmented data

%%%Case 1: both ground station and transponder data available for KF
%%%estimates of aircraft A and B
ZR = 100; % m
ER = 100; % m
mvns_aug = zeros(size([20:80],1));
%%% Iterating
for kk=20:80   
    mu = muaugHistfilt(:,kk);
    mu = [(mu(1)-mu(5)) (mu(3)-mu(7))]';
    P = PaugHistfilt(:,:,kk);
    P12 = P(1:4,5:8);
    P21 = P(5:8,1:4);
    P11 = P(1:4,1:4);
    P22 = P(5:8,5:8);
    P = P11 + P22 - P12 - P21;
    P = 0.5 * (P + P');
    P = [P(1,1) P(1,3); P(3,1) P(3,3)];
    %%% Calculating probablity of collision
    mvns_aug(kk-19,1) = mvncdf(-[ZR;ER],[ZR;ER],mu,P);
end
figure
hold on
plot(tvec(20:80),mvns_aug,'-kx','markersize',14,'linewidth',2)
%%%Case 2: only transponder data available for KF
%%%estimates of aircraft A and B
mvns_d = zeros(size([20:80],1));
for kk=20:80   
    mu = muaugHistfilt_ydo(:,kk);
    mu = [(mu(1)-mu(5)) (mu(3)-mu(7))]';
    P = PaugHistfilt_ydo(:,:,kk);
    P12 = P(1:4,5:8);
    P21 = P(5:8,1:4);
    P11 = P(1:4,1:4);
    P22 = P(5:8,5:8);
    P = P11 + P22 - P12 - P21;
    P = 0.5 * (P + P');
    P = [P(1,1) P(1,3); P(3,1) P(3,3)];
    %%% Calculating probablity of collision
    mvns_d(kk-19,1) = mvncdf(-[ZR;ER],[ZR;ER],mu,P);
end

plot(tvec(20:80),mvns_d,'--ro','markersize',5,'linewidth',2)
PlotBoi2('Time, [sec]','Probability of Collision',18)




