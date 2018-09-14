clear
clc
close all
load hw6problem5data.mat 
% Includes: dt, q_1, q_2, Rlin, uTraj, xLtrue, xRtrue0, XRTrueHist
%           Yfulllin, YpartialLID, YPartiallin
q_1 = 1;
q_2 = q_1;
F = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
F = [F, zeros(4,24); zeros(24,4), eye(24)]; % Modified F, 28x28
G = [.5*(dt^2) 0; dt 0; 0 .5*(dt^2); 0 dt];
G = [G; zeros(24,2)]; % Modified G, 28x2
C = [(1/3)*(dt^3) .5*(dt^2); .5*(dt^2) dt];
H = [-1 0 0 0; 0 0 -1 0]; % H = 2x4
H = [H;H;H;H;H;H;H;H;H;H;H;H]; % H = 24x4
H = [H, eye(24)]; % H = 24x28
Q = blkdiag(q_1*C, q_2*C);
Q = blkdiag(Q,zeros(24,24)); % Resizing Q, 28x28
R = blkdiag(Rlin,Rlin,Rlin,Rlin,Rlin,Rlin,Rlin,Rlin,Rlin,Rlin,Rlin,Rlin);
% Yfulllin = 200 sets of stacked and ordered 2x1 measurements to each LM
% xLtrue = true ordered positions of the landmarks (stack of 1x2 vecs)
% uTraj = control input history (u) (acceleration)
% q_1,q_2 = constant process noise intensities
% Rlin = measurement noise covariance
% xRTrueHist = State true values
xLTrue = reshape(xLtrue,[],1);
tvec = 0:199;
mu0a = [xRtrue0;xLTrue];
P0a = eye(28);

muaHistfilt = zeros(28,length(tvec));
    muaHistfilt(:,1) = mu0a;
PaHistfilt = zeros(28,28,length(tvec));
    PaHistfilt(:,:,1) = P0a;
figure(1)
hold all

for kk = 1:length(tvec)-1
    %%% 5iv)
    Ymeas = Ypartiallin{1,kk}; % [E;N;E;N;E;N]
    Ynum  = YpartialLID{1,kk}; % [1;  3;  12]
    Hn = [-1 0 0 0; 0 0 -1 0];
    Hn = repmat(Hn,size(Ynum,1),1);
    Hr = zeros(size(Hn,1),24);
    for i = 1:size(Ynum,1) % Looping through sensor numbers
        n = Ynum(i);
%         nn = n+1;
        Hr(2*i-1:2*i,2*n-1:2*n) = eye(2);
    end
    Hn = [Hn,Hr];
    R = Rlin;
    for i = 1:size(Ynum,1)-1
        R = blkdiag(R,Rlin);
    end

    if size(Ymeas) ~= [0 0] % If there is a measurement
    %%% Prediction Step
    mua_kkp1_minus = F*muaHistfilt(:,kk) + G*uTraj(:,kk); % 28x1
    Pa_kkp1_minus = F*PaHistfilt(:,:,kk)*F' + Q; % 28x28
    
    %%% Kalman Gain
    Kkkp1 = Pa_kkp1_minus * Hn' / (Hn*Pa_kkp1_minus*Hn' + R);

    %%% Measurement Update
%     ykkp1 = Yfulllin(:,kk); % 24x1
    ykkp1 = Ymeas;
    mua_kkp1_plus = mua_kkp1_minus + Kkkp1*(ykkp1 - Hn*mua_kkp1_minus);
    Pa_kkp1_plus = (eye(28) - Kkkp1*Hn)*Pa_kkp1_minus;
    muaHistfilt(:,kk+1) = mua_kkp1_plus;
    PaHistfilt(:,:,kk+1) = Pa_kkp1_plus;
    end
   
    if size(Ymeas) == [0 0] % If there is not a measurement
    %%% Prediction Step
    mua_kkp1_minus = F*muaHistfilt(:,kk) + G*uTraj(:,kk); % 28x1
    Pa_kkp1_minus = F*PaHistfilt(:,:,kk)*F' + Q; % 28x28
    
    muaHistfilt(:,kk+1) = mua_kkp1_minus;
    PaHistfilt(:,:,kk+1) = Pa_kkp1_minus;
    end
    
    
%     %% 5ii
%     %%% Prediction Step
%     mua_kkp1_minus = F*muaHistfilt(:,kk) + G*uTraj(:,kk); % 28x1
%     Pa_kkp1_minus = F*PaHistfilt(:,:,kk)*F' + Q; % 28x28
%     
%     %%% Kalman Gain
%     Kkkp1 = Pa_kkp1_minus * H' / (H*Pa_kkp1_minus*H' + R); % 28x24
%     
%     %%% Measurement Update
%     ykkp1 = Yfulllin(:,kk); % 24x1
%     mua_kkp1_plus = mua_kkp1_minus + Kkkp1*(ykkp1 - H*mua_kkp1_minus);
%     Pa_kkp1_plus = (eye(28) - Kkkp1*H)*Pa_kkp1_minus;
%     muaHistfilt(:,kk+1) = mua_kkp1_plus;
%     PaHistfilt(:,:,kk+1) = Pa_kkp1_plus;
%     
%     
    if mod(kk,10)==0 || kk == 1
        [Xell,Yell] = calc_gsigma_ellipse_plotpoints...
        ([muaHistfilt(1,kk),muaHistfilt(3,kk)]',...
         [PaHistfilt(1,1,kk),PaHistfilt(1,3,kk);PaHistfilt(3,1,kk),PaHistfilt(3,3,kk)],...
          2,100);
        plot(Xell,Yell,'b.-') %%plot new ellipse
        
    end
end

%Plotting Measurements
for i = 1:12
    p3 = plot(Yfulllin(2*i-1,:)+muaHistfilt(1,:), Yfulllin(2*i,:)+muaHistfilt(3,:),'c.');
end
p1 = plot(muaHistfilt(1,:),muaHistfilt(3,:),'linewidth',2);
p2 = plot(xRTrueHist(1,:),xRTrueHist(3,:),'linewidth',2);
p4 = plot(xLtrue(1,:),xLtrue(2,:),'kx','linewidth',2,'markersize',8);
legend([p1 p2 p4 p3],'KF','True','Landmarks','Measurements')
PlotBoi2('East Position, [m]','North Position, [m]',18)

%plot state errors
error_Hist = muaHistfilt(1:4,2:200)-xRTrueHist(:,1:199);
figure(2)

for ss=1:4
    
    eval(['subplot(41',num2str(ss),')'])
    plot(error_Hist(ss,:),'b'),hold on
    plot(2*sqrt(squeeze(PaHistfilt(ss,ss,:)))','b--'),
    plot(-2*sqrt(squeeze(PaHistfilt(ss,ss,:)))','b--')
    legend('Estimated State Error','2-Sigma Error Bounds')
    PlotBoi2('Time [sec]','Error',12)
%     if ss == 1
%         title('q = 1')
%     end
end



