clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cal0 = mean(raw_compressed);
cal1 = mean(raw_extended);

cal=[cal0,cal1];
delta = [0,.03026];

p=polyfit(cal,displacement,1);

fit = @(x) (p(1)*x + p(2)); %%% enter voltage, get displacement



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Building
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = [0:1099];

stroke1 = fit(raw_h1_pot);
stroke2 = fit(raw_h2_pot);
stroke3 = fit(raw_h3_pot);
stroke4 = fit(raw_h4_pot);
stroke5 = fit(raw_h5_pot);

acc1 = raw_h1_acc;
acc2 = raw_h2_acc;
acc3 = raw_h3_acc;
acc4 = raw_h4_acc;
acc5 = raw_h5_acc;

figure(1)
hold all
plot(t,stroke1)
plot(t,stroke2)
plot(t,stroke3)
plot(t,stroke4)
plot(t,stroke5)
legend('1','2','3','4','5')


figure(2)
plot(t,acc1)
figure(3)
plot(t,acc2)
figure(4)
plot(t,acc3) %%% the faster it's going, the  less powerful the force of friction - which is why each higher drop
                % shows a higher magnitude of acceleration
figure(5)
plot(t,acc4)
figure(6)
plot(t,acc5)
%legend('1','2','3','4','5')
















