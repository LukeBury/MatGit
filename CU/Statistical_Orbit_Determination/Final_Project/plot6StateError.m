%%% Inputs
% 1) Data set 1 [nt x 6]
% 2) Data set 2 [nt x 6]
function plot6StateError(Times,data1, data2)
figure
subplot(3,2,1)
plot(Times,data1(:,1) - data2(:,1),'.')
subplot(3,2,3)
plot(Times,data1(:,2) - data2(:,2),'.')
subplot(3,2,5)
plot(Times,data1(:,3) - data2(:,3),'.')

subplot(3,2,2)
plot(Times,data1(:,4) - data2(:,4),'.')
subplot(3,2,4)
plot(Times,data1(:,5) - data2(:,5),'.')
subplot(3,2,6)
plot(Times,data1(:,6) - data2(:,6),'.')

end