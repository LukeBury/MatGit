% m = total number of measurements
function [ RMS ] = calculateRMS(n, eis, R, m)
temp = 0;
for k = 1:n
   temp = temp + eis(:,k)'*R*eis(:,k);
end
RMS = (temp/m)^(1/2);

end

