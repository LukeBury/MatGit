clear
clc

% 
% m = [3 5 7 4 3 8];
% 
% sum1=0;
% sum2=0;
% 
% for k=1:length(m)
%     sum1 = sum1+m(k);
% end
% 
% for k=1:length(m)
%     sum2 = sum2+m(k);
% end
% 
% prod1 = sum1*sum2
% 
% sum3=0;
% 
% for i=1:length(m)
%     for k=1:length(m)
%         sum3 = sum3 + m(i)*m(k);
%     end
% end
% 
% sum3
% 
% sum4=0;
% 
% for i=1:length(m)
%     for k=1:length(m)
%         sum4 = sum4 + sqrt(m(i))*sqrt(m(k));
%     end
% end
% 
% sum4

% 
% 
% r = [2 7 5 3 8 8];
% summ1 = 0;
% summ2 = 0;
% for i = 1:length(m)
%     summ1 = summ1 + m(i)*r(i)*r(i);
%     summ2 = summ2 + r(i)*r(i);
% end



a = [9 3 2 5 6 7 8 3];
n = length(a);
sum1 = 0;
sum3 = 0;

for i = 1:n
    for j = 1:n
        if j == i
            continue
        end
        sum1 = sum1 + a(i)*a(j);
    end
end


for i = 1:(n-1)
    for j = (i+1):n
        sum3 = sum3 + a(i)*a(j);
    end
end


sum1
sum3



