function [ D ] = bin2oct( B )
% Binary to Octal Conversion
% Input   : B     B is a row vector representing a binary sequence
% Output  : D     Binary input in octal format
[t,S]=size(B);
B1=fliplr(B);
K=1;
for i=1:3:S
    if((i+2)<=S)
        b=B1(1,i:i+2);
        Dec=b(1,1)*2^0+b(1,2)*2^1+b(1,3)*2^2;
        Conv(1,K)=Dec;
        K=K+1;
    
    
    elseif((i+2)<=S-1)
        b=B1(1,i:i+1);
        Dec=b(1,1)*2^0+b(1,2)*2^1;
        Conv(1,K)=Dec;
        K=K+1;
    else
        b=B1(1,i:i);
        Dec=b(1,1)*2^0;
        Conv(1,K)=Dec;
        K=K+1;
    end
    
end
D=fliplr(Conv);

% Putting into octal format
D = num2str(D);
D(isspace(D)) = '';
D = str2num(D);
end


