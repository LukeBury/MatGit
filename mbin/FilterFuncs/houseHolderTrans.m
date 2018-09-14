%%% Inputs:
% 1) Matrix to be transformed
%%% Outputs:
% 2) Transformed matrix
function [ A ] = houseHolderTrans( A )
n = size(A,2)-1;
m = size(A,1)-n;
for k = 1:n
    temp = 0;
    u = zeros(n+m,1);
    for i = k:m+n
        temp = temp + A(i,k)^2;
    end
    sig = sign(A(k,k)) * sqrt(temp);
    u(k) = A(k,k) + sig;
    A(k,k) = -sig;
    for i = k+1:m+n
        u(i,1) = A(i,k);
    end
    B = 1/(sig*u(k,1));
    for j = k+1:n+1
        temp = 0;
        for i = k:m+n
            temp = temp + u(i,1)*A(i,j);
        end
        gamma = B * temp;
        for i = k:m+n
            A(i,j) = A(i,j) - gamma*u(i,1);
        end
        for i = k+1:m+n
            A(i,k) = 0;
        end
    end

end
end

