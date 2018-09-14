function [ decimal ] = octal_repr( x )
%OCTAL_REPR Coverts binary vector to octal string representation
% LSB is first bit

N = length(x);
N = (floor(N / 3) + 1) * 3;
x(N) = 0;  % pads array with 0s at end
x = reshape(x, 3, N / 3)';
M = length(x);
digits = zeros(M, 1);
for i = 1:M
%     digits(i) = binvec2dec(x(i, :));
    str_x = num2str(x(i, :));
    digits(i) = bin2dec(str_x);

end
decimal = 0;
for i = 1:M
    decimal = decimal + digits(i) * 10^(i-1);
end

end

