function [z] = circular_correlation(x, y)
    N = length(x);
    z = zeros(N, 1);
    for n = 1:N
        z(n) = sum(x .* circshift(y, -n));
    end
end
