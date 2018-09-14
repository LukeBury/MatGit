% Made for out-of-loop
%%% Inputs
% 1) Length(Times)
% 2) Pre-Fit Residuals [2xn]
% 3) Htildes [2x6xn]
% 4) xHats [6xn]
%%% Outputs
% 1) Post-Fit Residuals [2xn]
function [ ei ] = calcPostFitResiduals(n, yi, H, X)
if size(H,2) ~= size(X,1)
    warning('States must be in column-vector form')
    return
end
ei = zeros(size(yi,1),n);
if size(X,2) == 1 % Batch
    for k = 1:n
        ei(:,k) = yi(:,k) - H(:,:,k)*X;
    end
else              % Non-Batch
    for k = 1:n 
        ei(:,k) = yi(:,k) - H(:,:,k)*X(:,k);
    end
end

end

