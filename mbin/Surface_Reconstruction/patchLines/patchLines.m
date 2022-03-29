function [Xp,Yp,Zp] = patchLines(x1,x2)
%PATCHLINES creates patch data to create 3D surface between two 3D lines
%   [Xp,Yp,Zp] = PATCHLINES(x1,x2) creates n-1 patches where 
%   n = size(x1,2). Data can then be plot using patch(Xp,Yp,Zp,'b')
%   Ian Elliott
if (size(x1, 1) ~= size(x2, 1))
    error('Lines must be equal length')
end

% Get length of lines
n = size(x1, 1);

% Create patch points
Xp = nan(4, n - 1);
Yp = nan(4, n - 1);
Zp = nan(4, n - 1);
for ii = 1:(n - 1)  
    r1 = x1(ii, 1:3);
    r2 = x1(ii + 1, 1:3);
    r3 = x2(ii, 1:3);
    r4 = x2(ii + 1, 1:3);
    
    Xp(:, ii) = [r1(1) r2(1) r4(1) r3(1)]';
    Yp(:, ii) = [r1(2) r2(2) r4(2) r3(2)]';
    Zp(:, ii) = [r1(3) r2(3) r4(3) r3(3)]';
end

end
