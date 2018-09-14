clear
clc
% 
% b = 0.5;
% 
% a = 0.5;
% 
% for kk = 1:1000
%     fprintf('-------%2.0f\n',kk)
%     M = 1/(1-a)
%     a = b/M
%     
% end

% ------------------------------------------------------------------------
%%% 4.6
% ------------------------------------------------------------------------

matA = [-2 -1; 3 2];
[vA, eA] = eig(matA);
dA = det(matA);

matB = [2 0; 0 2];
[vB, eB] = eig(matB);
dB = det(matB);

matC = [-5 -2; 5 1];
[vC, eC] = eig(matC);
dC = det(matC);

matD = [2 1; 1 2];
[vD, eD] = eig(matD);
dD = det(matD);

matE = [7 -10; 5 -8];
[vE, eE] = eig(matE);
dE = det(matE);

matF = [3 1; -1 1];
[vF, eF] = eig(matF);
dF = det(matF);

matG = [-5 1; -6 0];
[vG, eG] = eig(matG);
dG = det(matG);

matH = [1 0; -2 -1];
[vH, eH] = eig(matH);
dH = det(matH);

% No determinants are = 0


