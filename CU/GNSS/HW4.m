clear
clc
close all

%% =======================================================================
%%% Problem 2-1
%%%-----------------------------------------------------------------------

% Assigning Initial Shift Registers
G1 = ones(1,10);
G2 = ones(1,10);

% Assigning phase selector
phaseSelector = [2,6];

% Initializing CA Code
code1 = [];

% Define chip length
chipLength = 1023;

for cycle = 1:16
    if cycle ~= 1
        % Update shift registers
        G1 = G1gen(G1);
        G2 = G2gen(G2);
    end
    
    %%% G1
    G1out = G1(end)

    %%% G2
    G2out = mod(G2(phaseSelector(1)) + G2(phaseSelector(2)),2);
    
    %%% Updating code
    newbit = mod(G1out + G2out,2);
    code1 = [code1, newbit];
end

code1(1:16)'
% code1(end-16:end)
return

code2 = [];
for cycle = chipLength+1:2*chipLength
    if cycle ~= 1
        % Update shift registers
        G1 = G1gen(G1);
        G2 = G2gen(G2);
    end
    
    %%% G1
    G1out = G1(end);

    %%% G2
    G2out = mod(G2(phaseSelector(1)) + G2(phaseSelector(2)),2);
    
    %%% Updating code
    newbit = mod(G1out + G2out,2);
    code2 = [code2, newbit];
end

% code2(1:16)
% code2(end-16:end)


%% =======================================================================
%%% Problem 2-2
%%%-----------------------------------------------------------------------
% % Converting 0 -> 1, and 1 -> -1
% code1(code1==1)= -1;
% code1(code1==0)= 1;
% 
% R19 = 0;
% for cycle = 1:length(code1)
%     i = cycle-1;
%     R19 = R19 + 
% R19 = (1/1023)



%% =======================================================================
%%% Functions
%%%-----------------------------------------------------------------------
function [newregister] = G1gen(G1)
newbit = mod(G1(3) + G1(10),2);
% if G1(3) + G1(10) == 2 || G1(3) + G1(10) == 0
%     newbit = 0;
% else
%     newbit = 1;
% end
newregister = [newbit, G1(1:9)];
end

function [newregister] = G2gen(G2)
newbit = mod(G2(2) + G2(3) + G2(6) + G2(8) + G2(9) + G2(10),2);
newregister = [newbit, G2(1:9)];
end














