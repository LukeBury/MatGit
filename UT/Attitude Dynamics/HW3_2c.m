function HW3_2
clear
clc

alpha = 11; % November - degrees
beta  = 18; % degrees
gamma =  8; % '92

vA =[1;1;1]; % vector to be rotated

% Creating each of the rotation matrices
C1 = [1 0 0;0 cosd(beta) sind(beta);0 -sind(beta) cosd(beta)];
C2 = [cosd(alpha) 0 -sind(alpha);0 1 0;sind(alpha) 0 cosd(alpha)];
C3 = [cosd(gamma) sind(gamma) 0;-sind(gamma) cosd(gamma) 0;0 0 1];

% Combining rotation matrices 
C = C3 * C1 * C2;

phi = acosd((trace(C)-1)/2); %degrees

avec = [C(2,3)-C(3,2);C(3,1)-C(1,3);C(1,2)-C(2,1)];
mag_avec = norm(avec);

a = avec./mag_avec;

vB = qn_LGB(a,phi,vA);

%confirming that magnitude remains the same
if .99999999 < norm(vA) == norm(vB) < 1.00000001
    fprintf('The magnitude has remained constant\n');
else
    fprintf('The rotation has failed\n')
end

fprintf('\nAlso, Luke Bury rulez\n')
end

function [vB] = qn_LGB(a,phi,vA)

mag_a = norm(a); %acquiring magnitude of vector a

eig_a = a./mag_a; %turning a into unit eigenaxis

qBA = zeros(4,1); %building quaternion qBA

%%% establishing each piece of the quaternion
qBA(1) = eig_a(1)*sind(phi/2);
qBA(2) = eig_a(2)*sind(phi/2);
qBA(3) = eig_a(3)*sind(phi/2);
qBA(4) = cosd(phi/2);


%%% vB = qBA x Va x qBAs
%%% Preparing components for qBA x Va
Va = [vA(1);vA(2);vA(3);0];


eq = [qBA(1);qBA(2);qBA(3)];
eVa = [Va(1);Va(2);Va(3)];

nq = qBA(4);
nVa = Va(4);


%%% qBA x Va
qBAxVa = [eq.*nVa + eVa.*nq - cross(eq,eVa); nq*nVa - dot(eq,eVa)];

%%% Preparing components for (qBA x Va) x qBAs - where (qBA x Va) = [e1;n1]
qBAs= [-qBA(1);-qBA(2);-qBA(3);qBA(4)];

e1 = [qBAxVa(1);qBAxVa(2);qBAxVa(3)];
eqs = [qBAs(1);qBAs(2);qBAs(3)];

n1 = qBAxVa(4);
nqs = qBAs(4);

%%% (qBA x Va) x qBAs - where answer = vB

vB = [e1.*nqs + eqs.*n1 - cross(e1,eqs); n1*nqs - dot(e1,eqs)];

vB
end




