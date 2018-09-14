clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
% ========================================================================
%%% Run switches
% ========================================================================
run_DQMethod = 1;
run_Quest    = 1;
run_OLAE     = 1;

% ========================================================================
%%% Establishing Vectors
% ========================================================================
%%% Body
v1_B = [0.8273;  0.5541; -0.0920];   v1_B = v1_B./norm(v1_B);
v2_B = [-0.8285; 0.5522; -0.0955];  v2_B = v2_B./norm(v2_B);
v3_B = [0.2155;  0.5522;  0.8022];    v3_B = v3_B./norm(v3_B);
v4_B = [0.5570; -0.7442; -0.2884];  v4_B = v4_B./norm(v4_B);

%%% Inertial
v1_N = [-0.1517; -0.9669;  0.2050];  v1_N = v1_N./norm(v1_N);
v2_N = [-0.8393;  0.4494; -0.3044];  v2_N = v2_N./norm(v2_N);
v3_N = [-0.0886; -0.5856; -0.8000]; v3_N = v3_N./norm(v3_N);
v4_N = [0.8814;  -0.0303;  0.5202];   v4_N = v4_N./norm(v4_N);

% ========================================================================
%%% TRIAD Method
% ========================================================================
% -------------------------------------
%%% Assigning T frame from sensors
% -------------------------------------
%%% Set from v2
t1_Bv2 = v1_B;
t2_Bv2 = cross(v1_B, v2_B); t2_Bv2 = t2_Bv2./norm(t2_Bv2);
t3_Bv2 = cross(t1_Bv2,t2_Bv2);

BTv2 = [t1_Bv2, t2_Bv2, t3_Bv2];

%%% Set from v3
t1_Bv3 = v1_B;
t2_Bv3 = cross(v1_B, v3_B); t2_Bv3 = t2_Bv3./norm(t2_Bv3);
t3_Bv3 = cross(t1_Bv3,t2_Bv3);

BTv3 = [t1_Bv3, t2_Bv3, t3_Bv3];

%%% Set from v4
t1_Bv4 = v1_B;
t2_Bv4 = cross(v1_B, v4_B); t2_Bv4 = t2_Bv4./norm(t2_Bv4);
t3_Bv4 = cross(t1_Bv4,t2_Bv4);

BTv4 = [t1_Bv4, t2_Bv4, t3_Bv4];


% -------------------------------------
%%% Assigning T frame from inertial knowledge
% -------------------------------------
%%% Set from v2
t1_Nv2 = v1_N;
t2_Nv2 = cross(v1_N, v2_N); t2_Nv2 = t2_Nv2./norm(t2_Nv2);
t3_Nv2 = cross(t1_Nv2, t2_Nv2);

NTv2 = [t1_Nv2, t2_Nv2, t3_Nv2];

%%% Set from v3
t1_Nv3 = v1_N;
t2_Nv3 = cross(v1_N, v3_N); t2_Nv3 = t2_Nv3./norm(t2_Nv3);
t3_Nv3 = cross(t1_Nv3, t2_Nv3);

NTv3 = [t1_Nv3, t2_Nv3, t3_Nv3];

%%% Set from v4
t1_Nv4 = v1_N;
t2_Nv4 = cross(v1_N, v4_N); t2_Nv4 = t2_Nv4./norm(t2_Nv4);
t3_Nv4 = cross(t1_Nv4, t2_Nv4);

NTv4 = [t1_Nv4, t2_Nv4, t3_Nv4];
% -------------------------------------
%%% Computing BN
% -------------------------------------
BN_TRIAD = BTv2 * NTv2';
BNv2 = BN_TRIAD;
BNv3 = BTv3 * NTv3';
BNv4 = BTv4 * NTv4';

BN_TRIAD-BNv3;
BN_TRIAD-BNv4;




% ========================================================================
%%% Davenport's q-method
% ========================================================================
if run_DQMethod
w1 = 2;
w2 = 1;
w3 = 1;

%%% Components of K matrix
B = w1*v1_B*v1_N' + w2*v2_B*v2_N' + w3*v3_B*v3_N';
sig = trace(B);
S = B + B';
Z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];

%%% Assign K matrix
K = [sig, Z';...
     Z, S-sig*eye(3)];
 
%%% Pulling max eigenvalue and vector from K
[v,e] = eig(K);
[max_eVal, index] = max(diag(e));
max_eVec = v(:,index);

%%% Mapping these Euler parameters (max eigenvector) to DCM
[BN_Dq] = quat2DCM(max_eVec);

B_diff = BN_Dq*BN_TRIAD';
diff_Dq = acos(.5*(trace(B_diff)-1))
diff_Dq_deg = diff_Dq*180/pi

% diff_Dq = acos(.5*(trace(BN_TRIAD)-1)) - acos(.5*(trace(BN_Dq)-1))

end % Davenport-q

% ========================================================================
%%% QUEST Method
% ========================================================================
if run_Quest
w1 = 2;
w2 = 1;
w3 = 1;

%%% Components of K matrix
B = w1*v1_B*v1_N' + w2*v2_B*v2_N' + w3*v3_B*v3_N';
sig = trace(B);
S = B + B';
Z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];

%%% Assign K matrix
K = [sig, Z';...
     Z, S-sig*eye(3)];
 
%%% Finding max eigenvalue
lam = w1+w2+w3;
err = 10;
while err > 1e-8
    h = .001;
    f1 = det(K - lam*eye(4));
    f2 = det(K - (lam+h)*eye(4));
    fp = (f2-f1)/h;
    lam_new = lam - f1/fp;
    lam = lam_new;
    err = abs(f1);
end

%%% Solving for CRP (q vec)
q = inv((lam+sig)*eye(3)-S)*Z;

%%% Solving for DCM
[BN_Quest] = crp2DCM(q);

B_diff = BN_Quest*BN_TRIAD';
diff_Quest = acos(.5*(trace(B_diff)-1))
diff_Quest_deg = diff_Quest*180/pi

% diff_Quest = acos(.5*(trace(BN_TRIAD)-1)) - acos(.5*(trace(BN_Quest)-1))

end % QUEST

% ========================================================================
%%% OLAE Method
% ========================================================================
if run_OLAE
w1 = 2;
w2 = 1;
w3 = 1;

%%% Creating summation and difference matrices
s1 = v1_B + v1_N; s2 = v2_B + v2_N; s3 = v3_B + v3_N;
d1 = v1_B - v1_N; d2 = v2_B - v2_N; d3 = v3_B - v3_N; 

st1 = [0 -s1(3) s1(2); s1(3) 0 -s1(1); -s1(2) s1(1) 0];
st2 = [0 -s2(3) s2(2); s2(3) 0 -s2(1); -s2(2) s2(1) 0];
st3 = [0 -s3(3) s3(2); s3(3) 0 -s3(1); -s3(2) s3(1) 0];

%%% Stacking summation and difference matrices
d = [d1;d2;d3];
S = [st1;st2;st3];

%%% Creating weighting matrix
W = blkdiag(w1*eye(3),w2*eye(3),w3*eye(3));

%%% Solving for CRP attitude description
q = inv(S'*W*S)*S'*W*d;

%%% Finding DCM
[BN_OLAE] = crp2DCM(q);

B_diff = BN_OLAE*BN_TRIAD';
diff_OLAE = acos(.5*(trace(B_diff)-1))
diff_OLAE_deg = diff_OLAE*180/pi

% diff_OLAE = acos(.5*(trace(BN_TRIAD)-1)) - acos(.5*(trace(BN_OLAE)-1))
end % OLAE

















