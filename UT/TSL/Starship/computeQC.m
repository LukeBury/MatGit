function [ qC ] = computeQC( q0, qDes )

% calculate intermediate quaternion after rotation about axis 1
qE = quaternionError(qDes, q0);
anglesE = q_to_euler(qE);
e1 = [1;0;0];
qE1 = [e1*sin(anglesE(1)/2); cos(anglesE(1)/2)];
qC1 = R_to_q(q_to_R(qE1)'*q_to_R(q0));

% calculate intermediate quaternion after rotation about axis 2
qE = quaternionError(qDes, qC1);
anglesE = q_to_euler(qE);
e2 = [0;1;0];
qE2 = [e2*sin(anglesE(2)/2); cos(anglesE(2)/2)];
qC2 = R_to_q(q_to_R(qE2)'*q_to_R(qC1));

% calculate intermediate quaternion after rotation about axis 3
qE = quaternionError(qDes, qC2);
anglesE = q_to_euler(qE);
e3 = [0;0;1];
qE3 = [e3*sin(anglesE(3)/2); cos(anglesE(3)/2)];
qC3 = R_to_q(q_to_R(qE3)'*q_to_R(qC2));

qC = [qC1 qC2 qC3]; % vector of intermediate commanded quaternions

end

