%Rodrigues rotation about given axis and angle

function [vout,rotu]=rotuv(v,u,theta)
%INPUTS
%v(3x1): vector to be rotated
%u(3x1): vector axis of rotation
%theta(1): (rad)angle of rotation about u axis

%OUTPUTS
%vout(3x1): rotated vector
%rotu(3x3): rotation matrix
w=[ 0 -u(3) u(2);
   u(3) 0  -u(1);
  -u(2) u(1) 0 ];
rotu=eye(3)+w*sin(theta)+w*w*(1-cos(theta));
vout=rotu*v;
end