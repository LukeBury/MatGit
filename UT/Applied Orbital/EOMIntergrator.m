function [ dY ] = EOMIntergrator( t,Y )
dY = zeros(6,1);

% Unpack the state vector
y = Y(1:3);
dy = Y(4:6);
% dv=zeros(2,3);
% dv(1)=v;
% dv(2)=v_dot;


u=398600.4415;

% Dynamics of the system
% Here's where the fun happens!
ddy = (-u/(norm(y)^3))*y;

% Output the DERIVATIVE of the state
dY(1:3) = dy;
dY(4:6) = ddy;

end


