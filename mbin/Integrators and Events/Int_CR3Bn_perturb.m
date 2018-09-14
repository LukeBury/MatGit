function [dX] = Int_CR3Bn_perturb(t,X,u,R1,R2,J21,J22,J31,J32,J41,J42)
%%% For numerical integration in the normalized CR3BP with J2 of each body
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          u - mass ratio of CR3BP system
%          R1 - radius of primary body
%          R2 - radius of secondary body
%          J21 - J2 of primary body
%          J22 - J2 of secondary body
%          J31 - J3 of primary body
%          J32 - J3 of secondary body
%          J41 - J4 of primary body
%          J42 - J4 of secondary body
% =======================================================================
%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Define position of bodies
x1 = -u;
x2 = 1-u;

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x+u)^2 + y^2 + z^2);
r2 = sqrt((x+u-1)^2 + y^2 + z^2);

%%% Determine perturbations to be included and calculate effect
% gamma21 = 0; gamma22 = 0; gamma21_z = 0; gamma22_z = 0;
% gamma31 = 0; gamma32 = 0; gamma31_z = 0; gamma32_z = 0;
% gamma41 = 0; gamma42 = 0; gamma41_z = 0; gamma42_z = 0;
% for kk = 1:length(perturbations)
%     if isequal(perturbations{kk},'J21') == 1
%         gamma21 = 3*(1-u)*R1^2*J21*(5*z^2-r1^2)/(2*r1^7);
%         gamma21_z = 3*(1-u)*R1^2*J21*(5*z^2-3*r1^2)/(2*r1^7);
%         
%     elseif isequal(perturbations{kk},'J22') == 1
%         gamma22 = 3*u*R2^2*J22*(5*z^2-r2^2)/(2*r2^7);
%         gamma22_z = 3*u*R2^2*J22*(5*z^2-3*r2^2)/(2*r2^7);
%     
%     elseif isequal(perturbations{kk},'J31') == 1
%         gamma31 = R1^3*J31*(1-u)*5*z*(7*z^2-3*r1^2)/(2*r1^9);
%         gamma31_z = R1^3*J31*(1-u)*(3*r1^4-30*r1^2*z^2+35*z^4)/(2*r1^9);
%     
%     elseif isequal(perturbations{kk},'J32') == 1
%         gamma32 = R2^3*J32*u*5*z*(7*z^2-3*r2^2)/(2*r2^9);
%         gamma32_z = R2^3*J32*u*(3*r2^4-30*r2^2*z^2+35*z^4)/(2*r2^9);
%     
%     elseif isequal(perturbations{kk},'J41') == 1
%         gamma41 = R1^4*J41*(1-u)*15*(r1^4-14*r1^2*z^2+21*z^4)/(8*r1^11);
%         gamma41_z = R1^4*J41*(1-u)*5*(15*r1^4-70*r1^2*z^2+63*z^4)/(8*r1^11);
%     
%     elseif isequal(perturbations{kk},'J42') == 1
%         gamma42 = R2^4*J42*u*15*(r2^4-14*r2^2*z^2+21*z^4)/(8*r2^11);
%         gamma42_z = R2^4*J42*u*5*(15*r2^4-70*r2^2*z^2+63*z^4)/(8*r2^11);
%     end
% end

%%% Calculating effects of J2, J3, J4
gamma21 = 3*(1-u)*R1^2*J21*(5*z^2-r1^2)/(2*r1^7);
gamma21_z = 3*(1-u)*R1^2*J21*(5*z^2-3*r1^2)/(2*r1^7);

gamma22 = 3*u*R2^2*J22*(5*z^2-r2^2)/(2*r2^7);
gamma22_z = 3*u*R2^2*J22*(5*z^2-3*r2^2)/(2*r2^7);

gamma31 = R1^3*J31*(1-u)*5*z*(7*z^2-3*r1^2)/(2*r1^9);
gamma31_z = R1^3*J31*(1-u)*(3*r1^4-30*r1^2*z^2+35*z^4)/(2*r1^9);

gamma32 = R2^3*J32*u*5*z*(7*z^2-3*r2^2)/(2*r2^9);
gamma32_z = R2^3*J32*u*(3*r2^4-30*r2^2*z^2+35*z^4)/(2*r2^9);

gamma41 = R1^4*J41*(1-u)*15*(r1^4-14*r1^2*z^2+21*z^4)/(8*r1^11);
gamma41_z = R1^4*J41*(1-u)*5*(15*r1^4-70*r1^2*z^2+63*z^4)/(8*r1^11);

gamma42 = R2^4*J42*u*15*(r2^4-14*r2^2*z^2+21*z^4)/(8*r2^11);
gamma42_z = R2^4*J42*u*5*(15*r2^4-70*r2^2*z^2+63*z^4)/(8*r2^11);

%%% Equations of Motion
ddx = 2*dy + x + ((1-u)/r1^3-gamma21 - gamma31 - gamma41)*(x1-x) + (u/r2^3 - gamma22 - gamma32 - gamma42)*(x2-x);
ddy = -2*dx + y - y*((1-u)/r1^3 + u/r2^3 - gamma21 - gamma22 - gamma31 - gamma32 - gamma41 - gamma42);
ddz = z*(-(1-u)/r1^3 - u/r2^3 + gamma21_z + gamma22_z + gamma41_z + gamma42_z) + (gamma31_z+gamma32_z);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end
