function L = collinearEquilibriumPoints_ZH_C22s(prms)
%EquilibriumPoints calculates the equilibrium points of a CR3BP system.
%   EquilibriumPoints inputs a mass ration for a CR3BP system and a scalar
%   or vector of integer identifies for each of the L1-5 equilibrium points
%   to return. Collinear equilibrium points are solved for using a
%   Newton-Raphson method of root finding.
%
%   L = EquilibriumPoints(MU,POINTS) inputs the mass ratio, MU, and POINTS,
%   a single integer identifying which equilibrium point to calculate. The
%   function returns L, a 1x3 vector corresponding to the x, y and z 
%   location of the equilibrium point in the non-dimensional, rotating  
%   frame. If POINTS is input as a vector of integers of size n identifying 
%   which equilibrium points to calculate. The function returns L, an nx3  
%   matrix corresponding to the z and y locations of each identified 
%   equilibrium point in the  non-dimensional, rotating frame.
%
%   L = EquilibriumPoints(MU) returns the coordinates of all five
%   equilibrium points.
%   
%   Ian Elliott - 09/12/2017
%   Luke Bury   - 03/28/2018
%

% Initialize column vectors
x_colVec = zeros(3,1);
y_colVec = zeros(3,1);
z_colVec = zeros(3,1);

%%% Unpacking some parameters (for brevity/clarity)
mu  = prms.u;

if isfield(prms, 'R1')
    R1 = prms.R1;
end
if isfield(prms, 'R2')
    R2 = prms.R2;
end

if isfield(prms, 'J2p')
    J2p = prms.J2p;
    dUdx_J2p     = @(x) -3*(1-mu)*(R1^2)*J2p*(x+mu) / (2*abs((x+mu)^5));
    d2Udx2_J2p   = @(x) 6*(1-mu)*(R1^2)*J2p / abs((x+mu)^5);
else
    dUdx_J2p   = @(x) 0;
    d2Udx2_J2p = @(x) 0;
end

if isfield(prms, 'J2s')
    J2s = prms.J2s;
    dUdx_J2s     = @(x) -3*mu*(R2^2)*J2s*(x-1+mu) / (2*abs((x-1+mu)^5));
    d2Udx2_J2s   = @(x) 6*mu*(R2^2)*J2s / abs((x-1+mu)^5);
else
    dUdx_J2s   = @(x) 0;
    d2Udx2_J2s = @(x) 0;
end

%%% J3 has no effect
% % % if isfield(prms, 'J3p')
% % %     dUdx_J3p     = @(mu, R1, J3p, x) ;
% % %     d2Udx2_J3p   = @(mu, R1, J3p, x) ;
% % % else
% % %     dUdx_J3p   = @(mu, R1, J2p, x) 0;
% % %     d2Udx2_J3p = @(mu, R1, J2p, x) 0;
% % % end
% % % 
% % % if isfield(prms, 'J3s')
% % %     gamma_J3s = 5*mu*(prms.R2^3)*prms.J3s*X(3)*(7*X(3)^2 - 3*r2^2) / (2*r2^9);
% % % else
% % %     gamma_J3s = 0;
% % % end

if isfield(prms, 'J4p')
    J4p = prms.J4p;
    dUdx_J4p     = @(x) 15*(1-mu)*(R1^4)*J4p*(x+mu) / (8*abs((x+mu)^7));
    d2Udx2_J4p   = @(x) -45*(1-mu)*(R1^4)*J4p / (4*abs((x+mu)^7));
else
    dUdx_J4p   = @(x) 0;
    d2Udx2_J4p = @(x) 0;
end

if isfield(prms, 'J4s')
    J4s = prms.J4s;
    dUdx_J4s     = @(x) 15*mu*(R2^4)*J4s*(x-1+mu) / (8*abs((x-1+mu)^7));
    d2Udx2_J4s   = @(x) -45*mu*(R2^4)*J4s / (4*abs((x-1+mu)^7));
else
    dUdx_J4s   = @(x) 0;
    d2Udx2_J4s = @(x) 0;
end

%%% J5 has no effect
% % % if isfield(prms, 'J5p')
% % %     gamma_J5p = 7*(1-u)*(prms.R1^5)*J5p*X(3)*(15*r1^4 - 90*X(3)^2*r1^2 + 99*X(3)^4) * (x+u) / (8*r1^13);
% % % else
% % %     gamma_J5p = 0;
% % % end
% % % 
% % % if isfield(prms, 'J5s')
% % %     gamma_J5s = 7*mu*(prms.R2^5)*prms.J5s*X(3)*(15*r2^4 - 90*X(3)^2*r2^2 + 99*X(3)^4) / (8*r2^13);
% % % else
% % %     gamma_J5s = 0;
% % % end

if isfield(prms, 'J6p')
    J6p = prms.J6p;
    dUdx_J6p     = @(x) -35*(1-mu)*(R1^6)*J6p*(x+mu) / (16*abs((x+mu)^9));
    d2Udx2_J6p   = @(x) 35*(1-mu)*(R1^6)*J6p / (2*abs((x+mu)^9));
else
    dUdx_J6p   = @(x) 0;
    d2Udx2_J6p = @(x) 0;
end

if isfield(prms, 'J6s')
    J6s = prms.J6s;
    dUdx_J6s     = @(x) -35*mu*(R2^6)*J6s*(x-1+mu) / (16*abs((x-1+mu)^9));
    d2Udx2_J6s   = @(x) 35*mu*(R2^6)*J6s / (2*abs((x-1+mu)^9));
else
    dUdx_J6s   = @(x) 0;
    d2Udx2_J6s = @(x) 0;
end

if isfield(prms, 'C22s')
    dUdx_C22s     = @(x) -9*mu*R2*R2*prms.C22s*(x-1+mu) / abs((x-1+mu)^5);
    d2Udx2_C22s   = @(x) 36*mu*R2*R2*prms.C22s / abs((x-1+mu)^5);
else
    dUdx_C22s   = @(x) 0;
    d2Udx2_C22s = @(x) 0;
end

dUdx   = @(x) prms.n*prms.n*x - (1-mu)*(x+mu)/abs(x+mu)^3 - mu*(x-1+mu)/abs(x-1+mu)^3 + dUdx_J2p(x) +   dUdx_J2s(x) +   dUdx_J4p(x) +   dUdx_J4s(x) +   dUdx_J6p(x) +   dUdx_J6s(x)   + dUdx_C22s(x);
d2Udx2 = @(x) prms.n*prms.n + 2*(1-mu)/abs(x+mu)^3 + 2*mu/abs(x-1+mu)^3 +               d2Udx2_J2p(x) + d2Udx2_J2s(x) + d2Udx2_J4p(x) + d2Udx2_J4s(x) + d2Udx2_J6p(x) + d2Udx2_J6s(x) + d2Udx2_C22s(x);

points = 1:3;
for ii = 1:3
    switch points(ii)
        % L1, between primaries
        case 1      
            x_colVec(ii) = NewtonMethod(dUdx, d2Udx2, 1/2-mu, 1e-15, 500);
            y_colVec(ii) = 0;
            
        % L2, right of p2   
        case 2
            x_colVec(ii) = NewtonMethod(dUdx, d2Udx2, (1-mu)+1e-5, 1e-15, 500); 
            y_colVec(ii) = 0;
            
        % L3, left of p1
        case 3
            x_colVec(ii) = NewtonMethod(dUdx, d2Udx2, -(1-mu), 1e-15, 500); 

            y_colVec(ii) = 0;
        
    end
end
    
% Return cooridinates of equilibrium points
L = [x_colVec,y_colVec,z_colVec];

if length(unique(L(:,1)))~=3
    warning('Don''t have 3 unique solutions')
end

if sign(d2Udx2(x_colVec(1))) == -1
    warning('L1 value is local maxima not minima')
end
if sign(d2Udx2(x_colVec(2))) == -1
    warning('L2 value is local maxima not minima')
end
if sign(d2Udx2(x_colVec(3))) == -1
    warning('L3 value is local maxima not minima')
end

end








