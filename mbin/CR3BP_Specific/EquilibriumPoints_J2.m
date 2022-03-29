function L = EquilibriumPoints_J2(u,n,J21,J22,R1_n,R2_n,points)
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
if nargin < 6, points = 1:3; end
numPoints = length(points);

% Initialize column vectors
x = zeros(numPoints,1);
y = zeros(numPoints,1);
z = zeros(numPoints,1);

% Collinear equilibrium points 
dUdx = @(x) n*n*x - (1-u)*(x+u)/abs(x+u)^3 - u*(x-1+u)/abs(x-1+u)^3  - 3*J21*R1_n^2*(1-u)*(x+u)/(2*abs(x+u)^5)- 3*J22*R2_n^2*u*(x-1+u)/(2*abs(x-1+u)^5);
d2Udx2 = @(x) n*n + 2*(1-u)/abs(x+u)^3 + 2*u/abs(x-1+u)^3 + 6*J21*R1_n^2*(1-u)/abs(x+u)^5 + 6*J22*R2_n^2*u/abs(x-1+u)^5;

for ii = 1:numPoints
    switch points(ii)
        % L1, between primaries
        case 1      
            x(ii) = NewtonMethod(dUdx, d2Udx2, 1/2-u);
            y(ii) = 0;
            
        % L2, right of p2   
        case 2
            x(ii) = NewtonMethod(dUdx, d2Udx2, (1-u)+1e-5); 
            y(ii) = 0;
            
        % L3, left of p1
        case 3
            x(ii) = NewtonMethod(dUdx, d2Udx2, 1.25*(-u)); 
            y(ii) = 0;
        
% %         % L4, upper triangular point
% %         case 4
% %             x(ii) = 1/2 - u;
% %             y(ii) = sqrt(3)/2;
% %         
% %         % L5, lower triangular point
% %         case 5
% %             x(ii) = 1/2 - u;
% %             y(ii) = -sqrt(3)/2;
    end
end
    
% Return cooridinates of equilibrium points
L = [x,y,z];
end

