function [x_L, y_L] = CR3BP_newton_raph(x0, y0, v)
% This function finds the closest CR3BP Lagrange point for a set of intial
% conditions and v

z = 0;
y = y0; % Intial guess for y
x = x0; % Intial guess for x

err = 1;
i = 1;

while err > 1e-10
    
    r1 = sqrt((x(i) + v)^2 + y(i)^2 + z^2);
    r2 = sqrt((x(i) - 1 + v)^2 + y(i)^2 + z^2);
    
    f_x = x(i) - (1 - v)*(x(i) + v)/r1^3 - v*(x(i) - 1 + v)/r2^3;
    f_diff_x = (v - 1)/((v + x(i))^2 + y(i)^2 + z^2)^(3/2) - v/((v + x(i) - 1)^2 + y(i)^2 + z^2)^(3/2) + (3*v*(2*v + 2*x(i) - 2)*(v + x(i) - 1))/(2*((v + x(i) - 1)^2 + y(i)^2 + z^2)^(5/2)) - (3*(2*v + 2*x(i))*(v + x(i))*(v - 1))/(2*((v + x(i))^2 + y(i)^2 + z^2)^(5/2)) + 1;

    
    f_y = y(i) - (1 - v)*y(i)/r1^3 - v*y(i)/r2^3;
    f_diff_y = (v - 1)/((v + x(i))^2 + y(i)^2 + z^2)^(3/2) - v/((v + x(i) - 1)^2 + y(i)^2 + z^2)^(3/2) - (3*y(i)^2*(v - 1))/((v + x(i))^2 + y(i)^2 + z^2)^(5/2) + (3*v*y(i)^2)/((v + x(i) - 1)^2 + y(i)^2 + z^2)^(5/2) + 1;

    
    % x iteration
    x(i+1) = x(i) - f_x/f_diff_x;
    
    % y iteration
    y(i+1) = y(i) - f_y/f_diff_y;
    
    
    err_x = abs(x(i+1) - x(i));
    err_y = abs(y(i+1) - y(i));
    err = err_x + err_y;
    
    i = i + 1;
    
end

x_L = x(end);
y_L = y(end);

end

