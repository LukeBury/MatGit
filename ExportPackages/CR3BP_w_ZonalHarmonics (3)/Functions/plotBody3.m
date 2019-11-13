function plotBody3( radius, position, color, transp)
%%% Plot 3D spherical body
%%% Inputs
%         radius - radius of body
%         position - 3D position vector for body center
%         color - [1x3] color code for body
%         transp - transparency of object (0 = invisible, 1 = full color)

[x,y,z] = sphere;
surf(x*radius+position(1) ,y*radius+position(2), z*radius+position(3));
colormap(color) 
if nargin == 4
    alpha(transp)
end
view(3);
grid on

end
% rotate(h,direction,alpha)