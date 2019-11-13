function plotBodyTexture3( radius, position , img, alph)
%%% Plot 3D spherical body
%%% Inputs
%         radius - radius of body
%         position - 3D position vector for body center
%         img - imread of image file ... imread('picture.jpg')
%         alph - alpha value for body (1 = solid, 0 = invisible)
[x,y,z] = sphere(100);
warp(x*radius+position(1) ,y*radius+position(2), z*radius+position(3),flipud(img))
set(gca,'YDir','normal')
grid on

if nargin == 4
    alpha(alph)
end
end

