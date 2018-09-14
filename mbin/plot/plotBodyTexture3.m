function plotBodyTexture3( radius, position , img)
%%% Plot 3D spherical body
%%% Inputs
%         radius - radius of body
%         position - 3D position vector for body center
%         img - imread of image file ... imread('picture.jpg')

[x,y,z] = sphere;
warp(x*radius+position(1) ,y*radius+position(2), z*radius+position(3),flipud(img))
set(gca,'YDir','normal')
grid on
end

