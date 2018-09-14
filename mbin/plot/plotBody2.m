function plotBody2( radius, position, facecolor, edgecolor, lw , alpha)
%%% Plot 2D Body
%%% Inputs
%         radius - radius of body
%         position - 2D position vector for body center
%         facecolor - [1x3] color code for body fill **** [[1x3],0] for transparent!
%         edgecolor - [1x3] color code for body edge
%         linewidth - desired width of line
%         alpha - transparency of object (0 = invisible, 1 = full color)

if nargin == 5
    rectangle('Position',[position(1)-radius,position(2)-radius,2*radius,2*radius],'Curvature',[1,1], 'FaceColor',facecolor,'EdgeColor',edgecolor,'linewidth',lw);
elseif nargin == 6
    rectangle('Position',[position(1)-radius,position(2)-radius,2*radius,2*radius],'Curvature',[1,1], 'FaceColor',[facecolor, alpha],'EdgeColor',edgecolor,'linewidth',lw);
end
end

