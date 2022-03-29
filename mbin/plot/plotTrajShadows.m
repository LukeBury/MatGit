function plotTrajShadows(Xin, lw, color, varargin)
%%% Description
% Plot shadows of a trajectory to help show 3D shape. When calling this
% function, provide a trajectory (Xin), a linewidth (lw), and a color for
% the shadows (color). After that, separated by commas, list any desired
% axes for a shadow along with a position along that axis for the shadow to
% appear. A full call to this function may look like:
%   plotTrajShadows(X_traj, lw, color, 'x', 1.01, 'y', 0.02, 'z', 0.025)
%       
% ------------------------------------------------------------------------
%%% Inputs
% Xin      - [nx6 or nx3] state history of trajectory
% lw       - [scalar] linewidth for shadows
% color    - [1x3] color for shadows. Grey = [1, 1, 1].*(127/255)
% varargin - {nx1} 
%             -(optional) specify axis and position for shadow (separated 
%              by commas)
%              ex: varargin(1) = 'x', varargin(2) = 1.01
%             -(optional) Entering 'bodyshadow' as an additional argument
%              will turn on projected shadow of a body of a given radius
%              and x-axis position
%              ex: varargin(1) = 'bodyshadow', varargin(2) = [x-pos, radius]
%
% ------------------------------------------------------------------------
%%% Outputs
% N/A
% ------------------------------------------------------------------------
% Created: 1/26/21
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Loop through variable-length input to determine if the body shadow is
%%% desired
bodyshadow_flag = false;
for arg_index = 1:length(varargin)
    if isequal(lower(varargin{arg_index}), 'bodyshadow')
        bodyshadow_flag = true;
        
        body_xpos     = varargin{arg_index+1}(1);
        radius = varargin{arg_index+1}(2);
        
        circle_vertices = nsidedpoly(1000, 'Center', [0, 0], 'Radius', radius);
        break
    end
end

%%% Loop through the variable-length input argument list and look for axis
%%% shadow specfications
for arg_index = 1:length(varargin)
    
    %%% If the current argument is 'x', plot the x-axis shadow
    if isequal(varargin{arg_index}, 'x')
        x_location = varargin{arg_index+1};
        plot3(ones(size(Xin,1),1).*x_location, Xin(:,2), Xin(:,3),'linewidth',lw,'color',color)
        
        if bodyshadow_flag
            fill3(ones(size(circle_vertices.Vertices,1),1)*x_location, circle_vertices.Vertices(:,1), circle_vertices.Vertices(:,2), [0.4980, 0.4980, 0.4980], 'Facealpha', 0.4, 'Edgealpha', 0)
        end
    %%% If the current argument is 'y', plot the y-axis shadow
    elseif isequal(varargin{arg_index}, 'y')
        y_location = varargin{arg_index+1};
        plot3(Xin(:,1), ones(size(Xin,1),1).*y_location, Xin(:,3),'linewidth',lw,'color',color)
        
        if bodyshadow_flag
            fill3(circle_vertices.Vertices(:,1) + body_xpos, ones(size(circle_vertices.Vertices,1),1)*y_location, circle_vertices.Vertices(:,2), [0.4980, 0.4980, 0.4980], 'Facealpha', 0.4, 'Edgealpha', 0)
        end
    %%% If the current argument is 'zy', plot the z-axis shadow
    elseif isequal(varargin{arg_index}, 'z')
        z_location = varargin{arg_index+1};
        plot3(Xin(:,1), Xin(:,2), ones(size(Xin,1),1).*z_location,'linewidth',lw,'color',color)
        
        if bodyshadow_flag
            fill3(circle_vertices.Vertices(:,1) + body_xpos, circle_vertices.Vertices(:,2), ones(size(circle_vertices.Vertices,1),1)*z_location, [0.4980, 0.4980, 0.4980], 'Facealpha', 0.4, 'Edgealpha', 0)
        end
    end
end

end % function






% 
% 



