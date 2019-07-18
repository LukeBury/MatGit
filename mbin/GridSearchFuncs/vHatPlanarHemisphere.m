function [vHats] = vHatPlanarHemisphere(n_points,centerAxis,plane)
%%% Description
%       Returns a set of vHat vectors that populate a planar hemisphere
%       
% --------------------------------------------------------------
%%% Inputs
%       n_points   - [1x1] desired number of vHat vectors 
%       centerAxis - Axis for planar vHats to be centered on ('x','-x','y','-y','z', or '-z')
%       plane      - Desired plane, ('xy','yz', or 'xz')
% --------------------------------------------------------------
%%% Outputs
%       vHats - [nx3] vHat vectors in hemisphere 
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% 
% ------------------------------------
%%% Preallocating
vHats = zeros(n_points,3);

%%% Defining rotation angles
rotAngles_rad = linspace(-pi/2,pi/2,n_points);

%%% Defining basis vectors
x_vec = [1; 0; 0];
y_vec = [0; 1; 0];
z_vec = [0; 0; 1];

%%% Populating planar hemisphere
if     isequal(centerAxis,'x') && isequal(plane,'xy')
    for kk = 1:n_points
        vHats(kk,:) = R3(x_vec, rotAngles_rad(kk));
    end
elseif isequal(centerAxis,'x') && isequal(plane,'xz')
    for kk = 1:n_points
        vHats(kk,:) = R2(x_vec, rotAngles_rad(kk));
    end
elseif isequal(centerAxis,'-x') && isequal(plane,'xy')
    for kk = 1:n_points
        vHats(kk,:) = R3(-x_vec, rotAngles_rad(kk));
    end
elseif isequal(centerAxis,'-x') && isequal(plane,'xz')
    for kk = 1:n_points
        vHats(kk,:) = R2(-x_vec, rotAngles_rad(kk));
    end
    
    
elseif isequal(centerAxis,'y') && isequal(plane,'xy')
    for kk = 1:n_points
        vHats(kk,:) = R3(y_vec, rotAngles_rad(kk));
    end
elseif isequal(centerAxis,'y') && isequal(plane,'yz')
    for kk = 1:n_points
        vHats(kk,:) = R1(y_vec, rotAngles_rad(kk));
    end
elseif isequal(centerAxis,'-y') && isequal(plane,'xy')
    for kk = 1:n_points
        vHats(kk,:) = R3(-y_vec, rotAngles_rad(kk));
    end
elseif isequal(centerAxis,'-y') && isequal(plane,'yz')
    for kk = 1:n_points
        vHats(kk,:) = R1(-y_vec, rotAngles_rad(kk));
    end
    
    
elseif isequal(centerAxis,'z') && isequal(plane,'xz')
    for kk = 1:n_points
        vHats(kk,:) = R2(z_vec, rotAngles_rad(kk));
    end
elseif isequal(centerAxis,'z') && isequal(plane,'yz')
    for kk = 1:n_points
        vHats(kk,:) = R1(z_vec, rotAngles_rad(kk));
    end
elseif isequal(centerAxis,'-z') && isequal(plane,'xz')
    for kk = 1:n_points
        vHats(kk,:) = R2(-z_vec, rotAngles_rad(kk));
    end
elseif isequal(centerAxis,'-z') && isequal(plane,'yz')
    for kk = 1:n_points
        vHats(kk,:) = R1(-z_vec, rotAngles_rad(kk));
    end
else
    warning('Input error')
end

end
