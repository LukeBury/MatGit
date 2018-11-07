function [impactAngle] = calcImpactAngle(r,v,unit)
%%% Calculates impact angle of spacecraft
%%% Inputs:
%          r - body-centric position vector [3x1]
%          v - velocity vector at impact [3x1]
%%% Outputs:
%          impactAngle - impact ange wrt surface (0-90 / 0-pi)
%%% =====================================================================
%%% Angle between velocity and surface
    A = r;
    B = v;
    impactAngle = acos(dot(A,B)/(norm(A)*norm(B)));
    if 0 <= impactAngle && impactAngle <= pi/2
        impactAngle = pi/2 - impactAngle;
    elseif pi/2 < impactAngle && impactAngle <= pi
        impactAngle = impactAngle - pi/2;
    end

    if isequal(unit, 'radians')
        return;
    elseif isequal(unit, 'degrees')
        impactAngle = impactAngle*180/pi;
    else
        warning('bad unit input')
    end

end
