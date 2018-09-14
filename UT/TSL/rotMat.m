
function [mat] = rotMat(ang,type)



if(type==1)
    mat = [1 0 0;
           0 cosd(ang) sind(ang);
           0 -sind(ang) cosd(ang)];
elseif(type==2)
    mat = [cosd(ang) 0 -sind(ang);
       0 1 0;
       sind(ang) 0 cosd(ang)];
else
    mat = [cosd(ang) sind(ang) 0;
       -sind(ang) cosd(ang) 0;
       0 0 1];
end
    

%declare as a double (ang) and an int (type)