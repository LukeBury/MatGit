clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
% ------------------------------------------------------------------------
%%% 3.11
% ------------------------------------------------------------------------
% r1 = 10;
% r2 = -15;
% r3 = 20;
% 
% dr1 = 2;
% dr2 = 1;
% dr3 = 0;
% dr = [dr1; dr2; dr3];
% 
% B_inv = [-sind(r2) 0 1;...
%     cosd(r2)*sind(r3), cosd(r3), 0;...
%     cosd(r2)*cosd(r3), -sind(r3), 0];
% 
% w_B = B_inv * dr
% 
% E321 = [cosd(r1)*cosd(r2), cosd(r2)*sind(r1), -sind(r2);...
%     sind(r3)*sind(r2)*cosd(r1)-cosd(r3)*sind(r1), sind(r3)*sind(r2)*sind(r1)+cosd(r3)*cosd(r1), sind(r3)*cosd(r2);...
%     cosd(r3)*sind(r2)*cosd(r1)+sind(r3)*sind(r1), cosd(r3)*sind(r2)*sind(r1)-sind(r3)*cosd(r1), cosd(r3)*cosd(r2)];
% 
% w_N = (E321')*w_B


% ------------------------------------------------------------------------
%%% 3.12
% ------------------------------------------------------------------------

t0 = 0; dt = .001; tf = 60;
time = dt:dt:tf;
deg2rad = pi/180;
angs = zeros(3,length(time)+1);
angs(:,1) = [40; 30; 80].*deg2rad; % rad

for ti = 1:length(time)
    t = time(ti);
    w = [sin(.1*t); 0.01; cos(.1*t)].*(20*deg2rad);
    
    r1 = angs(1,ti);
    r2 = angs(2,ti);
    r3 = angs(3,ti);
    mat = (1/cos(r2)).*[0, sin(r3), cos(r3);...
        0, cos(r3)*cos(r2), -sin(r3)*cos(r2);...
        cos(r2), sin(r3)*sin(r2), cos(r3)*sin(r2)];
    rate = mat*w;
    angs(:,ti+1) = angs(:,ti) + rate.*dt;
end
% Converting back to degrees
angs = angs./deg2rad;

% Plotting
colors = get_color_palettes();
figure; %hold all
subplot(2,1,1); hold all
title('Calculated Values')
plot(time,angs(1,1:end-1),'linewidth',2,'color',colors.sch.d3_1(1,:),'markersize',18);
plot(time,angs(2,1:end-1),'linewidth',2,'color',colors.sch.d3_1(2,:),'markersize',18);
plot(time,angs(3,1:end-1),'linewidth',2,'color',colors.sch.d3_1(3,:),'markersize',18);
legend('yaw','pitch','roll')
PlotBoi2('Time, s','Angle, °',14)

for t = 1:length(time)
    % Adjusting Yaw (-180 180)
    if angs(1,t+1) > 180
        angs(1,t+1) = angs(1,t+1) - 360;
    elseif angs(1,t+1) < -180
        angs(1,t+1) = angs(1,t+1) + 360;
    end
    
    % Adjusting Pitch (-90 90)
    if angs(2,t+1) > 180
        angs(2,t+1) = angs(2,t+1) - 180;
    elseif angs(1,t+1) < -180
        angs(2,t+1) = angs(2,t+1) + 180;
    end
    
    % Adjusting Roll (-180 180)
    if angs(3,t+1) > 180
        angs(3,t+1) = angs(3,t+1) - 360;
    elseif angs(1,t+1) < -180
        angs(3,t+1) = angs(3,t+1) + 360;
    end
end

subplot(2,1,2); hold all
title('Adjusted Values')
plot(time,angs(1,1:end-1),'.','linewidth',1,'color',colors.sch.d3_1(1,:),'markersize',6);
plot(time,angs(2,1:end-1),'.','linewidth',1,'color',colors.sch.d3_1(2,:),'markersize',6);
plot(time,angs(3,1:end-1),'.','linewidth',1,'color',colors.sch.d3_1(3,:),'markersize',6);
legend('yaw','pitch','roll')
PlotBoi2('Time, s','Angle, °',14)
ylim([-180 180])
yticks([-180:45:180])

% ------------------------------------------------------------------------
%%% 3.13
% ------------------------------------------------------------------------
% 
% deg2rad = pi/180;
% r1 = -30*deg2rad;
% r2 = 40*deg2rad;
% r3 = 20*deg2rad;
% 
% % 3-1-3 DCM
% C = [cos(r3)*cos(r1) - sin(r3)*cos(r2)*sin(r1), cos(r3)*sin(r1)+sin(r3)*cos(r2)*cos(r1), sin(r3)*sin(r2);...
%     -sin(r3)*cos(r1)-cos(r3)*cos(r2)*sin(r1), -sin(r3)*sin(r1)+cos(r3)*cos(r2)*cos(r1), cos(r3)*sin(r2);...
%     sin(r2)*sin(r1), -sin(r2)*cos(r1), cos(r2)];
% phi = acos(.5*(C(1,1) + C(2,2) + C(3,3) - 1));
% phiP = phi - 2*pi;
% ehat = (1/(2*sin(phi))) * [C(2,3)-C(3,2); C(3,1)-C(1,3); C(1,2)-C(2,1)];
% 
% % Euler Parameters (quaternion set)
% B0 = cos(phi/2);
% B1 = ehat(1)*sin(phi/2);
% B2 = ehat(2)*sin(phi/2);
% B3 = ehat(3)*sin(phi/2);
% % CRP
% q1 = B1/B0;
% q2 = B2/B0;
% q3 = B3/B0;
% % MRP
% s1 = B1/(1+B0)
% s2 = B2/(1+B0)
% s3 = B3/(1+B0)

% ------------------------------------------------------------------------
%%% 3.14
% ------------------------------------------------------------------------

% r = 45*pi/180;
% e = 1/sqrt(3);
% g = 1-cos(r);
% 
% ee = e*e*g;
% C = [ee+cos(r), ee+e*sin(r), ee-e*sin(r);...
%     ee-e*sin(r), ee+cos(r), ee+e*sin(r);...
%     ee+e*sin(r), ee-e*sin(r), ee+cos(r)]
% r1 = atan(C(1,2)/C(1,1))
% atan2(C(1,2),C(1,1))
% r2 = -asin(C(1,3))
% r3 = atan(C(2,3)/C(3,3))
% atan2(C(2,3),C(3,3))
% 
% r1*180/pi
% r2*180/pi
% r3*180/pi


% ------------------------------------------------------------------------
%%% 3.15
% ------------------------------------------------------------------------
% e = [1/sqrt(3);1/sqrt(3);1/sqrt(3)];
% etilde = [0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
% phi = 45*pi/180;
% expm(-phi*etilde)
% eye(3,3).*cos(phi) - sin(phi).*etilde + (1-cos(phi))*e*e'

% ------------------------------------------------------------------------
%%% 3.18
% ------------------------------------------------------------------------
% syms B0a B1a B2a B3a B0b B1b B2b B3b real
% BN = [B0a^2+B1a^2-B2a^2-B3a^2,   2*(B1a*B2a+B0a*B3a),       2*(B1a*B3a-B0a*B2a);...
%        2*(B1a*B2a-B0a*B3a),   B0a^2-B1a^2+B2a^2-B3a^2,     2*(B2a*B3a+B0a*B1a);...
%        2*(B1a*B3a+B0a*B2a),     2*(B2a*B3a-B0a*B1a),    B0a^2-B1a^2-B2a^2+B3a^2];
% 
% FB = [B0b^2+B1b^2-B2b^2-B3b^2,   2*(B1b*B2b+B0b*B3b),       2*(B1b*B3b-B0b*B2b);...
%        2*(B1b*B2b-B0b*B3b),   B0b^2-B1b^2+B2b^2-B3b^2,     2*(B2b*B3b+B0b*B1b);...
%        2*(B1b*B3b+B0b*B2b),     2*(B2b*B3b-B0b*B1b),    B0b^2-B1b^2-B2b^2+B3b^2];
% 
% 
% FBBN = FB*BN;
% 
% B = [B0b -B1b -B2b -B3b;...
%      B1b  B0b  B3b -B2b;...
%      B2b -B3b  B0b  B1b;...
%      B3b  B2b -B1b  B0b];
% B = B*[B0a; B1a; B2a; B3a];
% B0 = B(1); B1 = B(2); B2 = B(3); B3 = B(4);
% 
% FN = [B0^2+B1^2-B2^2-B3^2,   2*(B1*B2+B0*B3),       2*(B1*B3-B0*B2);...
%        2*(B1*B2-B0*B3),   B0^2-B1^2+B2^2-B3^2,     2*(B2*B3+B0*B1);...
%        2*(B1*B3+B0*B2),     2*(B2*B3-B0*B1),    B0^2-B1^2-B2^2+B3^2];
% 
% 
% % B0a = sqrt(.15); B0b = sqrt(.1);
% % B1a = sqrt(.25); B1b = sqrt(.2);
% % B2a = sqrt(.32); B2b = sqrt(.3);
% % B3a = sqrt(.28); B3b = sqrt(.4);
% % 
% % 
% % vpa(subs(FBBN),2)
% % vpa(subs(FN),2)

% ------------------------------------------------------------------------
%%% 3.19
% ------------------------------------------------------------------------
% syms b0 b1 b2 b3 w1 w2 w3
% 
% db1 = (2*w3*(b0*b2-b1*b3) + w1*(b0^2-b1^2-b2^2+b3^2) + w1*(b0^2-b1^2+b2^2-b3^2) - 2*w2*(b1*b2+b0*b3))/(4*b0) - (2*(b2*b3+b0*b1)-2*(b2*b3-b0*b1))*(2*(-b1*w1-b2*w2-b3*w3))/(16*b0^2);
% db1
% simplify(db1)
% collect(db1, [b0, b1, b2, b3])















