clear
clc
close all

Rj = 71492;     %km
uj = 126686534; %km^3/s^2
a1 = 4282371;   %km
a2 = 4281985.3; %km
e1 = 0.7997;
e2 = 0.7998;

E = [.0042:-.00005:0]; %Eccentric Anomalies of ESC
E = E;
f = 2*atan(sqrt((1+e1)/(1-e1)).*tan(E./2)); %True Anomalies of ESC

r1_mag = a1.*(1-e1.*cos(E));

r1 = zeros(3,85);
vi = zeros(1,85);
% for i = 1:85
%     r1(1,i) = r1_mag(i).*sin(f(i));
%     r1(2,i) = r1_mag(i).*cos(f(i));
% end
for i = 1:85
    r1(1,i) = r1_mag(i).*sin(f(i)); %r1x
    r1(2,i) = r1_mag(i).*cos(f(i)); %r1y
    
    vi(i) = sqrt(uj*(2/norm(r1(:,i))-1/a1)); %velocity magnitude
end

% r2 = [0;856482.457;0]; % vectory to ganymede impact site - km
r2 = [0;8.5776e+05;0]; % vectory to ganymede impact site - km
dT = [72999:-60:1]; %Possible dTs



% lambert(r1vec,r2vec,tf,m,GRADE,muC), returns
[v1,v2,ef]= lambert(r1(:,1),r2,dT(1),0,0,uj)

% dV = zeros(84,1217);
% v = zeros(84,1217,3);
% for r = 1:84     % various r1 position vectors
%     for t = 1:1217 % various ToF values
%         [v1,v2,ef]= lambert(r1(:,r),r2,dT(t),0,0,uj);
%         
%         Vvec = vi(r)*((r1(:,r+1)-r1(:,r))/norm(r1(:,r+1)-r1(:,r))); 
%         dV(r,t) = norm(v1-Vvec); 
%         v(r,t,:) = v1;
%         
%     end
%         
% end


% for i =1:85
%     hold all
%     plot3(r2(1),r2(2),r2(3),'x','markersize',10,'linewidth',1)
%     plot3(r1(1,i),r1(2,i),r1(3,i),'o','markersize',10,'linewidth',1)
%     drawnow
% end


