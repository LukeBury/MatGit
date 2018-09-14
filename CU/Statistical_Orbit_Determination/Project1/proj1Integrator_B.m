% ------------------------------------------------------------------------
%%% Numerical Integrating Reference Orbit and STM
% ------------------------------------------------------------------------
function [ dY ] = proj1Integrator_B(t,Y,Pphi,c,JD0,Am)
    %%% Size dY to fit state and all of reshaped (n^2,1) STM
    dY = zeros(16+16^2,1);
    
    %%% Unpack state
    x = Y(1); % km
    y = Y(2); % km
    z = Y(3); % km
    dx = Y(4); % km/s
    dy = Y(5); % km/s
    dz = Y(6); % km/s
    Cr = Y(7);
    uE = Y(8);
    uS = Y(9);
    AU = Y(10);
    pb1 = Y(11);
    pb2 = Y(12);
    pb3 = Y(13);
    dpb1 = Y(14);
    dpb2 = Y(15);
    dpb3 = Y(16);

    %%% Reshape (n^2,1) stm to (n,n)
    stm = reshape(Y(17:end),16,16);
    
    %%% Describe Sun Location
    [rSun, vSun, mu_p] = Ephem(JD0 + t/86400,3,'EME2000'); % km, km/s, km^3/s^2
    rSun = -rSun;
    xs = rSun(1); % km
    ys = rSun(2); % km
    zs = rSun(3); % km

    %%% Build A matrix and evaluate at current state
    A = zeros(16,16);
    A(1:3,4:6) = eye(3,3);
    % Asym(4,1)
    A(4,1) = (3*uE*x*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - uS*(1/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (3*abs(x - xs)*sign(x - xs)*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2)) - uE/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) + (AU^2*Am*Cr*Pphi)/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (3*AU^2*Am*Cr*Pphi*abs(x - xs)*sign(x - xs)*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(4,2)
    A(4,2) = (3*uE*x*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(y - ys)*sign(y - ys)*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(y - ys)*sign(y - ys)*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(4,3)
    A(4,3) = (3*uE*x*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(z - zs)*sign(z - zs)*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(z - zs)*sign(z - zs)*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(4,7)
    A(4,7) = (AU^2*Am*Pphi*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));
    % Asym(4,8)
    A(4,8) = -x/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2);
    % Asym(4,9)
    A(4,9) = - (x - xs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - xs/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2);
    % Asym(4,10)
    A(4,10) = (2*AU*Am*Cr*Pphi*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));
    % Asym(5,1)
    A(5,1) = (3*uE*y*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(x - xs)*sign(x - xs)*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(x - xs)*sign(x - xs)*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(5,2)
    A(5,2) = (3*uE*y*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - uS*(1/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (3*abs(y - ys)*sign(y - ys)*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2)) - uE/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) + (AU^2*Am*Cr*Pphi)/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (3*AU^2*Am*Cr*Pphi*abs(y - ys)*sign(y - ys)*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(5,3)
    A(5,3) = (3*uE*y*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(z - zs)*sign(z - zs)*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(z - zs)*sign(z - zs)*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(5,7)
    A(5,7) = (AU^2*Am*Pphi*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));
    % Asym(5,8)
    A(5,8) = -y/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2);
    % Asym(5,9)
    A(5,9) = - (y - ys)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - ys/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2);
    % Asym(5,10)
    A(5,10) = (2*AU*Am*Cr*Pphi*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));
    % Asym(6,1)
    A(6,1) = (3*uE*z*abs(x)*sign(x))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(x - xs)*sign(x - xs)*(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(x - xs)*sign(x - xs)*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(6,2)
    A(6,2) = (3*uE*z*abs(y)*sign(y))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) + (3*uS*abs(y - ys)*sign(y - ys)*(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2) - (3*AU^2*Am*Cr*Pphi*abs(y - ys)*sign(y - ys)*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(6,3)
    A(6,3) = (3*uE*z*abs(z)*sign(z))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(5/2) - uS*(1/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (3*abs(z - zs)*sign(z - zs)*(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2)) - uE/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) + (AU^2*Am*Cr*Pphi)/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (3*AU^2*Am*Cr*Pphi*abs(z - zs)*sign(z - zs)*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2));
    % Asym(6,7)
    A(6,7) = (AU^2*Am*Pphi*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));
    % Asym(6,8)
    A(6,8) = -z/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2);
    % Asym(6,9)
    A(6,9) = - (z - zs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - zs/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2);
    % Asym(6,10)
    A(6,10) = (2*AU*Am*Cr*Pphi*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2));

    %%% Calculate new STM
    stm_dot = A*stm;
        
    %%% Creating X-dot
    % Spacecraft velocities
    dY(1:3) = [dx; dy; dz];
    %%% Using u and SRP as dynamics
    % EQM(4)
    dY(4) = (AU^2*Am*Cr*Pphi*(x - xs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (uE*x)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - uS*((x - xs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) + xs/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2));
    % EQM(5)
    dY(5) = (AU^2*Am*Cr*Pphi*(y - ys))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (uE*y)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - uS*((y - ys)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) + ys/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2));
    % EQM(6)
    dY(6) = (AU^2*Am*Cr*Pphi*(z - zs))/(c*(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2)) - (uE*z)/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(3/2) - uS*((z - zs)/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) + zs/(abs(xs)^2 + abs(ys)^2 + abs(zs)^2)^(3/2));
    %%% Cr has no dynamics
    dY(7:16) = zeros(10,1);

    % Filling in reshaped (7^2,1) STM to state
    dY(17:end) = reshape(stm_dot,256,1);
end