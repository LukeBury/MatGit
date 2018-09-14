%%% Inputs (full states, [1x6])
% 1) Satellite wrt Europa (Europa-centered inertial)
% 2) Europa wrt Jupiter (JCI)
% 3) Jupiter wrt Sun (HCI)
% 4) Earth wrt Sun (HCI)
% 5) Station wrt Earth (ECI)
% 6) Rotation rate of Europa [1x1] (rad/s)
%%% Outputs
% 1) Htilde Matrix [2x6]
function [ Htilde ] = calculateHtilde_EurCI(sE,EJ,JSun,EaSun,stEa)
%%% Parsing States into components
xsE = sE(1); ysE = sE(2); zsE = sE(3); dxsE = sE(4); dysE = sE(5); dzsE = sE(6);
xEJ = EJ(1); yEJ = EJ(2); zEJ = EJ(3); dxEJ = EJ(4); dyEJ = EJ(5); dzEJ = EJ(6);
xJSun = JSun(1); yJSun = JSun(2); zJSun = JSun(3); dxJSun = JSun(4); dyJSun = JSun(5); dzJSun = JSun(6);
xEaSun = EaSun(1); yEaSun = EaSun(2); zEaSun = EaSun(3); dxEaSun = EaSun(4); dyEaSun = EaSun(5); dzEaSun = EaSun(6);
xstEa = stEa(1); ystEa = stEa(2); zstEa = stEa(3); dxstEa = stEa(4); dystEa = stEa(5); dzstEa = stEa(6);

%%% Assigning Htilde
Htilde = zeros(2,6);
Htilde(1,1) = (2*xEJ - 2*xEaSun + 2*xJSun + 2*xsE - 2*xstEa)/(2*((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2));
Htilde(1,2) = (2*yEJ - 2*yEaSun + 2*yJSun + 2*ysE - 2*ystEa)/(2*((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2));
Htilde(1,3) = (2*zEJ - 2*zEaSun + 2*zJSun + 2*zsE - 2*zstEa)/(2*((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2));
Htilde(2,1) = (dxEJ - dxEaSun + dxJSun + dxsE - dxstEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2) - (((dxEJ - dxEaSun + dxJSun + dxsE - dxstEa)*(xEJ - xEaSun + xJSun + xsE - xstEa) + (dyEJ - dyEaSun + dyJSun + dysE - dystEa)*(yEJ - yEaSun + yJSun + ysE - ystEa) + (dzEJ - dzEaSun + dzJSun + dzsE - dzstEa)*(zEJ - zEaSun + zJSun + zsE - zstEa))*(2*xEJ - 2*xEaSun + 2*xJSun + 2*xsE - 2*xstEa))/(2*((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(3/2));
Htilde(2,2) = (dyEJ - dyEaSun + dyJSun + dysE - dystEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2) - (((dxEJ - dxEaSun + dxJSun + dxsE - dxstEa)*(xEJ - xEaSun + xJSun + xsE - xstEa) + (dyEJ - dyEaSun + dyJSun + dysE - dystEa)*(yEJ - yEaSun + yJSun + ysE - ystEa) + (dzEJ - dzEaSun + dzJSun + dzsE - dzstEa)*(zEJ - zEaSun + zJSun + zsE - zstEa))*(2*yEJ - 2*yEaSun + 2*yJSun + 2*ysE - 2*ystEa))/(2*((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(3/2));
Htilde(2,3) = (dzEJ - dzEaSun + dzJSun + dzsE - dzstEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2) - (((dxEJ - dxEaSun + dxJSun + dxsE - dxstEa)*(xEJ - xEaSun + xJSun + xsE - xstEa) + (dyEJ - dyEaSun + dyJSun + dysE - dystEa)*(yEJ - yEaSun + yJSun + ysE - ystEa) + (dzEJ - dzEaSun + dzJSun + dzsE - dzstEa)*(zEJ - zEaSun + zJSun + zsE - zstEa))*(2*zEJ - 2*zEaSun + 2*zJSun + 2*zsE - 2*zstEa))/(2*((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(3/2));
Htilde(2,4) = (xEJ - xEaSun + xJSun + xsE - xstEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2);
Htilde(2,5) = (yEJ - yEaSun + yJSun + ysE - ystEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2);
Htilde(2,6) = (zEJ - zEaSun + zJSun + zsE - zstEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2);
% Htilde(2,1) = (dxEJ - dxEaSun + dxJSun + dxsE - dxstEa - wEur*ysE + wEur*(yEJ - yEaSun + yJSun + ysE - ystEa))/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2) - (((xEJ - xEaSun + xJSun + xsE - xstEa)*(dxEJ - dxEaSun + dxJSun + dxsE - dxstEa - wEur*ysE) + (yEJ - yEaSun + yJSun + ysE - ystEa)*(dyEJ - dyEaSun + dyJSun + dysE - dystEa + wEur*xsE) + (dzEJ - dzEaSun + dzJSun + dzsE - dzstEa)*(zEJ - zEaSun + zJSun + zsE - zstEa))*(2*xEJ - 2*xEaSun + 2*xJSun + 2*xsE - 2*xstEa))/(2*((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(3/2));
% Htilde(2,2) = (dyEJ - dyEaSun + dyJSun + dysE - dystEa - wEur*(xEJ - xEaSun + xJSun + xsE - xstEa) + wEur*xsE)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2) - (((xEJ - xEaSun + xJSun + xsE - xstEa)*(dxEJ - dxEaSun + dxJSun + dxsE - dxstEa - wEur*ysE) + (yEJ - yEaSun + yJSun + ysE - ystEa)*(dyEJ - dyEaSun + dyJSun + dysE - dystEa + wEur*xsE) + (dzEJ - dzEaSun + dzJSun + dzsE - dzstEa)*(zEJ - zEaSun + zJSun + zsE - zstEa))*(2*yEJ - 2*yEaSun + 2*yJSun + 2*ysE - 2*ystEa))/(2*((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(3/2));
% Htilde(2,3) = (dzEJ - dzEaSun + dzJSun + dzsE - dzstEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2) - (((xEJ - xEaSun + xJSun + xsE - xstEa)*(dxEJ - dxEaSun + dxJSun + dxsE - dxstEa - wEur*ysE) + (yEJ - yEaSun + yJSun + ysE - ystEa)*(dyEJ - dyEaSun + dyJSun + dysE - dystEa + wEur*xsE) + (dzEJ - dzEaSun + dzJSun + dzsE - dzstEa)*(zEJ - zEaSun + zJSun + zsE - zstEa))*(2*zEJ - 2*zEaSun + 2*zJSun + 2*zsE - 2*zstEa))/(2*((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(3/2));
% Htilde(2,4) = (xEJ - xEaSun + xJSun + xsE - xstEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2);
% Htilde(2,5) = (yEJ - yEaSun + yJSun + ysE - ystEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2);
% Htilde(2,6) = (zEJ - zEaSun + zJSun + zsE - zstEa)/((xEJ - xEaSun + xJSun + xsE - xstEa)^2 + (yEJ - yEaSun + yJSun + ysE - ystEa)^2 + (zEJ - zEaSun + zJSun + zsE - zstEa)^2)^(1/2);

end
