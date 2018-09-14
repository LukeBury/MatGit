
%%% Workspace
% Pvo       - Pressure Transducer Pressure Output
% Pstd      - Standard Deviation of ^^^
% Svo       - Strain Gage Output
% Svostd    - Standard Deviation of ^^^

%%% Build Measurement Numbers;
x = [1:13];
%%% Poisson's Ratios
vpoly = .37;

%%% Young's Moduli
Epoly = 354e3;

%%% Gage Factor
G = 2.11; % +- .5%
dG = .01055;
%%% Polycarbonate Tube Dimensions
ODpoly = 4;             % Outer Diameter
dODpoly = .05;
IDpoly = 3.75;          % Inner Diameter
dODpoly = .05;
rpoly = IDpoly/2;       % Radius
drpoly = .025;
tpoly = (ODpoly-IDpoly)/2;  % Thickness
dtpoly = .1;

%%%%%%%%%%%%% Polycarbonate - Pressure based on strain
Svr = (Svo-Svo(1))./(100*5);
dSvr = Svostd./500;
Ehs = -2 .* Svr ./ (G*((vpoly+1)-Svr.*(vpoly-1))); % hoop strain
for i=[1:13]
dEhs(i) = sqrt((2*Svr(i)*dG/(G^2*(vpoly*(-Svr(i))+vpoly+Svr(i)+1)))^2 + (2*(vpoly+1)*dSvr(i)/(G*(vpoly*(-Svr(i))+vpoly+Svr(i)+1)^2))^2);
end
dEhs; % error on hoop strain

Ppoly = Ehs.*Epoly*tpoly/rpoly; % pressure from strain
for k=[1:13]
dPpoly(k) = sqrt((Ehs(k)*tpoly*dEhs(k)/rpoly)^2+(Ehs(k)*Epoly*dtpoly/rpoly)^2+(Ehs(k)*tpoly*Epoly*drpoly/(rpoly^2))^2);
end
dPpoly;

%%% Polycarbonate - Pressure Transducer Calculations
P = (Pvo-Pvo(1))*1000;

expectedP = [0,10,20,30,40,50,60,50,40,30,20,10,0];
hold all
plot(x,Ppoly,'b','linewidth',2)
plot(x,P,'g','linewidth',2)
plot(x,expectedP,'k','linewidth',2)
legend('Strain Gage','Pressure Transducer','Expected')
xlabel('Measurement Number','FontName','Cambria','Fontsize',18)
ylabel('Pressure (psig)','FontName','Cambria','Fontsize',18)
errorbar(x,P,Pstd, 'rx','linewidth',2)
errorbar(x,Ppoly,dPpoly, 'rx','linewidth',2)
%%% Polycarbonate - Strain from Gages
Ehs;
%%% Polycarbonate - Strain based on Pressure transducer
Ehp = rpoly .* P ./ (Epoly * tpoly);

% figure
% hold all
% plot(x,Ehs)
% plot(x,Ehp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%       Pop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

valum = .33;
Ealum = 1e7;
Od = 2.579; %in      % outder diameter non pressurized
tc = .00709; %in      % thickness of can
dtc = 3.94e-4;
rc = (Od-2*tc)/2;   % inner radius of can
drc = (3.94e-4)/2+.00709;


Before = -0.028286828;
After  = -4.493549706;

Vrb = 0;
Vra = (After-Before)/(500*5);
dVra = 0.0034265/(500*5);


Eb = 0;
Ea = -4*Vra/(G*(1+2*Vra));
dEa = sqrt((4*dVra/(G*(2*Vra+1)^2))^2+(4*Vra*dG/(G^2*(2*Vra+1)))^2);

Pcan = Ea*Ealum*tc/(rc * (1 - valum/2))
dPcan = sqrt((dEa*Ealum*tc/(rc*(1-valum/2)))^2+(dtc*Ea*Ealum/(rc*(1-valum/2)))^2+(drc*Ea*Ealum*tc/(rc^2*(valum-2)))^2);












