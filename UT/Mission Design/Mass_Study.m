clear
clc
close all

tof = [11.45:.01:16]'; % Time of Flight, Days
%%%------------------------------------------------------------------------
%%%                        Battery Variables
%%%------------------------------------------------------------------------

Pidle = 1.13; %W

Li_ion_Dens = 125; %W-hr/kg
LiSO2_Dens  = 350; %W-hr/kg
LiS_Dens    = 500; %W-hr/kg

% Bmass = tof.*(Pidle*24)./Li_ion_Dens;      %kg - Li-Ion
% Bmass = tof.*(Pidle*24)./LiSO2_Dens;         %kg - LISO2
Bmass = tof.*(Pidle*24)./LiS_Dens;

%%%------------------------------------------------------------------------
%%%                        Fuel Variables
%%%------------------------------------------------------------------------

% From tof vs dv fit:
A = -35.54;
B = 35.66;

Ve = 9.81*300/1000;

mf = exp((tof-B)/(A*Ve));%% Mass Fraction - only works for flights on fit!!!
% m20 = mf.*(20+Bmass)-20-Bmass;
m20 = mf.*(30+Bmass)-30-Bmass; %really m30 b/c 30 is much more accurate than 20 (final = 32.7)
%%%------------------------------------------------------------------------
%%%                        Solution Specific
%%%------------------------------------------------------------------------
mf14 = exp((14.37-B)/(A*Ve)); %mass fraction for tof=14.37 days

Bmass14li = 14.37*(Pidle*24)./Li_ion_Dens;     %kg - Li-Ion
Bmass14ls = 14.37*(Pidle*24)./LiSO2_Dens;      %kg - LISO2
Bmass14lis= 14.37*(Pidle*24)./LiS_Dens;        %kg - LiS

%Vector of all reasonable *dry* weights at tof=14.37 solution
mi = [12:.01:36]'; %guess of dry mass w/ no battery (kg)
n=size(mi);

mfuel14li=zeros(size(mi));
for i=1:n
    mfuel14li(i) = (mi(i)+Bmass14li)*mf14-mi(i)-Bmass14li;
end
% mfuel14(81)

mfuel14ls=zeros(size(mi));
for i=1:n
    mfuel14ls(i) = (mi(i)+Bmass14ls)*mf14-mi(i)-Bmass14ls;
end

mfuel14lis=zeros(size(mi));
for i=1:n
    mfuel14lis(i) = (mi(i)+Bmass14lis)*mf14-mi(i)-Bmass14lis;
end

%vector of all total masses based on mi range above for tof=14.37 soltn
mwet14li = Bmass14li + mfuel14li;
mwet14ls = Bmass14ls + mfuel14ls;
mwet14lis= Bmass14lis+ mfuel14lis;

%%%------------------------------------------------------------------------
%%%                        Summation
%%%------------------------------------------------------------------------

mtot = Bmass+m20;

p=polyfit(tof,mtot,1);
% p =  -0.0207    7.9283

%%%------------------------------------------------------------------------
%%%                        Plotting
%%%------------------------------------------------------------------------
figure(2) %%%Proves mass/day is much higher for fuel, so we should minimize
          %%%Fuel. This is why 14.37 has a dotted line, it's the best
          %%%Solution for this
          %%% Currently, Li-Ion is the battery type, but the trend doesn't
          %%% change because Li-Ion is the heaviest battery we're
          %%% considering. If it's better than dV, than the others are too
hold all
plot(tof,m20,'b','linewidth',3);
plot(tof,Bmass,'r','linewidth',3);
plot(tof,mtot,'m','linewidth',3);
% plot(14.37,mtot14,'x','markersize',15)
% plot([14.37 14.37],get(gca,'ylim'),'--k','linewidth',1)
grid on
PlotBoi('TOF, days','Mass, kg')
legend('Propellent','Battery','Battery+Propellent','ToF = 14.37 Days')


figure(3) %%%At tof=14.37 solution, input *dry* mass vs total s/c mass, LiIon
hold on
plot(mi,mwet14li,'m','linewidth',3)
plot(mi,mwet14li+mi,'k','linewidth',3)
plot([10:.1:40],[41:.001/300:41.001],'--r','linewidth',1)
plot([26.25 26.25],get(gca,'ylim'),'--r','linewidth',1)
grid on
PlotBoi('Input Dry* Mass, kg','Total Mass, kg')
legend('Fuel + Li-ion Battery','SKY Total Mass','SKY Mass Limit')
%max *dry* mass for Li-Ion battery SKY = 26.25 kg


figure(4) %%%At tof=14.37 solution, input *dry* mass vs total s/c mass, LiS02
hold on
plot(mi,mwet14ls,'m','linewidth',3)
plot(mi,mwet14ls+mi,'k','linewidth',3)
plot([10:.1:40],[41:.001/300:41.001],'--r','linewidth',1)
plot([28.26 28.26],get(gca,'ylim'),'--r','linewidth',1)
grid on
PlotBoi('Input Dry* Mass, kg','Total Mass, kg')
legend('Fuel + LiSO2 Battery','SKY Total Mass','SKY Mass Limit')
%max *dry* mass for LiSO2 battery SKY = 28.26 kg

figure(5) %%%At tof=14.37 solution, input *dry* mass vs total s/c mass, LiS
hold on
plot(mi,mwet14lis,'m','linewidth',3)
plot(mi,mwet14lis+mi,'k','linewidth',3)
plot([10:.1:40],[41:.001/300:41.001],'--r','linewidth',1)
plot([28.59 28.59],get(gca,'ylim'),'--r','linewidth',1)
grid on
PlotBoi('Input Dry* Mass, kg','Total Mass, kg')
legend('Fuel + LiS Battery','SKY Total Mass','SKY Mass Limit')
%max *dry* mass for LiS battery SKY = 28.59 kg

%%
figure (6)
hold all
plot(mi,mwet14li+mi,'linewidth',2)
plot(mi,mwet14ls+mi,'linewidth',2)
plot(mi,mwet14lis+mi,'linewidth',2)
plot([10:.1:40],[41:.001/300:41.001],'--r','linewidth',1)
plot([32.66 32.66],get(gca,'ylim'),'--r','linewidth',1)
grid on
PlotBoi('Input Dry* Mass, kg','Total Mass, kg')
legend('Li-Ion','LiSO2','LiS')



















