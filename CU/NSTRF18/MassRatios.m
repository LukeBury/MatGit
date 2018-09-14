clear
clc
close all
addpath('/Users/lukebury/Documents/MATLAB/mbin')
addpath('/Users/lukebury/Documents/MATLAB/CU/bin/LowEnergy')

% %%% Get plot colors
% colors = get_color_palettes();

% Values from Wikipedia
%%% ----------------
earth.mass = 5.972e24;
moon.mass = 7.347e22;

%%% ----------------
mars.mass = 6.39e23;
phobos.mass = 1.0659e16;
deimos.mass = 1.476e15;

%%% ----------------
jupiter.mass = 1.898e27;
io.mass = 8.94e22;
europa.mass = 4.799e22;
ganymede.mass = 1.48e23;
callisto.mass = 1.08e23;

%%% ----------------
saturn.mass = 5.683e26;
enceladus.mass = 1.08e20;
mimas.mass = 3.8e19;
rhea.mass = 2.306e21;
titan.mass = 1.345e23;

%%% ----------------
uranus.mass = 8.681e25;
titania.mass = 3.42e21;
oberon.mass = 2.883e21;

%%% ----------------
neptune.mass = 1.024e26;
triton.mass = 2.14e22;

%%% ----------------
pluto.mass = 1.309e22;
charon.mass = 1.586e21;

%%% ----------------%%% ----------------%%% ----------------
e.m = moon.mass / (earth.mass + moon.mass);

m.p = phobos.mass / (mars.mass + phobos.mass);
m.d = deimos.mass / (mars.mass + deimos.mass);

j.i = io.mass / (jupiter.mass + io.mass);
j.e = europa.mass / (jupiter.mass + europa.mass);
j.g = ganymede.mass / (jupiter.mass + ganymede.mass);
j.c = callisto.mass / (jupiter.mass + callisto.mass);

s.e = enceladus.mass / (saturn.mass + enceladus.mass);
s.m = mimas.mass / (saturn.mass + mimas.mass);
s.r = rhea.mass / (saturn.mass + rhea.mass); 
s.t = titan.mass / (saturn.mass + titan.mass);

u.t = titania.mass / (uranus.mass + titania.mass);
u.o = oberon.mass / (uranus.mass + oberon.mass);

n.t = triton.mass / (neptune.mass + triton.mass);

p.c = charon.mass / (pluto.mass + charon.mass);

%%% ----------------%%% ----------------%%% ----------------
fprintf('Earth-Moon\t\t%1.3e\n\n',e.m)
fprintf('Mars-Phobos\t\t%1.3e\n',m.p)
fprintf('Mars-Deimos\t\t%1.3e\n\n',m.d)
fprintf('Jupiter-Io\t\t%1.3e\n',j.i)
fprintf('Jupiter-Europa\t\t%1.3e\n',j.e)
fprintf('Jupiter-Ganymede\t%1.3e\n',j.g)
fprintf('Jupiter-Callisto\t%1.3e\n\n',j.c)
fprintf('Saturn-Enceladus\t%1.3e\n',s.e)
fprintf('Saturn-Mimas\t\t%1.3e\n',s.m)
fprintf('Saturn-Rhea\t\t%1.3e\n',s.r)
fprintf('Saturn-Titan\t\t%1.3e\n\n',s.t)
fprintf('Uranus-Titania\t\t%1.3e\n',u.t)
fprintf('Uranus-Oberon\t\t%1.3e\n\n',u.o)
fprintf('Neptune-Triton\t\t%1.3e\n\n',n.t)
fprintf('Pluto-Charon\t\t%1.3e\n\n',p.c)



%%% ----------------%%% ----------------%%% ----------------
% figure; hold all
% plot(em,0,'o','color',colors.new.blue)
% plot([mp md],[0,0],'o','color',colors.new.red)
% plot([ji je jg jc],[0 0 0 0],'o','color',colors.new.orng)
% plot([se sm st],[0 0 0],'o','color',colors.new.grn)
% plot(nt,0,'o','color',colors.new.ltblue)
% plot(pc,0,'o','color',colors.new.mag)
% 
% 


% %%% ----------------%%% ----------------%%% ----------------
% % Lagrange Points
% Ls = EquilibriumPoints(j.e)





