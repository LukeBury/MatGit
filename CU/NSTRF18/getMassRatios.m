function [Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto] = getMassRatios()
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
Earth.moon = moon.mass / (earth.mass + moon.mass);

Mars.phobos = phobos.mass / (mars.mass + phobos.mass);
Mars.deimos = deimos.mass / (mars.mass + deimos.mass);

Jupiter.io = io.mass / (jupiter.mass + io.mass);
Jupiter.europa = europa.mass / (jupiter.mass + europa.mass);
Jupiter.ganymede = ganymede.mass / (jupiter.mass + ganymede.mass);
Jupiter.callisto = callisto.mass / (jupiter.mass + callisto.mass);

Saturn.enceladus = enceladus.mass / (saturn.mass + enceladus.mass);
Saturn.mimas = mimas.mass / (saturn.mass + mimas.mass);
Saturn.rhea = rhea.mass / (saturn.mass + rhea.mass); 
Saturn.titan = titan.mass / (saturn.mass + titan.mass);

Uranus.titania = titania.mass / (uranus.mass + titania.mass);
Uranus.oberon = oberon.mass / (uranus.mass + oberon.mass);

Neptune.triton = triton.mass / (neptune.mass + triton.mass);

Pluto.charon = charon.mass / (pluto.mass + charon.mass);

end




