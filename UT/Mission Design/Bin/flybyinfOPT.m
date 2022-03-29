function [dvd]=flybyinfOPT(rpmag,rsoi,vinfinhat,vinfinmag,vinfouthat,vinfoutmag,mu,optr)
% BEST CASE
% warning off
% Determines the powered gravity assist trajectory to a parking orbit of

%------INPUTS-----
% rpmag: periapsis radius
% rP(3): Position of arrival planet
% vinfinhat(3): vinf incoming dir vector
% vinfinmag: vinf incoming magnitude
% vinfouthat(3): vinf outgoing dir vector
% vinfoutmag: vinf outgoing magnitude
% mu: gravitational parameter of planet
% optr: optimization of rp
%------OUTPUTS-----
% dvd: Delta-V required for hyperbolic depature trajectory


phir=acos(vinfouthat'*vinfinhat); %dot prod
phimax=2*asin(1/(1+rpmag*min([vinfinmag vinfoutmag])^2/mu)); %maximum available turning angle
phisoi=2*asin(1/(1+rsoi*max([vinfinmag vinfoutmag])^2/mu)); %max available turning angle at rsoi

if (phimax<phir)||~optr
    dvd=sqrt(vinfinmag^2+vinfoutmag^2-2*vinfinmag*vinfoutmag*cos(phir-phimax));
else
    if phir>phisoi
    dvd=abs(vinfoutmag-vinfinmag);
    else
    dvd=sqrt(vinfinmag^2+vinfoutmag^2-2*vinfinmag*vinfoutmag*cos(phir-phisoi));
    end
end

end