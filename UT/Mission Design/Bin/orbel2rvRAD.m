function [r, v] = orbel2rvRAD(mu,a,e,i,o,w,TA)


rmag = a*(1-e^2)/(1+e*cos(TA));
vmag = sqrt(mu*(2/rmag-1/a));
FPA = atan2(e*sin(TA),1+e*cos(TA));
taw=TA+w;
coso=cos(o);
sino=sin(o);
cosi=cos(i);
sini=sin(i);
costaw=cos(taw);
sintaw=sin(taw);
    r = rmag*[(costaw*coso - sintaw*cosi*sino);
        (costaw*sino + sintaw*cosi*coso);
        (sintaw*sini)];
    taw=taw-FPA;
    costaw=cos(taw);
    sintaw=sin(taw);
    v = vmag*[(-sintaw*coso - costaw*cosi*sino);
        (-sintaw*sino + costaw*cosi*coso);
        (costaw*sini)];
end