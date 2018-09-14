function [r,v] = pleph2o(moon,jd)
muj = 1.26712764e+008;      %km3/s2  
dxtol=1E-7;

switch moon
    case 1      % IO (501)
        ac =  4.223264326785622E+05;         % [km]
        ec =  3.226427658892863E-03;
        ic =  2.245152452327303E+00;       % [deg]
        oc =  3.376623143769770E+02;       % [deg]
        wc =  7.188513760538278E+01;       % [deg]
        tac = 1.486567216517267E+02;        % [deg]
        jd0 = 2461771.5; % [days] 2028-01-01

    case 2      % Europa (502)
        ac =  6.714298550202193E+05;         % [km]
        ec =  9.124612688122750E-03;
        ic =  1.964512487994255E+00;       % [deg]
        oc =  3.274313364072149E+02;       % [deg]
        wc =  2.550843886415831E+02;       % [deg]
        tac = 2.577399892528971E+02;        % [deg]
        jd0 = 2461771.5; % [days] 2028-01-01

    case 3      % Ganymede (503)
        ac =  1.070181784520349E+06;         % [km]
        ec =  2.516490640664923E-03;
        ic =  2.348102576090548E+00;       % [deg]
        oc =  3.388041093837463E+02;       % [deg]
        wc =  1.648371670664958E+01;       % [deg]
        tac = 3.577388871943164E+02;        % [deg]
        jd0 = 2461771.5; % [days] 2028-01-01
    
    case 4      % Callisto (504)
        ac =  1.882631993643279E+06;         % [km]
        ec =  7.194661573100830E-03;
        ic =  1.949405019340739E+00;       % [deg]
        oc =  3.368681508016559E+02;       % [deg]
        wc =  3.423601275905439E+01;       % [deg]
        tac = 3.444350979779773E+02;        % [deg]
        jd0 = 2461771.5; % [days] 2028-01-01
    
    otherwise
        error('This moon is not defined for the pleph2o function');
end

    

% --------------------Keplerian elements------------------


    [rP,vP]=orbel2rv(muj,ac,ec,ic,oc,wc,tac);
    [r,v]=kepuv(muj,(jd-jd0)*86400,rP,vP,0,dxtol);
end