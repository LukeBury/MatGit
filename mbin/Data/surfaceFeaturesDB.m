function [surfaceFeatures] = surfaceFeaturesDB()

% Data from:
% https://planetarynames.wr.usgs.gov/nomenclature/SearchResults;jsessionid=26E3FC2425B5EE1CA341C1047FFF2C19

% ------------------------------------------------------------------------
%%% Europa
% ------------------------------------------------------------------------
surfaceFeatures.europa.name = 'Europa';
% features = {Name, [lat, lon](deg 0-360), diameter (km)}
surfaceFeatures.europa.features = ...
    {'Thera Macula',           [-46.7, 181.2],     95;...
     'Thrace Macula',          [-45.9, 172.1],     180.2;...
     'Puddle',                 [],                 NaN;...
     'Callanish',              [-16.7, 334.5],     107;...
     'Rhadamanthys Linea',     [19.3, 200.5],      1747;...
     'Galileo Plume',          [-5, 115],          NaN};

end




