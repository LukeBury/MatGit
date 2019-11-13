function [primary, secondary] = assignPrimaryAndSecondary_CR3BP(familyTag, bodies)
%%% Description
%       Assigns structs of body data to a primary and secondary body meant
%       for use in the CR3BP. 
%       
% ------------------------------------------------------------------------
%%% Inputs
%       familyTag - [str] A string containing the names of the desired
%                         primary and secondary bodies according to the
%                         convention 'Primary_Secondary'. For example,
%                         familyTag might equal
%                         "PO_ICs.Jupiter_Europa.CR3BP.L2_Vertical"
%       bodies    - [struct] A struct containing all the loaded bodies from
%                            body-data file. Ex:
%                            'bodies = getBodyData(mbinPath)'
% ------------------------------------------------------------------------
%%% Outputs
%       primary   - [struct] Struct of data for primary body
%       secondary - [struct] Struct of data for secondary body
% ------------------------------------------------------------------------
% Created: 09/26/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Assign Primary
% -------------------------------------------------
if contains(lower(familyTag),'earth_')
    primary = bodies.earth;
elseif contains(lower(familyTag),'mars_')
    primary = bodies.mars;
elseif contains(lower(familyTag),'jupiter_')
    primary = bodies.jupiter;
elseif contains(lower(familyTag),'saturn_')
    primary = bodies.saturn;
elseif contains(lower(familyTag),'uranus_')
    primary = bodies.uranus;
elseif contains(lower(familyTag),'neptune_')
    primary = bodies.neptune;
else
    warning('No primary')
    return
end

% -------------------------------------------------
%%% Assign Secondary
% -------------------------------------------------
if contains(lower(familyTag),'_moon')
    secondary = bodies.moon;
elseif contains(lower(familyTag),'_phobos')
    secondary = bodies.phobos;   
elseif contains(lower(familyTag),'_io')
    secondary = bodies.io;    
elseif contains(lower(familyTag),'_europa')
    secondary = bodies.europa;
elseif contains(lower(familyTag),'_ganymede')
    secondary = bodies.ganymede;
elseif contains(lower(familyTag),'_callisto')
    secondary = bodies.callisto;    
elseif contains(lower(familyTag),'_enceladus')
    secondary = bodies.enceladus;
elseif contains(lower(familyTag),'_titan')
    secondary = bodies.titan;    
elseif contains(lower(familyTag),'_triton')
    secondary = bodies.triton;    
elseif contains(lower(familyTag),'_cordelia')
    secondary = bodies.cordelia;
elseif contains(lower(familyTag),'_ophelia')
    secondary = bodies.ophelia;
else
    warning('No Secondary')
    return
end

end % function