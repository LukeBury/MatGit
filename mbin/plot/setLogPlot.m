%%% Set axis of current plot as log scale
function setLogPlot(logAxis) 
if isequal(logAxis,'x')
    set(gca, 'XScale', 'log')
elseif isequal(logAxis,'y')
    set(gca, 'YScale', 'log')
elseif isequal(logAxis,'z')
    set(gca, 'ZScale', 'log')
end
end

