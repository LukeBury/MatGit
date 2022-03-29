function [bifurcation_strings] = plot_BrouckeStabilityDiagram(alphas_input, betas_input, plot_diagram)
%%% Description
% Plots the broucke stability diagram for vectors of Broucke stability
% parameters, alpha and beta, and prints the location and type of
% bifurcations
%       
% ------------------------------------------------------------------------
%%% Inputs
% alphas_input - [1xn or nx1] Vector of Broucke stability parameter 'alpha' 
% betas_input - [1xn or nx1]  Vector of Broucke stability parameter 'beta' 
% plot_diagram - [logical] Whether or not to plot the diagram
% ------------------------------------------------------------------------
%%% Outputs
% (Process) Figure & print of bifurcation information
% bifurcation_strings - [cell array] Strings of bifurcation names
% ------------------------------------------------------------------------
% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Create the alpha-beta relationship lines that we want to watch for
%%% crosses
% -------------------------------------------------

if plot_diagram
    %%% Alpha vector to be used as independent variable in plot
    alphas_plot = linspace(min(alphas_input), max(alphas_input), 1500);

    %%% Five beta vectors that are dependent on alphas_IV
    betas_plot_pDoubling = 2.*alphas_plot - 2;
    betas_plot_tangent = -2.*alphas_plot - 2;
    betas_plot_pTripling = alphas_plot + 1;
    betas_plot_secondaryHopf = (alphas_plot.^2)/4 + 2;
    betas_plot_pQuadrupling = ones(1, length(alphas_plot)).*2;

    %%% Create the figure and plot the alpha-beta relationships
    figure; hold all
    p_tangent = plot(alphas_plot, betas_plot_tangent, 'b', 'linewidth',1);
    p_pDoubling = plot(alphas_plot, betas_plot_pDoubling, 'k', 'linewidth',1);
    p_pTripling = plot(alphas_plot, betas_plot_pTripling, 'r', 'linewidth',1);
    p_pQuadrupling = plot(alphas_plot, betas_plot_pQuadrupling, 'm', 'linewidth',1);
    p_secondaryHopf = plot(alphas_plot, betas_plot_secondaryHopf, ':k', 'linewidth',1);


    %%% Plot the actual alpha-beta data
    p_input = plot3(alphas_input, betas_input, linspace(1,length(alphas_input),length(alphas_input)), 'o','markersize',8,'markeredgecolor','k','markerfacecolor','r');
    
    PlotBoi2('$\alpha$','$\beta$', 26, 'LaTex')
    ylim([min(betas_input), max(betas_input)])
    xlim([min(alphas_input), max(alphas_input)])
    
    legend([p_tangent, p_pDoubling, p_pTripling, p_pQuadrupling, p_secondaryHopf, p_input], 'Tangent', 'Px2', 'Px3', 'Px4', 'Secondary Hopf', 'Input')
    
end

%%% Five beta vectors that are dependent on input alphas (from POs) to help
%%% find cross points
betas_atInputAlphas_pDoubling = 2.*alphas_input - 2;
betas_atInputAlphas_tangent = -2.*alphas_input - 2;
betas_atInputAlphas_pTripling = alphas_input + 1;
betas_atInputAlphas_secondaryHopf = (alphas_input.^2)/4 + 2;
betas_atInputAlphas_pQuadrupling = ones(length(alphas_input)).*2;

if plot_diagram
    %%% Plot of deltas from each bifurcation line
    figure('position',[943 52 654 947]);
    subplot(5,1,1); hold all
    plot(linspace(1,length(alphas_input),length(alphas_input)), betas_input - betas_atInputAlphas_pDoubling, 'b.', 'markersize',6)
    plot([0 length(alphas_input)], [0, 0], 'r', 'linewidth', 0.5)
    PlotBoi2('','P2 $\Delta \beta$',26,'LaTex')
    view(0,90)
    
    subplot(5,1,2); hold all
    plot(linspace(1,length(alphas_input),length(alphas_input)), betas_input - betas_atInputAlphas_pTripling, 'b.', 'markersize',6)
    plot([0 length(alphas_input)], [0, 0], 'r', 'linewidth', 0.5)
    PlotBoi2('','P3 $\Delta \beta$',26,'LaTex')
    view(0,90)
    
    subplot(5,1,3); hold all
    plot(linspace(1,length(alphas_input),length(alphas_input)), betas_input - betas_atInputAlphas_pQuadrupling, 'b.', 'markersize',6)
    plot([0 length(alphas_input)], [0, 0], 'r', 'linewidth', 0.5)
    PlotBoi2('','P4 $\Delta \beta$',26,'LaTex')
    view(0,90)
    
    subplot(5,1,4); hold all
    plot(linspace(1,length(alphas_input),length(alphas_input)), betas_input - betas_atInputAlphas_tangent, 'b.', 'markersize',6)
    plot([0 length(alphas_input)], [0, 0], 'r', 'linewidth', 0.5)
    PlotBoi2('','T $\Delta \beta$',26,'LaTex')
    view(0,90)
    
    subplot(5,1,5); hold all
    plot(linspace(1,length(alphas_input),length(alphas_input)), betas_input - betas_atInputAlphas_secondaryHopf, 'b.', 'markersize',6)
    plot([0 length(alphas_input)], [0, 0], 'r', 'linewidth', 0.5)
    PlotBoi2('Index','SH $\Delta \beta$',26,'LaTex')
    view(0,90)
    
end

%%% Initialize counts for various types of bifurcations
count_tangent           = 0;
count_periodDoubling    = 0;
count_periodTripling    = 0;
count_periodQuadrupling = 0;
count_secondaryHopf     = 0;

%%% In this loop, all we're doing is tracking whether our alpha-beta
%%% relationship crosses any of the other relationships (1-5). Whenever a
%%% cross happens, a bifurcation occurs. This events are tracked and the
%%% specific type of bifurcation is printed
bifurcation_strings = cell(0);
for kk = 1:length(betas_input)
    %%% Calculate current sign
    if kk == 1
        previousSign_pDoubling = sign(betas_input(kk) - betas_atInputAlphas_pDoubling(kk));
        previousSign_tangent = sign(betas_input(kk) - betas_atInputAlphas_tangent(kk));
        previousSign_pTripling = sign(betas_input(kk) - betas_atInputAlphas_pTripling(kk));
        previousSign_secondaryHopf = sign(betas_input(kk) - betas_atInputAlphas_secondaryHopf(kk));
        previousSign_pQuadrupling = sign(betas_input(kk) - betas_atInputAlphas_pQuadrupling(kk));
        continue
    else
        currentSign_pDoubling = sign(betas_input(kk) - betas_atInputAlphas_pDoubling(kk));
        currentSign_tangent = sign(betas_input(kk) - betas_atInputAlphas_tangent(kk));
        currentSign_pTripling = sign(betas_input(kk) - betas_atInputAlphas_pTripling(kk));
        currentSign_secondaryHopf = sign(betas_input(kk) - betas_atInputAlphas_secondaryHopf(kk));
        currentSign_pQuadrupling = sign(betas_input(kk) - betas_atInputAlphas_pQuadrupling(kk));
    end
    
    %%% Checks for whether signs have changed
    if (currentSign_pDoubling ~= previousSign_pDoubling) && (currentSign_pDoubling ~= 0)
        prevDiff    = betas_input(kk-1) - betas_atInputAlphas_pDoubling(kk-1);
        currentDiff = betas_input(kk) - betas_atInputAlphas_pDoubling(kk);
        
        if abs(prevDiff) <= abs(currentDiff)
            closerIndex = kk-1;
        else
            closerIndex = kk;
        end
        count_periodDoubling = count_periodDoubling + 1;
        bifurcation_strings = [bifurcation_strings; sprintf('%1dP2', count_periodDoubling)];
        fprintf('Bifurcation: %18s (%1dP2) ... PO index %1d-%1d ... Closer to %1d\n', 'Period Doubling', count_periodDoubling, kk-1, kk, closerIndex)        
    end

    if (currentSign_tangent ~= previousSign_tangent) && (currentSign_tangent ~= 0)
        prevDiff    = betas_input(kk-1) - betas_atInputAlphas_tangent(kk-1);
        currentDiff = betas_input(kk) - betas_atInputAlphas_tangent(kk);

        if abs(prevDiff) <= abs(currentDiff)
            closerIndex = kk-1;
        else
            closerIndex = kk;
        end
        count_tangent = count_tangent + 1;
        bifurcation_strings = [bifurcation_strings; sprintf('%1dT', count_tangent)];
        fprintf('Bifurcation: %18s (%1dT) ... PO index %1d-%1d ... Closer to %1d\n', 'Tangent', count_tangent, kk-1, kk, closerIndex)
    end

    if (currentSign_pTripling ~= previousSign_pTripling) && (currentSign_pTripling ~= 0)
        prevDiff    = betas_input(kk-1) - betas_atInputAlphas_pTripling(kk-1);
        currentDiff = betas_input(kk) - betas_atInputAlphas_pTripling(kk);

        if abs(prevDiff) <= abs(currentDiff)
            closerIndex = kk-1;
        else
            closerIndex = kk;
        end
        count_periodTripling = count_periodTripling + 1;
        bifurcation_strings = [bifurcation_strings; sprintf('%1dP3', count_periodTripling)];
        fprintf('Bifurcation: %18s (%1dP3) ... PO index %1d-%1d ... Closer to %1d\n', 'Period Tripling', count_periodTripling, kk-1, kk, closerIndex)
    end

    if (currentSign_secondaryHopf ~= previousSign_secondaryHopf) && (currentSign_secondaryHopf ~= 0)
        prevDiff    = betas_input(kk-1) - betas_atInputAlphas_secondaryHopf(kk-1);
        currentDiff = betas_input(kk) - betas_atInputAlphas_secondaryHopf(kk);

        if abs(prevDiff) <= abs(currentDiff)
            closerIndex = kk-1;
        else
            closerIndex = kk;
        end
        count_secondaryHopf = count_secondaryHopf + 1;
        bifurcation_strings = [bifurcation_strings; sprintf('%1dSH', count_secondaryHopf)];
        fprintf('Bifurcation: %18s (%1dSH) ... PO index %1d-%1d ... Closer to %1d\n', 'Secondary Hopf', count_secondaryHopf, kk-1, kk, closerIndex)
    end

    if (currentSign_pQuadrupling ~= previousSign_pQuadrupling) && (currentSign_pQuadrupling ~= 0)
        prevDiff    = betas_input(kk-1) - betas_atInputAlphas_pQuadrupling(kk-1);
        currentDiff = betas_input(kk) - betas_atInputAlphas_pQuadrupling(kk);

        if abs(prevDiff) <= abs(currentDiff)
            closerIndex = kk-1;
        else
            closerIndex = kk;
        end
        count_periodQuadrupling = count_periodQuadrupling + 1;
        bifurcation_strings = [bifurcation_strings; sprintf('%1dP4', count_periodQuadrupling)];
        fprintf('Bifurcation: %18s (%1dP4) ... PO index %1d-%1d ... Closer to %1d\n', 'Period Quadrupling', count_periodQuadrupling, kk-1, kk, closerIndex)
    end
    
    
    %%% Assigning current signs to previous signs
    previousSign_pDoubling = currentSign_pDoubling;
    previousSign_tangent = currentSign_tangent;
    previousSign_pTripling = currentSign_pTripling;
    previousSign_secondaryHopf = currentSign_secondaryHopf;
    previousSign_pQuadrupling = currentSign_pQuadrupling;
end
if (count_tangent + count_periodDoubling + count_periodTripling + count_periodQuadrupling + count_secondaryHopf) > 0
    fprintf('To print bifurcation list: fprintf(''%%s\\n'',bifurcation_strings{:})\n')
end
end % function