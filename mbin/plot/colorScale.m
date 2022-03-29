function [ colorMatrix ] = colorScale(colorStack_in, n_colors_out )
%%% Description
% For getting a matrix of length 'n' and width 3 that provides codes for
% color gradients evenly spaced between some input stack of colors. The
% spacing is done with linear interpolation between colors
%       
% ------------------------------------------------------------------------
%%% Inputs
% colorStack_in - [nx3] rows of desired rgb color codes
% n_colors_out - [scalar] desired number of color steps for output
% ------------------------------------------------------------------------
%%% Outputs
% colorMatrix - [nx3] matrix where each row is a color, scaled between the
%                     input colors provided
% ------------------------------------------------------------------------
% Created: 02/16/21
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Determine number of colors in the input stack
n_colors_in = size(colorStack_in, 1);

%%% Quick catch for cases where the output should just be the input
if n_colors_in == n_colors_out
    colorMatrix = colorStack_in;
    return
end

%%% Determine the number of sections (areas between color inputs)
n_sections = n_colors_in - 1;

%%% Determine base number of colors per section and the remainder
colorsPerSection_floor = floor(n_colors_out / n_sections);
remainderColors        = rem(n_colors_out, n_sections);

%%% Preallocate a cell specifying the number of colors per section so that
%%% n_colors_out is met
colorsPerSection       = cell(1, n_sections);

%%% Loop through each section number and assign the appropriate number of
%%% colors. This number starts at the floor calculated earlier, and '1' is
%%% added to each successive section unti n_colors_out is met. This way,
%%% all sections are within '1' of the same number
for kk = 1:n_sections
    colorsPerSection{kk} = colorsPerSection_floor;

    if kk <= remainderColors
        colorsPerSection{kk} = colorsPerSection{kk} + 1;
    end
end

%%% Just using this while writing the code as a test
if ~isequal(sum([colorsPerSection{:}]), n_colors_out)
    warning('Problem...')
end

%%% Preallocate the output matrix of colors
colorMatrix = NaN(n_colors_out, 3);

%%% Loop through each section to assign all the necessary colors
for section_i = 1:n_sections
    %%% Assign a starting index in the output matrix for the current
    %%% section
    if section_i == 1
        startingIndex = 1;
    elseif section_i > 1
        startingIndex = sum([colorsPerSection{1:(section_i-1)}]) + 1;
    end

    %%% Assigning intial and final colors for this section. Linear
    %%% interpolation between these two colors is used to populate the
    %%% current section
    c11 = colorStack_in(section_i,1);   c12 = colorStack_in(section_i,2);   c13 = colorStack_in(section_i,3);
    c21 = colorStack_in(section_i+1,1); c22 = colorStack_in(section_i+1,2); c23 = colorStack_in(section_i+1,3);
    
    %%% Assign each color in the current section with linear interpolation
    for colors_in_section_i = 1:colorsPerSection{section_i}
        %%% Determine the current index in the output matrix
        currentIndex = startingIndex + colors_in_section_i - 1;
        
        %%% Assign each of the three rgb color codes to the current row of
        %%% the output matrix based on linear interpolation between the two
        %%% colors of the current section
        colorMatrix(currentIndex,1) = c11 + (colors_in_section_i-1)*(c21-c11)/colorsPerSection{section_i};
        colorMatrix(currentIndex,2) = c12 + (colors_in_section_i-1)*(c22-c12)/colorsPerSection{section_i};
        colorMatrix(currentIndex,3) = c13 + (colors_in_section_i-1)*(c23-c13)/colorsPerSection{section_i};
    end
end

%%% Set the final color
colorMatrix(end,:) = colorStack_in(end,:);
    
end




