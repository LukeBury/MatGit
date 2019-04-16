function [ colors ] = get_colors(testColors)

%==========================================================================
%==========================================================================
%
% Returns a structure containing lots of fun colors to use for plotting.
% Can also be used to create a test plot to see what different colors look
% like together
%
%
% Author: Eric W. Ferguson, Luke Bury
% Date: 08/26/16
%
%
% INPUT:         Description                                          Units
%  testColor     (optional)...[nx3] matrix of colors you wish to see
%                plotted
%
% OUTPUT:       
%    
%  colors        structure of colors where each field is a           struct
%                pallete and field with a pallete are RGB colors                
%
%==========================================================================
%==========================================================================
colors = struct;
%==========================================================================
%% Standard Colors
%==========================================================================
colors.std.blue       = [50  100 200]./255;
colors.std.ltblue     = [125 216 255]./255;
colors.std.orange     = [255 127 14]./255;
colors.std.ltorange   = [255 187 120]./255;
colors.std.grn        = [44 160 44]./255;
colors.std.ltgrn      = [152 223 138]./255;
colors.std.drkgrn     = [15 87 20]./255;
colors.std.red        = [209   0   0]./255;
colors.std.ltred      = [255 152 150]./255;
colors.std.drkred     = [150 40 32]./255;
colors.std.purp       = [148 103 189]./255;
colors.std.ltpurp     = [197 176 213]./255;
colors.std.brown      = [140 86 75]./255;
colors.std.ltbrown    = [196 156 148]./255;
colors.std.pink       = [227 119 194]./255;
colors.std.ltpink     = [247 182 210]./255;
colors.std.grey       = [127 127 127]./255;
colors.std.ltgrey     = [199 199 199]./255;
colors.std.ylwgrn     = [188 189 34]./255;
colors.std.khaki      = [219 219 141]./255;
colors.std.turq       = [23 190 207]./255;
colors.std.ltturq     = [158 218 229]./255;
colors.std.mag        = [255  0 255]./255;
colors.std.maglt      = [229 194 237]./255;
colors.std.ylw        = [247 202   0]./255;
colors.std.black      = [0 0 0]./255;
colors.std.cyan       = [0 255 255]./255;
colors.std.white      = [255 255 255]./255;

%==========================================================================
%% Color Schemes
%==========================================================================
%%%-------------------------------------
%%% Risk color schemes
%%%-------------------------------------
colors.sch.r6 = [179 19 19;240 80 83;246 144 61;247 212 32;164 224 87;25 214 29]./255;
colors.sch.r9 = [107 11 11;179 19 19;240 80 83;246 127 61;247 180 32;247 212 32;217 205 32;164 224 87;25 214 29]./255;

%%%-------------------------------------
%%% Distinct color schemes
%%%-------------------------------------
%%% 3 colors
% redish, blueish, greenish
colors.sch.d3_1 = [239 71 111; 0 136 204; 4 173 74]./255; 

% redish, yellowish, blueish (from wild cherry background)
colors.sch.d3_2 = [255 81 104; 255 181 32; 107 165 211]./255; 

%%% 4 colors
% 4 distinct colors #1 (from my NSTRF '18)
colors.sch.d4_1 = [6 214 160; 38 84 124; 239 71 111; 255 196 61]./255; 
colors.sch.d4_2 = [0, 167, 225; 203, 133, 99; 243, 126, 93; 3, 25, 30]./255;
colors.sch.d4_3 = [213, 51, 105; 203, 173, 109; 138, 142, 145; 184, 212, 227]./255;

%%%-------------------------------------
%%% Similar color schemes
%%%-------------------------------------
%%% 3 colors
% Blues
colors.sch.s3_1 = [0 43 65;25 100 126;40 175 176]./255; 

%%% 5 colors
% Pinkish gradient with cyan finisher
colors.sch.s5_1 = [254, 184, 212; 229, 153, 215; 187, 163, 211; 146, 181, 216; 1, 248, 217]./255;

%%%-------------------------------------
%%% Color Gradients
%%%-------------------------------------
% (black-blueish greenish - yellow)
colors.g5_1 = [74 66 56; 77 83 89; 80 132 132; 121 201 158; 207 234 119]./255;
% (icy blues)
colors.g5_2 = [149 225 226; 65 211 224; 19 184 221; 18 104 132; 14 26 45]./255;


%==========================================================================
%% Color Testing
%==========================================================================
if nargin ~= 0
    figure
    Y = ones(2,size(testColors,1));
    h = area(Y);
    for kk = 1:size(testColors,1)
        h(kk).FaceColor = testColors(kk,:);
    end
end

end

