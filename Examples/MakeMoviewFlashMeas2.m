% Make a movie of the Itokawa Asteroid with given Measurement data already:
% 
% 
% Using VideoWriter 
% 
% Ann Dietrich
% 08/01/2014 
% Update: 08/10/2015
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% LOAD Measurements:
fn_orbmeas = 'OrbitIto_da300r120_X0ta1k_MeasFls64_FV32VR_5m.mat';

cd  '/Users/annabelle/Documents/a_docs/CU/LIDAR/Code2/OrbitsnMeas/';
load(fn_orbmeas);
    % Get: Xint, tint, XBFint, Obs_grid, sensor 
%cd ../

% Movie NaMe:
fn_movie = 'AstMovie_Itoa1ki90_VR1d5m_5fps.avi';



% Ending point?
tf = 86400; 
nt = length(tint); %find(tint==tf); %


% Figure set up
set(gcf, 'Color','white')  
set(gca, 'nextplot','replacechildren');
set(gca, 'CLim', [.10, 1]);
%set(gca, 'nextplot','replacechildren', 'Visible','off');


% New file: 
vidObj = VideoWriter(fn_movie);
vidObj.FrameRate = 5;
open(vidObj);


% Loop through R's V's and plot and make frames: 
%F(1:size(1:2:nt,2)) = struct('cdata',[], 'colormap',[]);
%im =1; 
for i = 1:nt 

    % time: 
    tc = tint(i); 
    
    % range:
    range = Obs_grid(i).range_noisy;
    
    % plot range 
    pcolor(range)
    axis tight
    
    title(sprintf('Range, Time: %.1f min',tc/60)); 
    h1 = colorbar; 
    ylabel(h1, 'Range (km)'); 
    shading flat 
    
    % Write each frame to the file.
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
       
    
end


 close(vidObj); 
 
 