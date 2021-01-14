%% _____________________________Parameters______________________________ %%
simulation_start_time = 0;
simulation_end_time = 3000;
simulation_time_step = 100;
% Cutoff radius in Agstrom
cutoff_bond_radius = 2;
% Network former labels
% This can also be a string input as shown next to the variable depending
% on your data
% Always ensure that 'O' or the code for oxygen (1) is the first element
network_former_labels = [1,2,3]; % [O, Si, Al]
% Modifier labels
% This can also be a string input as shown next to the variable depending
% on your data
modifier_labels = [4,5]; % [Ca, Na]
% Depth of surface to which constraint counting takes place in angstroms
surface_depth = 6;
% Resolution of the grid required (ngrid x ngrid)
ngrid = 10;

%% __________________________Count Constraints__________________________ %%

frames = 0; % Index the number of timesteps 
for n = simulation_start_time:simulation_time_step:simulation_end_time
    % Time series M.D data for a mixed modifier aluminosilicate glass
    % exposed to water on the top and bottom surfaces.
    % NOTE: This is not a necessary step as long as you have the following
    % data: 
    % #labels (nx1 double) - column vector with 
    % #coordinates (nx3 double) - columnvectors of x,y,z coordinates
    % #scale - (3x1 double) - Refer documentation for 
    % SurfMat.generate_coordinates by typing 
    % "help SurfMat.generate_coordinates" in the command window
    
    addpath('Aluminosilicate_example_data')
    filename = ['surface.', num2str(n), '.xyz'];
    
    % Import and Prepare Data
    [labels, coordinates, scale] = SurfMat.generate_coordinates(filename);
    
    % Delete Modifiers and NBOs %
    [labels, coordinates] = SurfMat.delete_modifiers(labels, coordinates, modifier_labels);
    [labels, coordinates, coordination] = SurfMat.get_coordination(labels,coordinates,cutoff_bond_radius,scale, network_former_labels);
    
    % Constraint counting %
    % Generate Constraint Density Map Data %
    [top_surface(:,:,((n/simulation_time_step)+1)), bottom_surface(:,:,((n/simulation_time_step)+1)), cs_x, cs_y] = ...
        SurfMat.count_surface_constraints(coordinates, coordination, scale);
    frames = frames + 1;
end

%% ____________________________Make a Movie!____________________________ %%

angs = string(char(197));
angsx = [string('x (') + angs + string(')')];
angsy = [string('y (') + angs + string(')')];
time = simulation_start_time:simulation_time_step:simulation_end_time;

%Contour Plot Movie Top surface
for i = 1:frames
    figure(1);
    contourf(cs_x,cs_y,top_surface(:,:,i));
    hold on ;
    tm = ['Time = ', num2str(time(i)), 'fs'];
    title(tm);
    xlabel(angsx);
    ylabel(angsy);
    set(gca, 'FontSize', 13, 'FontName', 'Times New Roman');
    title(tm,'FontSize', 16);
    colorbar;
    caxis([0 7]);
    axis square;
    F(i) = getframe(gcf);
    drawnow
end

% Save your movie! %
% create the video writer with 0.5 fps
writerObj = VideoWriter('glass1.avi');

% set the framerate
writerObj.FrameRate = 5;

% open the video writer
open(writerObj);

% write the frames to the video
clear i
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;
    writeVideo(writerObj, frame);
end

%close the writer object
close(writerObj);