close all 
clear all

% Define the filename
filename = '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/SQ_new/tmp.profile';  % Replace with your actual filename

% Open the file
fid = fopen(filename, 'r');
fgetl(fid);
fgetl(fid);

% Read the first line (timesteps, number-of-chunks, total-count)
firstLine = fgetl(fid);
info = str2num(firstLine);
timestep = info(1);
data = textscan(fid, '%f %f %f', 'HeaderLines', 0);

% Close the file
fclose(fid);

% Extract the data
chunk = data{1};
coord_1 = data{2};
ncount = data{3};
v_temp = data{4};

% Extract timesteps and number of chunks from the first line
info = str2num(firstLine);  % Convert the line to numeric array

if numel(info) >=3
    timestep =info(1);
else
    warning('Data line does not have enough elements.');
end
%timestep = info(1);
num_chunks = info(2);

% Determine the number of timesteps by counting the number of chunks in the file
num_timesteps = length(v_temp) / num_chunks;

% Reshape the data into a matrix where rows are timesteps and columns are chunks
v_temp_matrix = reshape(v_temp, [num_chunks, num_timesteps])';

% Generate time points for the x-axis (replace with your actual time points if available)
time_points = 1:num_timesteps;

% Plot the data
figure;
hold on;
for i = 1:num_chunks
    plot(time_points, v_temp_matrix(i, :), 'DisplayName', sprintf('Chunk %d', i));
end
hold off;

xlabel('Timestep');
ylabel('Temperature (v_temp)');
title('Temperature of Chunks Over Time');
legend('show');
grid on;

% Save the figure
saveas(gcf, 'temperature_plot.png');  % Save the plot as a PNG file