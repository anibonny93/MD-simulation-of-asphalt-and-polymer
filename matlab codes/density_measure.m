mass = 133813.156
filenames ={'/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/all oxidized/data.log', 
    '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/half_oxidized/data.log',
    '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/unoxidized/test.log',
    '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/Aromatics_oxidized/data.log',
    '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/Resins_oxidized/data.log'};
    
figure;
hold on;

for i =1:length(filenames)
    fid = fopen(filenames{i}, 'r');
    if fid == -1
        error('Cannot open file:%s', filename{i});
    end
    header_line =fgetl(fid);
    headersplit =strsplit(strtrim(header_line));
    volume_column = find(strcmp(headersplit,'Volume'));
    volumes = [];
    while ~feof(fid)
        line =fgetl(fid);
        if isempty(line) || contains (line,'Volume')
            continue;
        end
        tokens = strsplit(strtrim(line));
        if length(tokens) >= volume_column
            volume = str2double (tokens{volume_column});
        if ~isnan(volume)
            volumes = [volumes; volume];
        end
        end
    end
fclose(fid);
fprintf('Number of volume values read: %d\n',length(volumes));
if isempty(volumes)
    error('No volume values read from the file');
end


densities= ((mass*1.661)./volumes);
density_avg(i)= mean(densities(1002:2002));
time =1:length(densities);
plot(time,densities);
end

xlabel('Time step');
ylabel('Density (g/cm^3)');
title ('Density over time(ps)');
output_filename='density_over_time.txt';
legend('All oxidized','Half oxidized','Unoxidized','Aromatics oxidized','Resins oxidized');
hold off;
export_fig 'Density_over_time' -png -r300