% Open the input file for reading
fileID = fopen('/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2023/PE30/Complete_run/PE30_3.log', 'r');

% Open a temporary file for writing the modified content
tempID = fopen('/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2023/PE30/Complete_run/PE30_updated_3.log', 'w');

    % Check if the line contains the word "Step"
    [startIndex, endIndex, ~, ~, tokens] = regexp(line, 'Step (\d+)');
    if ~isempty(tokens)
        number = str2double(tokens{1});
        
        % Define the mapping
        mapping = [10 5; 11 6; 12 7; 13 8; 14 9];
        
        % Check if the number matches any mapping
        index = find(mapping(:, 1) == number, 1);
        if ~isempty(index)
            newNumber = mapping(index, 2);
            line = strrep(line, tokens{1}, num2str(newNumber));
        end
    end
    
    % Write the modified line to the temporary file
    fprintf(tempID, '%s\n', line);

% Close the input and temporary files
fclose(fileID);
fclose(tempID);

% Rename the temporary file to the original file
movefile('/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2023/PE30/Complete_run/PE30_updated_3.log', '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2023/PE30/Complete_run/PE30_3.log');
