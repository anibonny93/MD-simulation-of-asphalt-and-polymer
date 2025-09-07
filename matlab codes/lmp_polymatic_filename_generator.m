
%This will create the asphalt binder list

prefixes = {'APy','AT','BZB','DOCHN','PHPN','PH','QH','TIR','TBO'};
d_values = {'d1','d5','d10','d30'};
fileID = fopen('lmp_binder_2024.txt','w');
suffix ='mixed.lmp';
rep = 1;
for i = 1:length(prefixes)
    for j =1:length(d_values)
        for k=1:5
        filename =sprintf('%s_%s_%s %d ', prefixes{i}, d_values{j},suffix, rep);
        fprintf(fileID,'%s\n',filename);
        end
    end
end
fclose(fileID);
disp('File names populated successfully in lmp_binder_2024.txt');
clc
%This will create the polymer list

polymer = 'PE30';
fileID = fopen('lmp_polymer_2024.txt','w');
for i = 1:length(prefixes)
    for j =1:length(d_values)
        for k=1:5
        filename = sprintf('%s_%s_%s_%s %d ', polymer, prefixes{i}, d_values{j},suffix,rep);
        fprintf(fileID,'%s\n',filename);
        end
    end
end
fclose(fileID);
disp('File names populated successfully in lmp_polymer_2024.txt');


%This will create the output list

fileID = fopen('output_filenames_2024.txt','w');
%number = {'1','2','3','4','5'};

for i = 1:length(prefixes)
    for j =1:length(d_values)
        number = 1;
        for k=1:5
            filename = sprintf('%s_%s_%s_%d_%s', prefixes{i}, polymer, d_values{j}, number, suffix);
            fprintf(fileID,'%s\n',filename);
            number = number + 1;
        end
    end
end

fclose(fileID);


disp('File names populated successfully in output_filenames_2024.txt');




        