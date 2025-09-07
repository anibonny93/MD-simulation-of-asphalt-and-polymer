clear all
close all
path1='/home/anic/polymatic_v1.1/mix_files/';
cd(path1);
lmp_files=dir('*.lmp');
%lmp_files(1:2)= [];
total_files=numel(lmp_files);

for i=1:20
    file_name=lmp_files(i).name;
    fid_in=fopen(file_name, 'r');
    fgetl(fid_in);
    rest_of_file=fread(fid_in,'*uint8');
    fid_out=fopen(file_name, 'w');
    fprintf(fid_out, 'Random packing \n');
    fwrite(fid_out, rest_of_file);
    fclose(fid_in);
    fclose(fid_out);
end 