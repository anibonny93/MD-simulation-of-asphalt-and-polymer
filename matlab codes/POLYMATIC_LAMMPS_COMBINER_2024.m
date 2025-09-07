mix_files = '/home/anic/polymatic_v1.1/mix_files.txt';
fid =fopen(mix_files,'r');
if fid == -1
    error('Cannot open file:%s', mix_files);
end
filelocations = {};
while ~feof(fid)
    filelocation = strtrim(fgetl(fid));
    if ~isempty(filelocation)
        filelocations{end+1} = filelocation;
    end
end


% %POLYMATIC_LAMMPS_COMBINER_2023(filelist);