%filelist=ls;
%[m,n]=size(filelist);
%for i=1:m;
%    fid=fopen(filename(i,1:end));
%end
%path1='/home/anic/polymatic_v1.1/mix_files/';
%cd(path1);
%lmp_files=dir('*.lmp');
%lmp_files(1:2)= [];
%total_files=numel(lmp_files);

fid_in=fopen('/home/anic/polymatic_v1.1/APy_PP_3.lmp', 'r');
fid_out=fopen('/home/anic/polymatic_v1.1/APy_PP_3_new.lmp', 'w');

%file_name=lmp_files(i).name;
 
%fid_in=fopen(file_name, 'r');
%fid_out=fopen(file_name, 'w');
fgetl(fid_in);
rest_of_file=fread(fid_in,'*uint8');
fwrite(fid_out, rest_of_file);
fclose(fid_in);
fclose(fid_out);
 

%    ofile=fopen(file{i});
%    o_line=fgetl(ofile);

