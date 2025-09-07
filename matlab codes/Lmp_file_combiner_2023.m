close all 
%clear all
bfile=fopen('/home/anic/MATLAB/lmp_file_combiner_list/lmp_binder_small.txt');
pfile=fopen('/home/anic/MATLAB/lmp_file_combiner_list/lmp_polymer_small.txt');
ofile=fopen('/home/anic/MATLAB/lmp_file_combiner_list/output_filename_small.txt');
b_line=fgetl(bfile);
p_line=fgetl(pfile);
o_line=fgetl(ofile);
iteration=1:15;

for in=iteration
    binder=b_line;
    polymer=p_line;
    output=o_line;
    mix =['cd /home/anic/polymatic_v1.1/; ./pack.pl -i 512 ']
%fileID = fopen('AP_PE30_mix_file.txt','w')/
    counter=0;
    Limit=1;
    for i=1:512
        counter=counter+1;
           if counter<=Limit
               mix = [mix binder]
               %fprintf(fileID, binder);
           else
               mix = [mix polymer]
            %fprintf(fileID, polymer);
               counter=0;  
           end
    end
        mix = [mix '-l 180 -o ' output]
        system(mix)
b_line=fgetl(bfile);
p_line=fgetl(pfile);
o_line=fgetl(ofile);
end



