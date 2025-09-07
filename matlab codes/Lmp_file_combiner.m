close all 
clear all
fid = fopen('/home/anic/polymatic_v1.1/');
binder='DOCHN_1_mixed.lmp 1 ';
polymer='PE30_DOCHN_1_mixed.lmp 1 ';
mix =['cd /home/anic/polymatic_v1.1/; ./pack.pl -i 512 ']
%fileID = fopen('AP_PE30_mix_file.txt','w')
counter=0;
Limit=1;

    for i=1:512;
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
    mix = [mix '-l 180 -o DOCHN_PE30_v1_1.lmp']
    system(mix)