function [Results] = readlammpslog(filename)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


fid=fopen(filename);
stepcount=0;


x=fgetl(fid);
line_num=1;

while ~all(x==-1)
    if mod(line_num,10000)==0
        disp(['Reading line number ' num2str(line_num) '\n'])
    end
    %already read new x
    if startsWith(x,'---------------- Step') || startsWith(x, '------------ Step') 
        stepcount=stepcount+1;
        xsplit=split(x);
        Results.step(stepcount)=str2num(xsplit{3});
        Results.time(stepcount)=str2num(xsplit{7});
    elseif stepcount>0
        xsplit=split(x);
        for i=1:size(xsplit,1)
            if strcmp(xsplit{i},'TotEng')
                Results.TotEng(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'KinEng')
                Results.KinEng(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'Temp')
                Results.Temp(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'PotEng')
                Results.PotEng(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'E_bond')
                Results.E_bond(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'E_angle')
                Results.E_angle(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'E_dihed')
                Results.E_dihed(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'E_impro')
                Results.E_impro(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'E_vdwl')
                Results.E_vdwl(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'E_coul')
                Results.E_coul(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'E_long')
                Results.E_long(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'Press')
                Results.Press(stepcount)=str2num(xsplit{i+2});
            elseif strcmp(xsplit{i},'Volume')
                Results.Volume(stepcount)=str2num(xsplit{i+2});
            end
        end
    end
    x=fgetl(fid);
    if isempty(x)
        x=' ';
    end
    line_num=line_num+1;
end




