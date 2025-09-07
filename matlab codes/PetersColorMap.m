function [color]=PetersColorMap(i,type)
if exist('type','var')==0
    colors=[
        0 0 0;
        1 0 0;
        0 0 1;
        0 1 0;
        1 0 1
        1 0.5 0
        0 1 1];
    if i>size(colors,1)
        i=mod(i,size(colors,1));
        if i==0
            i=size(colors,1);
        end
    end
elseif strcmp(type,'original')
        colors=[
        0 0 0;
        1 0 0;
        0 0 1;
        0 1 0;
        1 0 1
        1 0.5 0
        0 1 1];
    if i>size(colors,1)
        i=mod(i,size(colors,1));
        if i==0
            i=size(colors,1);
        end
    end
elseif strcmp(type,'second')
    colors=[0 0 0;
        228 26 28;
        55 126 184;
        77 175 74;
        152 78 163;
        255 127 0;
        166 86 40;
        247 129 191;
        153 153 153]/255;
        if i>size(colors,1)
        i=mod(i,size(colors,1));
        if i==0
            i=size(colors,1);
        end
    end
elseif strcmp(type,'orange')
    colors=[
        255 29 29;
        255 85 65;
        255 51 0;
        255 118 0;
        255 155 0;
        255 191 4]/255;
elseif strcmp(type,'Diverging10ColorBrewer')
    colors=cbrewer('div','BrBG',10);
else
    error('invalid type');
end

color=colors(i,:);

end
