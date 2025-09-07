function [com,positions]=COM(positions,dim,masses)
%center of mass calculation
if size(positions,1)==1
%     warning('only 1 position found, returning that position')
    com=positions;
else

if exist('dim','var')==0
%     warning('no dimensions specified, assuming infinite boundaries')
    dim=[inf inf inf];
    dimflag=1;
else
    dimflag=0;
end

%     positions=pos(Chain,:); old
    
    distances=positions-positions(1,:);
    for i=1:3
       positions(distances(:,i)>dim(i)/2,i)=positions(distances(:,i)>dim(i)/2,i)-dim(i);
       positions(distances(:,i)<-dim(i)/2,i)=positions(distances(:,i)<-dim(i)/2,i)+dim(i);
    end
%     distances=sort(sum(distances.^2,2)).^0.5;
    
    if exist('masses','var')==0
        com=mean(positions);
    else
%         masses=mass(Chain,:); old
        for i=1:size(positions,1)
            masspos(i,:)=positions(i,:)*masses(i);
        end
        com=sum(masspos)/sum(masses);
    end
    if dimflag==0
        for i=1:3
            if com(i)>dim(i)
                com(i)=com(i)-dim(i);
            elseif com(i)<0
                com(i)=com(i)+dim(i);
            end
        end
    end
end
end