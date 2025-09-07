close all 
clear all

fid = fopen('/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2024/AAM1_thermal_conductivity/tmp.profile', 'r');
fgetl(fid);
fgetl(fid);
fgetl(fid);
%fgetl(fid);
counter=0;
while ~feof(fid)
    clear_line = fgetl(fid);
    if ischar(clear_line)
        info = str2double(strsplit(clear_line));
        
        
%         for i=1:20;
%         chunks = zeros(i, 1);
%         coord1 = zeros(i, 1);
%         ncount = zeros(i, 1);
%         v_temp = zeros(i, 1);
%         end
       counter=counter+1;
       timestep(counter) = info(1);
        for i = 1:20
            clear_line = fgetl(fid);
            data = str2double(strsplit(clear_line));
            chunks(i,counter) = data(2);
            coord1(i,counter) = data(3);
            ncount(i,counter) = data(4);
            v_temp(i,counter) = data(5);
        end
      
        disp(['Timestep: ', num2str(timestep(counter))]);
        disp(['Chunks: ', num2str(chunks(:,counter)')]);
        disp(['v_temp: ', num2str(v_temp(:,counter)')]);
      
    end
end

figure
hold on
for i=1:10:counter
    plot(v_temp(:,i));
end
% i=counter;
 plot(v_temp(1:12,i-4)./0.00199, '-o','LineWidth',2)
 xlabel('Chunks','FontSize',16);
 ylabel('Temperature','Fontsize',16);
 gca.FontSize = 16;

 max_value=v_temp(11,i-4)./0.00199;
 min_value=v_temp(1,i-4)./0.00199;
 x_max=11;
 x_min=1;
 line([0 12], [max_value max_value]) %(max_value,'r--','LineWidth',2);
 line([0 12], [min_value min_value]) 
 text(11,max_value,['(11, ',num2str(max_value),')'],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',14,'Color','black');
 text(1,min_value,['(1, ',num2str(min_value),')'],'VerticalAlignment','top','HorizontalAlignment','left','FontSize',14,'Color','black');
%text(x_max,max_value,['Max: ('num2ster(x_max)','num2str(max_value)')'],'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10,'Color','red');
 
%plot(v_temp(:,i),0.00199);

% [temp_value,chunk_value]=ginput(2);
% disp(['Point 1:(', num2str(temp_value(1)),',',num2str(chunk_value(1)),')']);
% disp(['Point 1:(', num2str(temp_value(2)),',',num2str(chunk_value(2)),')']);

fclose(fid);