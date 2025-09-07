filepath_1 = '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/unoxidized/saturates_aromatics_gofr.dat';
filepath_2 = '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/all oxidized/saturates_oxi_aromatics_gofr.dat';
%filepath_3= '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/unoxidized/saturates_aromatics_gofr.dat';
%filepath_4= '/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Asphalt_run/unoxidized/saturates_asphaltenes_gofr.dat';

data_1 = readtable(filepath_1);
data_2 = readtable(filepath_2);
%data_3= readtable(filepath_3);
%data_4= readtable(filepath_4);

variable_names = {'r','gofr','none'};
data_1.Properties.VariableNames = variable_names;
data_2.Properties.VariableNames = variable_names;
%data_3.Properties.VariableNames = variable_names;
%data_4.Properties.VariableNames = variable_names;

r1 = data_1.r;
rdf1 = data_1.gofr;
r2 = data_2.r;
rdf2 = data_2.gofr;


plot (r1, rdf1, 'LineWidth', 3);

hold on; 

plot (r2,rdf2, 'LineWidth', 3);

%hold off;

%r3 = data_3.r;
%rdf3 = data_3.gofr;
%r4 = data_4.r;
%rdf4 = data_4.gofr;

%plot (r3, rdf3, 'LineWidth', 3);

%hold on; 

%plot (r4,rdf4, 'LineWidth', 3);

hold off;

A=gca;

x = 0:0.1:45;


yValue = 1;
line([min(x), max(x)], [yValue, yValue], 'LineWidth', 2,  'Linestyle', ':', 'Color', 'k');

% A.YScale='log'
 axis([0 36 0.0 2.5])
        xlabel('r(Å)')
        ylabel('g(r)')
        legtext={'Saturates-Aromatics', 'Saturates-Oxidized Aromatics'};%, 'Saturates-Aromatics','Saturates-Asphaltenes'};
        leg=legend(legtext);
        h1=gca;
        set(gcf, 'units', 'inches', 'pos', [3 3 3.5 3])
        legend boxoff
      %  set(gca,'FontSize',16)
        set(gcf,'color','w');
      %  set(gca,'color','None');
        box on
        set(gcf, 'Position', [10, 10, 10, 10]);
        hLegend = legend('Location','east','Orientation','vertical');
        legend_title = 'Pair distribution of solubility classes';
        title(legend_title, 'FontSize', 18, 'FontWeight','normal');
       
        set(leg,'FontSize',16)
      %  set(h1,'TickLength',[.02 .01])
      %  set(h1,'XMinorTick','on')
        %export_fig 'gofr_comparison_APy_PE' -png -r300 -a1