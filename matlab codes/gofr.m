filepath_ini = '/home/anic/MATLAB/BZB_PE30_gofr_ini.dat';
%filepath_equi = '/home/anic/MATLAB/BZB_PE30_gofr_equi.dat';
%filepath_oxi_ini = '/home/anic/MATLAB/BZB_2_PE30_gofr_ini.dat';
%filepath_oxi_equi = '/home/anic/MATLAB/BZB_2_PE30_gofr_equi.dat';

data_ini = readtable(filepath_ini);
%data_equi = readtable(filepath_equi);
%data_oxi_ini = readtable(filepath_oxi_ini);
%data_oxi_equi = readtable(filepath_oxi_equi);

variable_names = {'r','gofr','none'};
data_ini.Properties.VariableNames = variable_names;
%data_equi.Properties.VariableNames = variable_names;
%data_oxi_ini.Properties.VariableNames = variable_names;
%data_oxi_equi.Properties.VariableNames = variable_names;

r1 = data_ini.r;
rdf1 = data_ini.gofr;
%r2 = data_equi.r;
%rdf2 = data_equi.gofr;


plot (r1, rdf1, 'LineWidth', 3);

hold on; 

%plot (r2,rdf2, 'LineWidth', 3);

%hold off;

%r3 = data_oxi_ini.r;
%rdf3 = data_oxi_ini.gofr;
%r4 = data_oxi_equi.r;
%rdf4 = data_oxi_equi.gofr;

%plot (r3, rdf3, 'LineWidth', 3);

%hold on; 

%plot (r4,rdf4, 'LineWidth', 3);

hold off;

A=gca;

x = 0:0.1:45;


yValue = 1;
line([min(x), max(x)], [yValue, yValue], 'LineWidth', 2,  'Linestyle', ':', 'Color', 'k');

% A.YScale='log'
 axis([0 45 0.6 1.05])
        xlabel('r(Å)')
        ylabel('g(r)')
     %   legtext={'Unoxidized minimization', 'Unoxidized equilibration', 'Oxidized minimization','Oxidized equilibration'};
        leg=legend(legtext);
        h1=gca;
        set(gcf, 'units', 'inches', 'pos', [3 3 3.5 3])
        legend boxoff
        set(gca,'FontSize',16)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        set(gcf, 'Position', [10, 10, 10, 10]);
        hLegend = legend('Location','east','Orientation','vertical');
       % legend_title = 'Oxidized Benzobisbenzothiophene-Oxidized Polyethylene g(r) plot';
        title(legend_title, 'FontSize', 18, 'FontWeight','normal');
       
        set(leg,'FontSize',16)
        set(h1,'TickLength',[.02 .01])
        set(h1,'XMinorTick','on')
        %export_fig 'gofr_comparison_APy_PE' -png -r300 -a1