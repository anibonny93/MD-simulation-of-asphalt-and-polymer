%plot(timeframe,chiall)
%hold on
%plot(timeframe_2,chiall_2)

%xlabel('1000 femtoseconds')
%ylabel('Chi value')
%export_fig 'Chi_comparison' -png -r800 -a1cd

%line_1= '#0072BD';
%line_2= '#D95319';
%line_3= '#EDB120';
%line_4= '#77AC30';
%color_1= sscanf(line_1(2:end),'%2x%2x%2x',[1 3])/255;

figure
plot(timeframe, chiall_APy_PE30_pppm,'r','MarkerSize',2)
hold on
plot(timeframe, chiall_APy_PE30_pppm_alt,'b','MarkerSize',2)
plot(timeframe, chiall_PE30_APy_pppm,'k','MarkerSize',2)
%plot(timeframe, chiall_PE30_AP_pppm_alt,'k','MarkerSize',2)

plot(timeframe,ones(length(timeframe),1)*mean(chiall_APy_PE30_pppm), 'r')
plot(timeframe,ones(length(timeframe),1)*mean(chiall_APy_PE30_pppm_alt),'b')
plot(timeframe,ones(length(timeframe),1)*mean(chiall_PE30_APy_pppm), 'k')
%plot(timeframe,ones(length(timeframe),1)*mean(chiall_PE30_AP_pppm_alt), 'k')

%chi_1=mean(chiall_PVC5_AP_pppm)
%chi_2=mean(chiall_PVC10_AP_pppm)
%chi_3=mean(chiall_PVC15_AP_pppm)
%chi_4=mean(chiall_PVC20_AP_pppm)
%chi_5=mean(chiall_PVC25_AP_pppm)
%plot(Tetramer_syndiotactic(:,2),Tetramer_syndiotactic(:,9),'r*','MarkerSize',2)
A=gca;
% A.YScale='log'
%  axis([0 5e6 -1e6 5e6 ])
xlabel('1000 femtoseconds')
        ylabel('Interaction parameter (\chi)')
        %legtext={'Discrete','Diblock','Triblock (Out)','Triblock (in)','Octamer','Dimer','Tetramer Alt','Triblock 2-8-2 Alt','Triblock 8-2-8 Alt'}
%         legtext={'Tetramer','Tetramer Alternating'}
        %legtext={'Tetramer','8-2-8','2-8-2','Tetramer Alternating','8-2-8 Alt','2-8-2 Alt'}
        %legtext={'Tetramer','Octamer','Dimer'}
        legtext={'APy packed first', 'Alternate packing' ,'PE packed first'} %'PE30-SQ'} %, 'AP-PE30 (0.75)', 'AP-PE30 (1.00)'}
        leg=legend(legtext)
        h1=gca;
        set(gcf, 'units', 'inches', 'pos', [25.1789 8 8.0947 6.7368])
        legend boxoff
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        set(leg,'FontSize',12);
        set(h1,'TickLength',[.02 .1])
        set(h1,'XMinorTick','on')
        export_fig 'Chi_comparison_APy_PE' -png -r800 -a1