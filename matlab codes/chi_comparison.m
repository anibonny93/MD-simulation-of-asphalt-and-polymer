%plot(timeframe,chiall)
%hold on
%plot(timeframe_2,chiall_2)

%xlabel('1000 femtoseconds')
%ylabel('Chi value')
%export_fig 'Chi_comparison' -png -r800 -a1cd


figure
plot(timeframe, chiall,'k','MarkerSize',2)
hold on
plot(timeframe, chiall_AP_PE30_pppm,'r','MarkerSize',2)
plot(timeframe,ones(length(timeframe),1)*mean_chi_AP_PE30, 'k')
plot(timeframe,ones(length(timeframe),1)*mean_chi_AP_PE30_pppm,'r')

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
        legtext={'AP-PE ewald','PE-AP pppm'}
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
        export_fig 'Chi_comparison_AP_PE_kspace' -png -r800 -a1