%plot(timeframe,chiall)
%hold on
%plot(timeframe_2,chiall_2)

%xlabel('1000 femtoseconds')
%ylabel('Chi value')
%export_fig 'Chi_comparison' -png -r800 -a1cd


figure
plot(timeframe, chiall,'k','MarkerSize',2)
hold on
plot(timeframe, chiall2,'r','MarkerSize',2)
%plot(timeframe, chiall_2,'b','MarkerSize',2)
%plot(timeframe, chiall_3,'g','MarkerSize',2)
%plot(timeframe, chiall_4,'m','MarkerSize',2)
plot(timeframe,ones(length(timeframe),1)*mean(chiall), 'k')
plot(timeframe,ones(length(timeframe),1)*mean(chiall2),'r')
%plot(timeframe,ones(length(timeframe),1)*mean(chiall_2), 'b')
%plot(timeframe,ones(length(timeframe),1)*mean(chiall_3), 'g')
%plot(timeframe,ones(length(timeframe),1)*mean(chiall_4), 'm')
chi_1=mean(chiall)
chi_2=mean(chiall2)
%chi_3=mean(chiall_3)
%chi_4=mean(chiall_4)
%chi_5=mean(chiall_5)
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
        legtext={'AP & PVC10', 'AP & PVC20'}
        %legtext={'Timestep size 0.10, chi value aditw', 'Timestep size 0.20, chi value '.format(chi_1), 'Timestep size 0.60', 'Timestep size 0.80', 'Timestep size 1' }
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
        export_fig 'Chi_comparison' -png -r800 -a1