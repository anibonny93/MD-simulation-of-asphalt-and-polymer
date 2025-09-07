PP=[0	0.71 0.35 0.58 0 0 0.7 0.88 0 0 1.5 0.7 0.6 0.7 0.71 0 0 0.49 0.59]


PB=[0	0.23 -0.06 0.12 0 0 0.175 0.3 0 0 0.77 0.1 0.14 0.15 0.18 0 0 0.09 0.16]


PE500=[0	0.184 0.082 0.1482 0 0 0.048 0.190 0 0 0.4307 0.0341 0.107 0.158 0.069 0 0 -0.0419 0.042]
PE500_err_high=[0 0.202  0.1098 0.1770 0 0 0.077 0.257 0 0 0.485 0.070 0.1302 0.174 0 0 -0.0146 0.069]
PE500_err_low=[0 0.1666 0.0557 0.1194 0 0 0.0189 0.1229 0 0 0.375 -0.00185 0.0844 0.142 0.0478 0 0 -0.0692 0.0148] 

PE500(1)=mean(PE500(2:4));
PE500(6)=mean(PE500(7:8));
PE500(10)=mean(PE500(11:15));
PE500(17)=mean(PE500(18:19));
PB(1)=mean(PB(2:4));
PB(6)=mean(PB(7:8));
PB(10)=mean(PB(11:15));
PB(17)=mean(PB(18:19));
PP(1)=mean(PP(2:4));
PP(6)=mean(PP(7:8));
PP(10)=mean(PP(11:15));
PP(17)=mean(PP(18:19));
%row1=400
%row2=500

types={'Asphaltenes'
'Asphaltene Phenol'
'Asphaltene Pyrrole'
'Asphaltene Thiophene'
'-'
'Aromatics'
'Dioctylcyclohexane Naphthalene'
'Perhydrophenanthrene Napththalene'
' -- '
'Resins'
'Benzobisbenzothiophene'
'Trimethylbenzeneoxane'
'Thioisorenieratane'
'Quinolinohopane'
'Pyridinohopane'
'  ---  '
'Saturates'
'Squalane'
'Hopane'
};

X=categorical(types);
X=reordercats(X,types);

%%%%%%%%%%%%%%%%%%%%%%%%%%% All
        X=categorical(types);
X=reordercats(X,types);
 A=[PE500;PP;PB]'

 figure
b=bar(X,A);


colors=[[0 1 0];[0 0 1];[1 0 0]]
% colors=[[0 0 1];[0.5 0.5 1]]
for k = 1:size(b,2)
    b(k).FaceColor = colors(k,:)
end


A=gca;
% A.YScale='log'
%  axis([0 600 0 20])

% xlabel('T (K)')
        ylabel('\chi')
        legtext={'AAA-1 Our Results','AAA-1 + 10%PS','Greenfield 2014','Experimental','Greenfield 2007'}
        legtext={'PE 500K','PP 500K','PB500K'}
        leg=legend(legtext)
        h1=gca;
        set(gcf, 'units', 'inches', 'pos', [3.1895 3.9368 17.505 10.916])

        legend boxoff
        set(gca,'FontSize',16)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        set(leg,'FontSize',14);
        set(h1,'TickLength',[.02 .1])
        set(h1,'XMinorTick','on')
        
        %          export_fig 'Rt_example' -png -r800 -a1
        
        ylim([-0.2 1.5])
