PP=[0	0.71 0.35 0.58 0 0 0.7 0.88 0 0 1.5 0.7 0.6 0.7 0.71 0 0 0.49 0.59]


PB=[0	0.23 -0.06 0.12 0 0 0.175 0.3 0 0 0.77 0.1 0.14 0.15 0.18 0 0 0.09 0.16]


PE500=[0	.132	-.026 0.038 0 0 0 0.18 0 0 0.405 0.03 0.07 0.055 0.088 0 0 0.01 0.03]


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
hold on
[ngroups,nbars]=size(A)

% x=nan(nbars,ngroups);
% for i=1:nbars
%     x(i,:)=b(i).EndPoints;
% end
groupwidth=min(0.8,nbars/(nbars+1.5));
Errors=[SE_low',SE_low',SE_low'];
for i=1:nbars
    x=(1:ngroups)-groupwidth/2+(2*i-1)*groupwidth/(2*nbars);
    errorbar(x+i,A(:,i),Errors(:,i),Errors(:,i),'k','linestyle','none');
end

% errorbar(x',A,0.1*ones(nabrs,ngroups),'k','linestyle','none')


% er1 = errorbar(X,PE500,SE_low,SE_high);
% er2= errorbar(X,PB,SE_low,SE_high);
% er3= errorbar(X,PP,SE_low,SE_high);


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
