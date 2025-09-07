PVC_AP=[-0.0596934882	0.0808422828	-0.0668618507	-0.0525251257
0.8166527459	0.0694112845	0.8104979809	0.8228075108
-0.1121568965	0.0821165273	-0.1194382475	-0.1048755455
-0.0824270015	0.0779805728	-0.0893416134	-0.0755123896
-0.0922155989	0.0634345187	-0.0978403983	-0.0865907995]


%row1=400
%row2=500

types={'PVC5 with AP'
'PVC10 with AP'
'PVC15 with AP'
'PVC20 with AP'
'PVC25 with AP'
};

X=categorical(types);
X=reordercats(X,types);

%%%%%%%%%%%%%%%%%%%%%%%%%%% All
        X=categorical(types);
X=reordercats(X,types); 
A=[PVC_AP(:,1)]

 figure
 hold on
b=bar(X,A)
colors=[[0 1 0];[0 0 1];[1 0 0]]
% colors=[[0 0 1];[0.5 0.5 1]]
colorsarray=[2 3 4 6 7 8];
for k = 1:size(b,2)
    b(k).FaceColor = PetersColorMap(colorsarray(k));
end

ErrCol=2;
Errors=[PVC_AP(:,ErrCol)]
p=0.025;
Nerr=4;
tval=tinv(1-p,Nerr);
Errors=[PVC_APht(:,ErrCol)*tval]

%%%%%%%%%%%% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(A);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',A,Errors,'k','linestyle','none');
%%%%%%%%%%%%%%%


% AA=gca;
% A.YScale='log'
%  axis([0 600 0 20])

% xlabel('T (K)')
        ylabel('\chi')
        legtext={'AAA-1 Our Results','AAA-1 + 10%PS','Greenfield 2014','Experimental','Greenfield 2007'}
        legtext={'PE','PB','PP (Ewald)','PP (PPPM)','PS (Old)', 'PVC'}
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
