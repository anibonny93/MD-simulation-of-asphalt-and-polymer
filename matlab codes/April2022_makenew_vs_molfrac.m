close all
%% Readd
PE=[0.0115
0.0175
0.0055
0
0.1141
0.1962
0.0269
0.1192
0
0.1502
0.2294
0.0711
0
0.1969
0.6902
0.0404
0.0830
0.1056
0.0651
]

PB=[0.0809
0.0500
0.1118
0
0.1445
0.2722
0.0075
0.0859
0
0.2595
0.3449
0.1740
0
0.3186
1.1418
0.0920
0.1095
0.1677
0.0820
];

PP=[0.0188
0.0014
0.0363
0
0.0961
0.2086
0.0392
0.0404
0
0.2067
0.3025
0.1108
0
0.2468
0.8599
0.0772
0.0736
0.1392
0.0838
];

PS=[0.3390
0.2777
0.4010
0
0.1224
0.2098
0.0753
0.0821
0
0.1169
0.1309
0.1028
0
0.1906
0.5937
0.0416
0.2493
0.0357
0.0327
];


% PP=PP_PPPM;
% PP(14,:)=PP_Ewald(14,:);
% PP(19,:)=PP_Ewald(19,:);
% PP(10,1)=mean(PP(11:15,1));
% PP(17,1)=mean(PP(18:19,1));
% PE(1)=PE(2);
% PB(1)=PB(2);
% PS(1)=PS(2);
% PP(1)=PP(2);
%row1=400
%row2=500

types={
'Saturates'
'Squalane'
'Hopane'
'-'
'Asphaltenes'
'Asphaltene Phenol'
'Asphaltene Pyrrole'
'Asphaltene Thiophene'
' -- '
'Aromatics'
'Perhydrophenanthrene Napththalene'
'Dioctylcyclohexane Naphthalene'
'  ---  '
'Resins'
'Benzobisbenzothiophene'
'Pyridinohopane'
'Quinolinohopane'
'Thioisorenieratane'
'Trimethylbenzeneoxane'
};

X=categorical(types);
X=reordercats(X,types);

x_asph=0.1725;
x_aro=0.1595;
x_res=0.3973;
x_sat=0.1112;

x=0:0.1:1;

x_array=repmat([x_sat;x_asph;x_aro;x_res],1,length(x));

%% Polyethylene

% Saturates
x_sat_array=x_array;
x_sat_array(1,:)=x;
x_sat_array=x_sat_array./(repmat(sum(x_sat_array,1),4,1));
chis_sat=[PE(3,1);PE(8,1);PE(12,1);PE(19,1)];
comb_chi_sat=x_sat_array'*chis_sat;
% Aspahltenes
x_asph_array=x_array;
x_asph_array(2,:)=x;
x_asph_array=x_asph_array./(repmat(sum(x_asph_array,1),4,1));
chis_asph=[PE(3,1);PE(8,1);PE(12,1);PE(19,1)];
comb_chi_asph=x_asph_array'*chis_asph;
% Aromatics
x_aro_array=x_array;
x_aro_array(3,:)=x;
x_aro_array=x_aro_array./(repmat(sum(x_aro_array,1),4,1));
chis_aro=[PE(3,1);PE(8,1);PE(12,1);PE(19,1)];
comb_chi_aro=x_aro_array'*chis_aro;
% Resins
x_res_array=x_array;
x_res_array(4,:)=x;
x_res_array=x_res_array./(repmat(sum(x_res_array,1),4,1));
chis_res=[PE(3,1);PE(8,1);PE(12,1);PE(19,1)];
comb_chi_res=x_res_array'*chis_res;

figure
hold on
plot(x_sat_array(1,:),comb_chi_sat,'k-','LineWidth',2)
plot(x_asph_array(2,:),comb_chi_asph,'b-','LineWidth',2)
plot(x_aro_array(3,:),comb_chi_aro,'r-','LineWidth',2)
plot(x_res_array(4,:),comb_chi_res,'g-','LineWidth',2)


xlabel('monomer mole fraction')
        ylabel('\chi_{PE}')
        legtext={'Saturates','Asphaltenes','Aromatics','Resins'}
        leg=legend(legtext,'Location','best')
        h1=gca;
        h1.XLim=[0 0.5]
        set(gcf, 'units', 'inches', 'pos', [0 0 8.0947 6.7368]*.6)
        legend boxoff
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        set(leg,'FontSize',12);
        set(h1,'TickLength',[.02 .1])
        set(h1,'XMinorTick','on')
        leg.Position=leg.Position
         export_fig 'PE_mol_frac' -png -r1000 -a4

%% Polybutadienne
% Saturates
x_sat_array=x_array;
x_sat_array(1,:)=x;
x_sat_array=x_sat_array./(repmat(sum(x_sat_array,1),4,1));
chis_sat=[PB(3,1);PB(8,1);PB(12,1);PB(19,1)];
comb_chi_sat=x_sat_array'*chis_sat;
% Aspahltenes
x_asph_array=x_array;
x_asph_array(2,:)=x;
x_asph_array=x_asph_array./(repmat(sum(x_asph_array,1),4,1));
chis_asph=[PB(3,1);PB(8,1);PB(12,1);PB(19,1)];
comb_chi_asph=x_asph_array'*chis_asph;
% Aromatics
x_aro_array=x_array;
x_aro_array(3,:)=x;
x_aro_array=x_aro_array./(repmat(sum(x_aro_array,1),4,1));
chis_aro=[PB(3,1);PB(8,1);PB(12,1);PB(19,1)];
comb_chi_aro=x_aro_array'*chis_aro;
% Resins
x_res_array=x_array;
x_res_array(4,:)=x;
x_res_array=x_res_array./(repmat(sum(x_res_array,1),4,1));
chis_res=[PB(3,1);PB(8,1);PB(12,1);PB(19,1)];
comb_chi_res=x_res_array'*chis_res;


figure
hold on
plot(x_sat_array(1,:),comb_chi_sat,'k-','LineWidth',2)
plot(x_asph_array(2,:),comb_chi_asph,'b-','LineWidth',2)
plot(x_aro_array(3,:),comb_chi_aro,'r-','LineWidth',2)
plot(x_res_array(4,:),comb_chi_res,'g-','LineWidth',2)


xlabel('monomer mole fraction')
        ylabel('\chi_{PB}')
        legtext={'Saturates','Asphaltenes','Aromatics','Resins'}
        leg=legend(legtext,'Location','best')
        h1=gca;   
        h1.XLim=[0 0.5]
        set(gcf, 'units', 'inches', 'pos', [0 0 8.0947 6.7368]*.6)
        legend boxoff
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        set(leg,'FontSize',12);
        set(h1,'TickLength',[.02 .1])
        set(h1,'XMinorTick','on')
leg.Position=leg.Position
         export_fig 'PB_mol_frac' -png -r1000 -a4
         
         
%% Polypropylene
% Saturates
x_sat_array=x_array;
x_sat_array(1,:)=x;
x_sat_array=x_sat_array./(repmat(sum(x_sat_array,1),4,1));
chis_sat=[PP(3,1);PP(8,1);PP(12,1);PP(19,1)];
comb_chi_sat=x_sat_array'*chis_sat;
% Aspahltenes
x_asph_array=x_array;
x_asph_array(2,:)=x;
x_asph_array=x_asph_array./(repmat(sum(x_asph_array,1),4,1));
chis_asph=[PP(3,1);PP(8,1);PP(12,1);PP(19,1)];
comb_chi_asph=x_asph_array'*chis_asph;
% Aromatics
x_aro_array=x_array;
x_aro_array(3,:)=x;
x_aro_array=x_aro_array./(repmat(sum(x_aro_array,1),4,1));
chis_aro=[PP(3,1);PP(8,1);PP(12,1);PP(19,1)];
comb_chi_aro=x_aro_array'*chis_aro;
% Resins
x_res_array=x_array;
x_res_array(4,:)=x;
x_res_array=x_res_array./(repmat(sum(x_res_array,1),4,1));
chis_res=[PP(3,1);PP(8,1);PP(12,1);PP(19,1)];
comb_chi_res=x_res_array'*chis_res;


figure
hold on
plot(x_sat_array(1,:),comb_chi_sat,'k-','LineWidth',2)
plot(x_asph_array(2,:),comb_chi_asph,'b-','LineWidth',2)
plot(x_aro_array(3,:),comb_chi_aro,'r-','LineWidth',2)
plot(x_res_array(4,:),comb_chi_res,'g-','LineWidth',2)

xlabel('monomer mole fraction')
        ylabel('\chi_{PP}')
        legtext={'Saturates','Asphaltenes','Aromatics','Resins'}
        leg=legend(legtext,'Location','best')
        h1=gca;
        h1.XLim=[0 0.5]
        set(gcf, 'units', 'inches', 'pos', [0 0 8.0947 6.7368]*.6)
        legend boxoff
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        set(leg,'FontSize',12);
        set(h1,'TickLength',[.02 .1])
        set(h1,'XMinorTick','on')
        leg.Position=leg.Position
         export_fig 'PP_mol_frac' -png -r1000 -a4
         
         
%% Polystyrnene

% Saturates
x_sat_array=x_array;
x_sat_array(1,:)=x;
x_sat_array=x_sat_array./(repmat(sum(x_sat_array,1),4,1));
chis_sat=[PS(3,1);PS(8,1);PS(12,1);PS(19,1)];
comb_chi_sat=x_sat_array'*chis_sat;
% Aspahltenes
x_asph_array=x_array;
x_asph_array(2,:)=x;
x_asph_array=x_asph_array./(repmat(sum(x_asph_array,1),4,1));
chis_asph=[PS(3,1);PS(8,1);PS(12,1);PS(19,1)];
comb_chi_asph=x_asph_array'*chis_asph;
% Aromatics
x_aro_array=x_array;
x_aro_array(3,:)=x;
x_aro_array=x_aro_array./(repmat(sum(x_aro_array,1),4,1));
chis_aro=[PS(3,1);PS(8,1);PS(12,1);PS(19,1)];
comb_chi_aro=x_aro_array'*chis_aro;
% Resins
x_res_array=x_array;
x_res_array(4,:)=x;
x_res_array=x_res_array./(repmat(sum(x_res_array,1),4,1));
chis_res=[PS(3,1);PS(8,1);PS(12,1);PS(19,1)];
comb_chi_res=x_res_array'*chis_res;


figure
hold on
plot(x_sat_array(1,:),comb_chi_sat,'k-','LineWidth',2)
plot(x_asph_array(2,:),comb_chi_asph,'b-','LineWidth',2)
plot(x_aro_array(3,:),comb_chi_aro,'r-','LineWidth',2)
plot(x_res_array(4,:),comb_chi_res,'g-','LineWidth',2)

xlabel('monomer mole fraction')
        ylabel('\chi_{PS}')
        legtext={'Saturates','Asphaltenes','Aromatics','Resins'}
        leg=legend(legtext,'Location','best')
        h1=gca;
        h1.XLim=[0 0.5]
        set(gcf, 'units', 'inches', 'pos', [0 0 8.0947 6.7368]*.6)
        legend boxoff
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        set(leg,'FontSize',12);
        set(h1,'TickLength',[.02 .1])
        set(h1,'XMinorTick','on')
leg.Position=leg.Position
         %export_fig 'PS_mol_frac' -png -a4
         export_fig 'PS_mol_frac' -png -r1000 -a4
