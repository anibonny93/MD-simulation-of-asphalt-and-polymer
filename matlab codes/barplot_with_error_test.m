%
%X_plot=1:19;
%PE500=[0	0.184 0.082 0.1482 0 0 0.048 0.190 0 0 0.4307 0.0341 0.107 0.158 0.069 0 0 -0.0419 0.042]';
%PE500_err_high=[0 0.202  0.1098 0.1770 0 0 0.077 0.257 0 0 0.485 0.070 0.1302 0.174  0.09116 0 0 -0.0146 0.069];
%PE500_err_low=[0 0.1666 0.0557 0.1194 0 0 0.0189 0.1229 0 0 0.375 -0.00185 0.0844 0.142 0.0478 0 0 -0.0692 0.0148];

%bar(X_plot,PE500)

%hold on 

%er = errorbar(X_plot,PE500,PE500_err_high,PE500_err_low);

%hold off 

PE500=[0	0.184 0.082 0.1482 0 0 0.048 0.190 0 0 0.4307 0.0341 0.107 0.158 0.069 0 0 -0.0419 0.042];  %% mean values
SE_low=[0 0.0144 0.0217 0.0232 0 0 0.0237 0.0541 0 0  0.04411 0.0289 0.0184 0.01296  0.01744 0 0 0.0220 0.0218]*tinv(0.975,4); %% standard error
SE_high=SE_low;

PE500(1)=mean(PE500(2:4));
PE500(6)=mean(PE500(7:8));
PE500(10)=mean(PE500(11:15));
PE500(17)=mean(PE500(18:19));

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

bar(X,PE500)
hold on 
er = errorbar(X,PE500,SE_low,SE_high);
hold off
