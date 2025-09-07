clc;close all;clear all;

%chi 1>>
phi =[0.1:0.001:0.9, log10(logspace(0.9,0.999999999,10000)) ]
k = 0.001985875;
T = 500;
N1_a = 4.87;
N2 = 1000;
X=[];
Y=[];
for chi = 0.00001:0.0001:0.30
%syms f(phi)
F = @(phi) ((phi.*log(phi))./N1_a+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
y = F(phi);

    
k=convhull(phi,y);
[phik,yk]=deal(phi(k),y(k));
j = diff(phik)>0;
[phik,yk]=deal(phik(j),yk(j));
interDist=vecnorm(  diff([phik(:),yk(:)],1,1) ,2,2);
[~,imax]=max(interDist);
x1=phik(imax); y1=yk(imax);      %the two tangent points (x1,y1) and (x2,y2)
x2=phik(imax+1); y2=yk(imax+1);
% %% Plot the results %%%
% plot(phi,y,'o-','linewidth',2);
% hold on
% pl=polyfit([x1,x2],[y1,y2],1);
% plot(phi,polyval(pl,phi))
% xlim([0,1])
% hold off
if isempty(x1)
    x1=NaN;
    x2=x1;
end
X=[X x1 x2];
Y=[Y chi chi];
end
%plot(X,Y,'o')
Constant_1=90.55;
T_1=Constant_1./Y;
%figure
%plot(X,T, 'o')
counter=1;
for i=1:length(X)
if isnan(X(counter)) || X(counter)<0.6
X(counter)=[];
T_1(counter)=[];
else
counter=counter+1;
end
end
D_1=[X;T_1]'
B_1=sortrows(D_1,1)


%chi_2>>
phi =[0.1:0.001:0.9, log10(logspace(0.9,0.999999999,10000)) ]
k = 0.001985875;
T = 500;
N1_b = 4.87;
N2 = 1000;
X=[];
Y=[];
for chi = 0.00001:0.0001:0.30;
%syms f(phi)
F = @(phi) ((phi.*log(phi))./N1_b+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
y = F(phi);

    
k=convhull(phi,y);
[phik,yk]=deal(phi(k),y(k));
j = diff(phik)>0;
[phik,yk]=deal(phik(j),yk(j));
interDist=vecnorm(  diff([phik(:),yk(:)],1,1) ,2,2);
[~,imax]=max(interDist);
x1=phik(imax); y1=yk(imax);      %the two tangent points (x1,y1) and (x2,y2)
x2=phik(imax+1); y2=yk(imax+1);
% %% Plot the results %%%
% plot(phi,y,'o-','linewidth',2);
% hold on
% pl=polyfit([x1,x2],[y1,y2],1);
% plot(phi,polyval(pl,phi))
% xlim([0,1])
% hold off
if isempty(x1)
    x1=NaN;
    x2=x1;
end
X=[X x1 x2];
Y=[Y chi chi];
end
%plot(X,Y,'o')
Constant_2=104.55;
T_2=Constant_2./Y;
%plot(X,T, 'o')
counter=1;
for i=1:length(X)
if isnan(X(counter)) || X(counter)<0.6
X(counter)=[];
T_2(counter)=[];
else
counter=counter+1;
end
end
D_2=[X;T_2]'
B_2=sortrows(D_2,1)



%chi_3>>
phi =[0.1:0.001:0.9, log10(logspace(0.9,0.999999999,10000)) ]
k = 0.001985875;
T = 500;
N1_c = 4.87;
N2 = 1000;
X=[];
Y=[];
for chi = 0.00001:0.0001:0.30;
%syms f(phi)
F = @(phi) ((phi.*log(phi))./N1_c+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
y = F(phi);

    
k=convhull(phi,y);
[phik,yk]=deal(phi(k),y(k));
j = diff(phik)>0;
[phik,yk]=deal(phik(j),yk(j));
interDist=vecnorm(  diff([phik(:),yk(:)],1,1) ,2,2);
[~,imax]=max(interDist);
x1=phik(imax); y1=yk(imax);      %the two tangent points (x1,y1) and (x2,y2)
x2=phik(imax+1); y2=yk(imax+1);
% %% Plot the results %%%
% plot(phi,y,'o-','linewidth',2);
% hold on
% pl=polyfit([x1,x2],[y1,y2],1);
% plot(phi,polyval(pl,phi))
% xlim([0,1])
% hold off
if isempty(x1)
    x1=NaN;
    x2=x1;
end
X=[X x1 x2];
Y=[Y chi chi];
end
%plot(X,Y,'o')
Constant_3=114.9;
T_3=Constant_3./Y;
%plot(X,T, 'o')
counter=1;
for i=1:length(X)
if isnan(X(counter)) || X(counter)<0.6
X(counter)=[];
T_3(counter)=[];
else
counter=counter+1;
end
end
D_3=[X;T_3]'
B_3=sortrows(D_3,1)




%chi_4>>
phi =[0.1:0.001:0.9, log10(logspace(0.9,0.999999999,10000)) ]
k = 0.001985875;
T = 500;
N1_d = 4.87;
N2 = 1000;
X=[];
Y=[];
for chi = 0.00001:0.0001:0.30;
%syms f(phi)
F = @(phi) ((phi.*log(phi))./N1_d+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
y = F(phi);

    
k=convhull(phi,y);
[phik,yk]=deal(phi(k),y(k));
j = diff(phik)>0;
[phik,yk]=deal(phik(j),yk(j));
interDist=vecnorm(  diff([phik(:),yk(:)],1,1) ,2,2);
[~,imax]=max(interDist);
x1=phik(imax); y1=yk(imax);      %the two tangent points (x1,y1) and (x2,y2)
x2=phik(imax+1); y2=yk(imax+1);
% %% Plot the results %%%
% plot(phi,y,'o-','linewidth',2);
% hold on
% pl=polyfit([x1,x2],[y1,y2],1);
% plot(phi,polyval(pl,phi))
% xlim([0,1])
% hold off
if isempty(x1)
    x1=NaN;
    x2=x1;
end
X=[X x1 x2];
Y=[Y chi chi];
end
%plot(X,Y,'o')
Constant_4=120.75;
T_4=Constant_4./Y;
%plot(X,T, 'o')
counter=1;
for i=1:length(X)
if isnan(X(counter)) || X(counter)<0.6
X(counter)=[];
T_4(counter)=[];
else
counter=counter+1;
end
end
D_4=[X;T_4]'
B_4=sortrows(D_4,1)



%chi_5>>
phi =[0.1:0.001:0.9, log10(logspace(0.9,0.999999999,10000)) ]
k = 0.001985875;
T = 500;
N1_d = 4.87;
N2 = 1000;
X=[];
Y=[];
for chi = 0.00001:0.0001:0.30;
%syms f(phi)
F = @(phi) ((phi.*log(phi))./N1_d+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
y = F(phi);

    
k=convhull(phi,y);
[phik,yk]=deal(phi(k),y(k));
j = diff(phik)>0;
[phik,yk]=deal(phik(j),yk(j));
interDist=vecnorm(  diff([phik(:),yk(:)],1,1) ,2,2);
[~,imax]=max(interDist);
x1=phik(imax); y1=yk(imax);      %the two tangent points (x1,y1) and (x2,y2)
x2=phik(imax+1); y2=yk(imax+1);
% %% Plot the results %%%
% plot(phi,y,'o-','linewidth',2);
% hold on
% pl=polyfit([x1,x2],[y1,y2],1);
% plot(phi,polyval(pl,phi))
% xlim([0,1])
% hold off
if isempty(x1)
    x1=NaN;
    x2=x1;
end
X=[X x1 x2];
Y=[Y chi chi];
end
%plot(X,Y,'o')
Constant_5=130.05;
T_5=Constant_5./Y;
%plot(X,T, 'o')
counter=1;
for i=1:length(X)
if isnan(X(counter)) || X(counter)<0.6
X(counter)=[];
T_5(counter)=[];
else
counter=counter+1;
end
end

D_5=[X;T_5]'
B_5=sortrows(D_5,1)

%plot(B(:,1),B(:,2),'-')
colors = {'black','red','green','blue','orange'}; 
%data=B;
figure
plot(1-B_1(:,1),B_1(:,2), 1-B_2(:,1),B_2(:,2), 1-B_3(:,1),B_3(:,2), 1-B_4(:,1),B_4(:,2), 1-B_5(:,1),B_5(:,2),'-')
colormap([0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 0 1]);
%plot(1-B_1(:,1),B_1(:,2),'-')
hold on 
xlabel('Monomer mole fraction')
        ylabel('Temperature (K)')
        legtext={'Day 0','Day 1','Day 5','Day 10', 'Day 30'}
        leg=legend(legtext,'Location','best')
        h1=gca;
        h1.XLim=[0 0.5]
        h1.YLim=[400 1200]
        set(gcf, 'units', 'inches', 'pos', [0 0 8.0947 6.7368]*.6)
        legend boxoff
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        set(leg,'FontSize',25);
        set(h1,'TickLength',[.02 .1])
        set(h1,'XMinorTick','on')
        %leg.Position=leg.Position
         export_fig 'bindoal_asphalt_AAA1' -png -r1-a4
