clc;close all;clear all;

%chi 1>>
phi =[0.1:0.001:0.9, log10(logspace(0.9,0.999999999,10000)) ]
k = 0.001985875;
T = 500;
N1 = 5;
N2 = 1000;
X=[];
Y=[];

for chi = 0.00001:0.0001:0.30;
%syms f(phi)
F = @(phi) ((phi.*log(phi))./N1+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
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
C=60;
T=C./Y
figure
%plot(X,T, 'o')

counter=1;
for i=1:length(X)
if isnan(X(counter)) || X(counter)<0.4
X(counter)=[];
T(counter)=[];
else
counter=counter+1;
end
end

D=[X;T]'
B=sortrows(D,1)

%plot(B(:,1),B(:,2),'-')

%data=B;
figure
plot(1-B(:,1),B(:,2),'-')
hold on 
%chi 1>>
phi =[0.1:0.001:0.9, log10(logspace(0.9,0.999999999,10000)) ]
k = 0.001985875;
T = 500;
N1 = 5;
N2 = 1000;
X=[];
Y=[];

for chi = 0.00001:0.0001:0.30;
%syms f(phi)
F = @(phi) ((phi.*log(phi))./N1+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
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
C=80;
T=C./Y
figure
%plot(X,T, 'o')

counter=1;
for i=1:length(X)
if isnan(X(counter)) || X(counter)<0.4
X(counter)=[];
T(counter)=[];
else
counter=counter+1;
end
end

D=[X;T]'
B=sortrows(D,1)

%plot(B(:,1),B(:,2),'-')

%data=B;
figure
plot(1-B(:,1),B(:,2),'-')
xlabel('Volume fraction')
        ylabel('Temperature (K)')
        %legtext={'Saturates','Asphaltenes','Aromatics','Resins'}
        %leg=legend(legtext,'Location','best')
        h1=gca;
        h1.XLim=[0 0.5]
        set(gcf, 'units', 'inches', 'pos', [0 0 8.0947 6.7368]*.6)
        legend boxoff
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        %set(leg,'FontSize',12);
        set(h1,'TickLength',[.02 .1])
        set(h1,'XMinorTick','on')
        leg.Position=leg.Position
         export_fig 'bindoal_test' -png -a4

