clc;close all;clear all;
phi =0.001:.001:0.999;
k = 0.001985875;
%T = 500;
%N1 = 10;
%N2 = 200;
N=2000;
chi = 0.01;
B=0.015;
A=0.025;
%syms f(phi)
T = @(phi) B./((log(phi./(1-phi))./((2.*phi-1).*N)))-A;
y = T(phi);

    
k=convhull(phi,y);
[phik,yk]=deal(phi(k),y(k));
j = diff(phik)>0;
[phik,yk]=deal(phik(j),yk(j));
interDist=vecnorm(  diff([phik(:),yk(:)],1,1) ,2,2);
[~,imax]=max(interDist);
x1=phik(imax); y1=yk(imax);      %the two tangent points (x1,y1) and (x2,y2)
x2=phik(imax+1); y2=yk(imax+1);
%% Plot the results %%%
plot(phi,y,'linewidth',2);
hold on
pl=polyfit([x1,x2],[y1,y2],1);
plot(phi,polyval(pl,phi))
xlim([0,1])
hold off
