
%phi=0:.001:1;
k = 0.001985875;
T = 500;
N1 = 1;
N2=2;
chi=2;
x=0:0.001:1;
syms x f(x) z(x)

f(x)=((x./N1).*log(x)+((1-x)./N2).*log(1-x)+chi.*x.*(1-x));
y=f(x);



k=convhull(x,y);
[xk, yk] = deal(x(k),y(k));
j = diff (xk)>0;
[xk, yk] = deal(xk(j),yk(j));
interdist=vecnorm(diff([xk(:),yk(:)],1,1) ,2,2);
[~,imax]=max(interdist);
x1=xk(imax); y1=yk(imax);
x2=xk(imax+1);y2=yk(imax+1);

close all
hold on
plot(x,polyval(p,x))
p1=polyfit([x1,x2],[y1,y2]);
plot(x,polyval(p1,x0))
xlim([0,3])
hold off


% df=diff(f,phi);
% 
% phi1=0.1:0.01:0.2;
% m=df(phi1);
% z1=f(phi1);
% z(phi)= z1+(m.*(phi-phi1));
% figure
% phi=0:0.01:1;
% 
% plot(phi,f(phi))
% hold on 
% plot(phi,z(phi))


%phi=0:.001:1;
%plot(phi,f(phi))




% 
% phi1=0:0.01:0.2;
% m=df(phi1)
% z1=f(phi1)
% z(phi)=z1+(m.*(phi-phi1))
% figure
% x=
% 
% 
% y=f(phi);
% plot(phi,y,'linewidth',2);
% hold on