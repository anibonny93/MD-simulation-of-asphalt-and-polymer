p=polyfit([0,1,2,3,0.5, 2.5],[0 0 0 0 -1, -10],4); %fake data
x=linspace(-1,4,10000);
y=polyval(p,x);
%% Calculation begins here %%%
k=convhull(x,y);
[xk,yk]=deal(x(k),y(k));
j = diff(xk)>0;
[xk,yk]=deal(xk(j),yk(j));
interDist=vecnorm(  diff([xk(:),yk(:)],1,1) ,2,2);
[~,imax]=max(interDist);
x1=xk(imax); y1=yk(imax);      %the two tangent points (x1,y1) and (x2,y2)
x2=xk(imax+1); y2=yk(imax+1);
%% Plot the results %%%
close all
hold on
plot(x,polyval(p,x))
 pl=polyfit([x1,x2],[y1,y2],1);
plot(x,polyval(pl,x))
xlim([0,3])
hold off