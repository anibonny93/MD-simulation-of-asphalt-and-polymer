clc;close all;clear all;

phi=0:.001:1;
k = 0.001985875;
T = 500;
N1 = 1;
N2=2;
chi =2;

f= @(phi) ((phi./N1).*log(phi)+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
y=f(phi);

dphi=mean(diff(phi));
dy=gradient(y,dphi);

for phiq1=find((phi>=0) & (phi<=0.2))
    slope_phi1=phi(phiq1);
    slope_y1=dy(phiq1);
    phiq2=find((phi>=0.8) & (phi<=1));
    slope_phi2=phi(phiq2);
    slope_y2=dy(phiq2);
    if slope_y1 == slope_y2
        plot ([slope_phi1,slope_y1], [slope_phi2,slope_y2])
    end
end
figure
plot(phi,dy,'-g')
hold on 
plot (slope_phi1, slope_y1, '-r', 'LineWidth',1)
plot (slope_phi2, slope_y2, '-r', 'LineWidth',1)
hold off
grid
legend ('Slope Of All Data', 'Slope of Data In Region Of Interest', 'Location', 'W')