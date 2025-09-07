clc;close all;clear all;
phi=0:.001:1;
k = 0.001985875;
T = 500;
N1 = 1;
N2=2;
chi =2;

f= @(phi) ((phi./N1).*log(phi)+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
y=f(phi);
plot(phi,y,'linewidth',2);
hold on

counter=1;
phi1= 0:0.01:0.3;
curve1=zeros(numel(phi1),2);
for phi1= phi1
    f= @(phi1) ((phi1./N1).*log(phi1)+((1-phi1)./N2).*log(1-phi1)+chi.*phi1.*(1-phi1))
    y1=f(phi1)
    curve1(counter,:)=[phi1 y1];
    counter=counter+1;
end
curve1=curve1';

counter=1;
phi2=0.7:0.01:1;
curve2=zeros(numel(phi2),2);
for phi2= phi2
    f= @(phi2) ((phi2./N1).*log(phi2)+((1-phi2)./N2).*log(1-phi2)+chi.*phi2.*(1-phi2))
    y1=f(phi2)
    curve2(counter,:)=[phi2 y1];
    counter=counter+1;
end
curve2=curve2';

P=InterX(curve1,curve2);




