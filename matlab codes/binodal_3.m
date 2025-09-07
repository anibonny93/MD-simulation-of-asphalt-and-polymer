clc;close all;clear all;
phi=0:.001:1;
k = 0.001985875;
T = 500;
N1 = 1;
N2=2;
chi =2;
%syms f(phi)
f= @(phi) ((phi./N1).*log(phi)+((1-phi)./N2).*log(1-phi)+chi.*phi.*(1-phi));
y=f(phi);
plot(phi,y,'linewidth',2);
hold on

for phi1 =0:0.001:0.2
    h= 0.000001;
    slope=(f(phi1+h)-f(phi1))./(h);
    tangent=slope.*(phi-phi1)+f(phi1);
%[v,w] =unique(tangent,'stable');
%duplicate_indices=setdiff(1:numel(tangent),w)
    %plot(phi,tangent,'linewidth',2);
    
end

for phi2 =0.8:0.001:1
    h= 0.000001;
    slope=(f(phi2+h)-f(phi2))./(h);
    tangent2=slope.*(phi-phi2)+f(phi2);
%[v,w] =unique(tangent,'stable');
%duplicate_indices=setdiff(1:numel(tangent),w)
    %plot(phi,tangent2,'linewidth',2);
end

for phi1=0.145695
    for phi2=0.9353
        tangent_m=f(phi1)+((f(phi1)-f(phi2))./(phi1-phi2)).*(phi-phi1);
    end
end
    plot(phi,tangent_m)

%     phi2=0.8:0.01:1;
% slope2=(f(phi2+h)-f(phi2))./(h)
% tangent2=slope2*(phi-phi2)+f(phi2);
% plot(phi,tangent2,'linewidth',2);
% 
% m=diff(f);
% phi1=0:0.1:2;
% tangent=m.*(phi-phi1)+f(phi1);
% 
% 
% phi=0:0.001:0.2;
% df(phi);
% f(phi);
% plot(phi,f(phi),df(phi),'r*');