clc;close all;clear all;
figure
hold on

i

for N2=10:-1:4
    N2;
    % close all
    x1=0.0001:0.0001:0.9999;
    k=0.001985875;
    N1=20000*2;
    % N2=7;
    chi=0.1223;
    T0=500;
    zdW=chi*k*T0;
    % zdW=1.5887;
    % zdW=0.3972;
    options=optimset('Display','final');
    
    counter=1;
    Ts=[100:0.01:1000]
    for T=Ts
        if mod(T,100)==0
            T
        end
        %     T
        %find sign changes
        s=d2Spinodal(x1,N1,N2,zdW,T);
        s1=find(s(1:end-1)>=0 & s(2:end)<0);
        s2=find(s(1:end-1)<=0 & s(2:end)>0);
        if length(s1)+length(s2) >2
            error('too many sign changes?????')
        elseif length(s1)+length(s2)==0 %all mixed
            X1(counter)=NaN;
            X2(counter)=NaN;
        else
            interval=mean([x1(s1),x1(s2)]);
            
            try
                X1(counter)=fzero(@(x)d2Spinodal(x,N1,N2,zdW,T),[0.000001 interval],options);
                X2(counter)=fzero(@(x)d2Spinodal(x,N1,N2,zdW,T),[interval 0.9999999]);
            catch
                %     Test1=d2Spinodal([0.000001 0.49999999],N1,N2,zdW,T)
                %     Test=d2Spinodal(x1,N1,N2,zdW,T);
                %     figure
                %     plot(x1,Test)
                X1(counter)=NaN;
                X2(counter)=NaN;
            end
        end
        counter=counter+1;
    end
    
    % figure
    hold on
    plot([X1 fliplr(X2)],[Ts fliplr(Ts)],'-','Color',PetersColorMap(N2-3))
    % plot(X2,Ts,'k-')
    
    
    %
    %
    % Test=d2Spinodal(x1,N1,N2,zdW,320);
    %
    % figure
    % plot(x1,Test)
    
end


A=gca;
% A.YScale='log'
%  axis([0 5e6 -1e6 5e6 ])

xlabel('monomer mole fraction')
ylabel('T (K)')
legtext=fliplr({'N_{Binder}=4','N_{Binder}=5','N_{Binder}=6','N_{Binder}=7','N_{Binder}=8','N_{Binder}=9','N_{Binder}=10'})
leg=legend(legtext)
h1=gca;
set(gcf, 'units', 'inches', 'pos', [3 3 3.5 3])

legend boxon
set(gca,'FontSize',10)
set(gcf,'color','w');
set(gca,'color','None');
box on
set(leg,'FontSize',9);
set(h1,'TickLength',[.02 .1])
set(h1,'XMinorTick','on')
leg.Location='eastoutside'
%         export_fig 'Rt_example' -png -r800 -a1

function output=d2Spinodal(x1,N1,N2,zdW,T)
k=0.001985875;
chi=zdW/(k*T);
output=1./(N1*x1)+1./(N2*(1-x1))-2*chi;
output=output*(k*T);
end