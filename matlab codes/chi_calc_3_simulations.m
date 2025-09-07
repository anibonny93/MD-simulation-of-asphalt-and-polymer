close all 
clear all

path_p=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/PS10_new/PS10_']; %%%Change for polymer

path_b=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/Oxidized_asphalts/TIR_oxidized/TIR_']; %%%Change for binder

path_mix=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/Oxidized_asphalt_PS/TIR_PS_oxidized/TIR_PS_oxidized_']; %%%Change for mix

rep=256;%input('Enter replication number:')
mpm=10;%input('Enter monomer per molecule for polymer:')
temperature=500;%input('Enter the temperature:')
in=3;%input('Enter the run number')

[mass_p] = ReadMassFromLAMMPSData([path_p '1_replicated.data'])
mass_p=mass_p*rep;
dcool_p=[temperature]';
Ucool_p=dcool_p;

[mass_b] = ReadMassFromLAMMPSData([path_b '1_replicated.data']) %'_b' denotes binder
mass_b=mass_b*rep;
dcool_b=[temperature]';
Ucool_b=dcool_b;
 
[mass_mix] = ReadMassFromLAMMPSData([path_mix '1_replicated.data'])
mass_mix=mass_mix*rep;
dcool_mix=[temperature];
Ucool_mix=dcool_mix;

chiall_1=0;
chiall_2=0;
chiall_3=0;
chiall_4=0;
chiall_5=0;

timeframe=800:10:5400;%input('put the timeframe:')


for i=1:in
    filename_p=[path_p num2str(i) '.log'];
    eval(['[R' num2str(i) '] = readlammpslog(filename_p)'])
    filename_b=[path_b num2str(i) '.log'];
    eval(['[B' num2str(i) '] = readlammpslog(filename_b)'])
    filename_mix=[path_mix num2str(i) '.log'];
    eval(['[M' num2str(i) '] = readlammpslog(filename_mix)'])
end

twidth=100;
counter=0;

for tstart=timeframe %loops over time
    tend=tstart+twidth;
    counter=counter+1;
    
    %This section is for the polymer
    
        dcool_p(2)=mean(R1.Volume(tstart:tend));
        Ucool_p(2)=mean(R1.PotEng(tstart:tend));
        dcool_p(3)=mean(R2.Volume(tstart:tend));
        Ucool_p(3)=mean(R2.PotEng(tstart:tend));
        dcool_p(4)=mean(R3.Volume(tstart:tend));
        Ucool_p(4)=mean(R3.PotEng(tstart:tend));
%        dcool_p(5)=mean(R4.Volume(tstart:tend));
 %       Ucool_p(5)=mean(R4.PotEng(tstart:tend));
  %      dcool_p(6)=mean(R5.Volume(tstart:tend));
   %     Ucool_p(6)=mean(R5.PotEng(tstart:tend));
        

        dcool_p;
        Ucool_p;
        
    dcool_p(:,2:end)=((dcool_p(:,2:end)/mass_p).^-1)*1.661; %'_p' signifies polymer
    dall_p=flipud(dcool_p);
    %dall_p(:,end+1)=mean(dall_p(:,2:end)')';
    %dall_p(:,end+1)=std(dall_p(:,2:end-1)')';
    %SE=dall_p(:,end)/sqrt(N);
    %dall_p(:,end+1)=-SE*tinv(0.975,N-1)+dall_p(:,N+2);
    %dall_p(:,end+1)=SE*tinv(0.975,N-1)+dall_p(:,N+2);

    %Uall_p=flipud(Ucool_p);
    %Uall_p(:,end+1)=mean(Uall_p(:,2:end)')';
    %Uall_p(:,end+1)=std(Uall_p(:,2:end-1)')';
    %SE=Uall_p(:,end)/sqrt(N);
    %Uall_p(:,end+1)=-SE*tinv(0.975,N-1)+Uall_p(:,N+2);
    %Uall_p(:,end+1)=SE*tinv(0.975,N-1)+Uall_p(:,N+2);

    %dall_p
    %Uall_p

    %dmean_p=dcool_p(2);
    %Umean_p=Ucool_p(2);
    mpm_p=[mpm]; %monomer per molecule for polymer
    vmpm(1)= (mass_p/rep)./(mpm_p.*dcool_p(2));%volume per mole of monomer for T 500K
    vmpm(2)= (mass_p/rep)./(mpm_p.*dcool_p(3));
    vmpm(3)= (mass_p/rep)./(mpm_p.*dcool_p(4));
 %   vmpm(4)= (mass_p/rep)./(mpm_p.*dcool_p(5));
  %  vmpm(5)= (mass_p/rep)./(mpm_p.*dcool_p(6));
    tm_p=rep.*mpm_p;%total number of monomers for polymer for T 500K

    
    %This section is for the binder

        dcool_b(2)=mean(B1.Volume(tstart:tend));
        Ucool_b(2)=mean(B1.PotEng(tstart:tend));
        dcool_b(3)=mean(B2.Volume(tstart:tend));
        Ucool_b(3)=mean(B2.PotEng(tstart:tend));
        dcool_b(4)=mean(B3.Volume(tstart:tend));
        Ucool_b(4)=mean(B3.PotEng(tstart:tend));
  %      dcool_b(5)=mean(B4.Volume(tstart:tend));
   %     Ucool_b(5)=mean(B4.PotEng(tstart:tend));
    %    dcool_b(6)=mean(B5.Volume(tstart:tend));
     %   Ucool_b(6)=mean(B5.PotEng(tstart:tend));
    
    dcool_b(:,2:end)=((dcool_b(:,2:end)/mass_b).^-1)*1.661;
    dall_b=flipud(dcool_b);
    
    %dall_b(:,end+1)=mean(dall_b(:,2:end)')';
    %%=dall_b(:,end)/sqrt(N);
    %dall_b(:,end+1)=-SE*tinv(0.975,N-1)+dall_b(:,N+2);
    %dall_b(:,end+1)=SE*tinv(0.975,N-1)+dall_b(:,N+2);

    %Uall_b=flipud(Ucool_b);
    %Uall_b(:,end+1)=mean(Uall_b(:,2:end)')';
    %Uall_b(:,end+1)=std(Uall_b(:,2:end-1)')';
    %SE=Uall_b(:,end)/sqrt(N);
    %%Uall_b(:,end+1)=-SE*tinv(0.975,N-1)+Uall_b(:,N+2);
    %Uall_b(:,end+1)=SE*tinv(0.975,N-1)+Uall_b(:,N+2);


% 
%     dcool_b
%     Ucool_b
% 
       %dmean_b=dcool_b(2)
       %Umean_b=Ucool_b(2)
       
        mpm_b(1)=(mass_b/rep)./(vmpm(1).*dcool_b(2));
        tm_b(1)=rep.*mpm_b(1);%total number of monomers for binder for T 500K 
        mpm_b(2)=(mass_b/rep)./(vmpm(2).*dcool_b(3));
        tm_b(2)=rep.*mpm_b(2);                    
        mpm_b(3)=(mass_b/rep)./(vmpm(3).*dcool_b(4));
        tm_b(3)=rep.*mpm_b(3);   
  %      mpm_b(4)=(mass_b/rep)./(vmpm(4).*dcool_b(5));
   %     tm_b(4)=rep.*mpm_b(4);   
    %    mpm_b(5)=(mass_b/rep)./(vmpm(5).*dcool_b(6));
     %   tm_b(5)=rep.*mpm_b(5);


     %This section is for the mix
     

        dcool_mix(2)=mean(M1.Volume(tstart:tend));
        Ucool_mix(2)=mean(M1.PotEng(tstart:tend));
        dcool_mix(3)=mean(M2.Volume(tstart:tend));
        Ucool_mix(3)=mean(M2.PotEng(tstart:tend));
        dcool_mix(4)=mean(M3.Volume(tstart:tend));
        Ucool_mix(4)=mean(M3.PotEng(tstart:tend));
 %       dcool_mix(5)=mean(M4.Volume(tstart:tend));
  %      Ucool_mix(5)=mean(M4.PotEng(tstart:tend));
   %     dcool_mix(6)=mean(M5.Volume(tstart:tend));
    %    Ucool_mix(6)=mean(M5.PotEng(tstart:tend));
  

    dcool_mix(:,2:end)=((dcool_mix(:,2:end)/mass_mix).^-1)*1.661;
    dall_mix=flipud(dcool_mix);
    %dall_mix(:,end+1)=mean(dall_mix(:,2:end)')';
    %dall_mix(:,end+1)=std(dall_mix(:,2:end-1)')';
    %SE=dall_mix(:,end)/sqrt(N);
    %dall_mix(:,end+1)=-SE*tinv(0.975,N-1)+dall_mix(:,N+2);
    %dall_mix(:,end+1)=SE*tinv(0.975,N-1)+dall_mix(:,N+2);

    %Uall_mix=flipud(Ucool_mix);
    %Uall_mix(:,end+1)=mean(Uall_mix(:,2:end)')';
    %Uall_mix(:,end+1)=std(Uall_mix(:,2:end-1)')';
    %SE=Uall_mix(:,end)/sqrt(N);
    %Uall_mix(:,end+1)=-SE*tinv(0.975,N-1)+Uall_mix(:,N+2);
    %Uall_mix(:,end+1)=SE*tinv(0.975,N-1)+Uall_mix(:,N+2);


%     dcool_mix
%     Ucool_mix

%     dmean_mix=dcool_mix(2)
%     Umean_mix=Ucool_mix(2)

%for mix 1

    Umix_pm(1) = Ucool_mix(2)./(tm_b(1)+tm_p);
    U_p_pm(1) = Ucool_p(2)./tm_p;
    U_b_pm(1) = Ucool_b(2)./tm_b(1);
    mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b(1))./mpm_b(1).^2);
    mass_1=((mass_p/rep)*tm_p)./mpm_p.^2;
    mass_2(1)=((mass_b/rep)*tm_b(1))./mpm_b(1).^2;
    
    PV_m(1)=(1./dcool_mix(2))*(1/1.661).*mass_mix_m;
    PV_1(1)=(1./dcool_p(2))*(1/1.661).*mass_1;
    PV_2(1)=(1./dcool_b(2))*(1/1.661).*mass_2(1);
    PV_mix_pm(1)=(0.00001458*PV_m(1))./(tm_p+tm_b(1));
    PV_1_pm(1)=(0.00001458*PV_1(1))./(tm_p);
    PV_2_pm(1)=(0.00001458*PV_2(1))./(tm_b(1));
    phi1(1)=(PV_1_pm(1))./(PV_1_pm(1)+PV_2_pm(1));
    phi2(1)=1-phi1(1);
    deltaU(1)=(Umix_pm(1).*(tm_p+tm_b(1))-(U_p_pm(1).*tm_p)-(U_b_pm(1).*tm_b(1)))./(tm_b(1)+tm_p);
    deltaPV(1)=(PV_mix_pm(1).*(tm_p+tm_b(1))-PV_1_pm(1).*tm_p-PV_2_pm(1).*tm_b(1));
    deltaH(1)=deltaU(1)+deltaPV(1);
    chi(1)=deltaH(1)./(phi1(1).*phi2(1)*temperature*0.001985875);
    
%for mix 2
    Umix_pm(2) = Ucool_mix(3)./(tm_b(2)+tm_p);
    U_p_pm(2) = Ucool_p(3)./tm_p;
    U_b_pm(2) = Ucool_b(3)./tm_b(2);
    mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b(2))./mpm_b(2).^2);
    mass_1=((mass_p/rep)*tm_p)./mpm_p.^2;
    mass_2(2)=((mass_b/rep)*tm_b(2))./mpm_b(2).^2;
    
    PV_m(2)=(1./dcool_mix(3))*(1/1.661).*mass_mix_m;
    PV_1(2)=(1./dcool_p(3))*(1/1.661).*mass_1;
    PV_2(2)=(1./dcool_b(3))*(1/1.661).*mass_2(2);
    PV_mix_pm(2)=(0.00001458*PV_m(2))./(tm_p+tm_b(2));
    PV_1_pm(2)=(0.00001458*PV_1(2))./(tm_p);
    PV_2_pm(2)=(0.00001458*PV_2(2))./(tm_b(2));
    phi1(2)=(PV_1_pm(2))./(PV_1_pm(2)+PV_2_pm(2));
    phi2(2)=1-phi1(2);
    deltaU(2)=(Umix_pm(2).*(tm_p+tm_b(2))-(U_p_pm(2).*tm_p)-(U_b_pm(2).*tm_b(2)))./(tm_b(2)+tm_p);
    deltaPV(2)=(PV_mix_pm(2).*(tm_p+tm_b(2))-PV_1_pm(2).*tm_p-PV_2_pm(2).*tm_b(2));
    deltaH(2)=deltaU(2)+deltaPV(2);
    chi(2)=deltaH(2)./(phi1(2).*phi2(2)*temperature*0.001985875);
    
%for mix 3
    Umix_pm(3) = Ucool_mix(4)./(tm_b(3)+tm_p);
    U_p_pm(3) = Ucool_p(4)./tm_p;
    U_b_pm(3) = Ucool_b(4)./tm_b(3);
    mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b(3))./mpm_b(3).^2);
    mass_1=((mass_p/rep)*tm_p)./mpm_p.^2;
    mass_2(3)=((mass_b/rep)*tm_b(3))./mpm_b(3).^2;
    
    PV_m(3)=(1./dcool_mix(4))*(1/1.661).*mass_mix_m;
    PV_1(3)=(1./dcool_p(4))*(1/1.661).*mass_1;
    PV_2(3)=(1./dcool_b(4))*(1/1.661).*mass_2(3);
    PV_mix_pm(3)=(0.00001458*PV_m(3))./(tm_p+tm_b(3));
    PV_1_pm(3)=(0.00001458*PV_1(3))./(tm_p);
    PV_2_pm(3)=(0.00001458*PV_2(3))./(tm_b(3));
    phi1(3)=(PV_1_pm(3))./(PV_1_pm(3)+PV_2_pm(3));
    phi2(3)=1-phi1(3);
    deltaU(3)=(Umix_pm(3).*(tm_p+tm_b(3))-(U_p_pm(3).*tm_p)-(U_b_pm(3).*tm_b(3)))./(tm_b(3)+tm_p);
    deltaPV(3)=(PV_mix_pm(3).*(tm_p+tm_b(3))-PV_1_pm(3).*tm_p-PV_2_pm(3).*tm_b(3));
    deltaH(3)=deltaU(3)+deltaPV(3);
    chi(3)=deltaH(3)./(phi1(3).*phi2(3)*temperature*0.001985875);

% %for mix 4
%     Umix_pm(4) = Ucool_mix(5)./(tm_b(4)+tm_p);
%     U_p_pm(4) = Ucool_p(5)./tm_p;
%     U_b_pm(4) = Ucool_b(5)./tm_b(4);
%     mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b(4))./mpm_b(4).^2);
%     mass_1=((mass_p/rep)*tm_p)./mpm_p.^2;
%     mass_2(4)=((mass_b/rep)*tm_b(4))./mpm_b(4).^2;
%     
%     PV_m(4)=(1./dcool_mix(5))*(1/1.661).*mass_mix_m;
%     PV_1(4)=(1./dcool_p(5))*(1/1.661).*mass_1;
%     PV_2(4)=(1./dcool_b(5))*(1/1.661).*mass_2(4);
%     PV_mix_pm(4)=(0.00001458*PV_m(4))./(tm_p+tm_b(4));
%     PV_1_pm(4)=(0.00001458*PV_1(4))./(tm_p);
%     PV_2_pm(4)=(0.00001458*PV_2(4))./(tm_b(4));
%     phi1(4)=(PV_1_pm(4))./(PV_1_pm(4)+PV_2_pm(4));
%     phi2(4)=1-phi1(4);
%     deltaU(4)=(Umix_pm(4).*(tm_p+tm_b(4))-(U_p_pm(4).*tm_p)-(U_b_pm(4).*tm_b(4)))./(tm_b(4)+tm_p);
%     deltaPV(4)=(PV_mix_pm(4).*(tm_p+tm_b(4))-PV_1_pm(4).*tm_p-PV_2_pm(4).*tm_b(4));
%     deltaH(4)=deltaU(4)+deltaPV(4);
%     chi(4)=deltaH(4)./(phi1(4).*phi2(4)*temperature*0.001985875);
%  
%  %for mix 5
%     Umix_pm(5) = Ucool_mix(6)./(tm_b(5)+tm_p);
%     U_p_pm(5) = Ucool_p(6)./tm_p;
%     U_b_pm(5) = Ucool_b(6)./tm_b(5);
%     mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b(5))./mpm_b(5).^2);
%     mass_1=((mass_p/rep)*tm_p)./mpm_p.^2;
%     mass_2(5)=((mass_b/rep)*tm_b(5))./mpm_b(5).^2;
%     
%     PV_m(5)=(1./dcool_mix(6))*(1/1.661).*mass_mix_m;
%     PV_1(5)=(1./dcool_p(6))*(1/1.661).*mass_1;
%     PV_2(5)=(1./dcool_b(6))*(1/1.661).*mass_2(5);
%     PV_mix_pm(5)=(0.00001458*PV_m(5))./(tm_p+tm_b(5));
%     PV_1_pm(5)=(0.00001458*PV_1(5))./(tm_p);
%     PV_2_pm(5)=(0.00001458*PV_2(5))./(tm_b(5));
%     phi1(5)=(PV_1_pm(5))./(PV_1_pm(5)+PV_2_pm(5));
%     phi2(5)=1-phi1(5);
%     deltaU(5)=(Umix_pm(5).*(tm_p+tm_b(5))-(U_p_pm(5).*tm_p)-(U_b_pm(5).*tm_b(5)))./(tm_b(5)+tm_p);
%     deltaPV(5)=(PV_mix_pm(5).*(tm_p+tm_b(5))-PV_1_pm(5).*tm_p-PV_2_pm(5).*tm_b(5));
%     deltaH(5)=deltaU(5)+deltaPV(5);
%     chi(5)=deltaH(5)./(phi1(5).*phi2(5)*temperature*0.001985875);

  
    
    chiall_1(counter)=chi(1);
    chiall_2(counter)=chi(2);
    chiall_3(counter)=chi(3);
 %   chiall_4(counter)=chi(4);
  %  chiall_5(counter)=chi(5); 
    
    %temp=[400;500]
    %chivalue=[temp chi1]'

%    pdenall_mix(counter)=dcool_mix(2);
%     den_p(counter)=dmean_p;
%     den_b(counter)=dmean_b;
%     den_mix(counter)=dmean_mix;
%     U_b(counter)=Umean_b;
%     U_p(counter)=Umean_p;
%     U_mix(counter)=Umean_mix;
    %densites
    %potential energies 
    %for polymer, binder, and mixture
    %
    %xlabel('T (K)');
    %ylabel('\chi');
end

mean_chi(1)=mean(chiall_1);
mean_chi(2)=mean(chiall_2);
mean_chi(3)=mean(chiall_3);
%mean_chi(4)=mean(chiall_4);
%mean_chi(5)=mean(chiall_5);

figure
    plot(timeframe,chiall_1,'k')
    xlabel('1000 femtoseconds');
    ylabel('/chi');

    hold on
    %yline(mean_chi(1), 'k')
    %errorbar(1,chiall_1, 'k')
    
    plot(timeframe,chiall_2, 'r')
    %yline(mean_chi(2), 'r')
   % errorbar(1,chiall_2, 'r')
%   plot(timeframe, mean_chi2, 'r')
    
    plot(timeframe,chiall_3, 'g')
    %yline(mean_chi(3), 'g')
    %errorbar(1,chiall_3, 'g')
%   plot(timeframe, mean_chi3, 'g')
%    plot(timeframe,chiall_4, 'b')
%    plot(timeframe,chiall_5, 'y')
    
   % plot(timeframe,chiall_single, 'b')
    %errorbar(1,mean_chi_PS_APy_pppm, 'b')
%   plot(timeframe, mean_chi_PS_APy_pppm, 'p')    
%     
% 
%         legtext={'Mix 1','Mix 2','Mix 3'}
%         %legtext={'PE','PB','PP','PS','PVC'}
%         leg=legend(legtext)
%         h1=gca;
%         set(gcf, 'units', 'inches', 'pos', [3.1895 3.9368 17.505 10.916])
% 
%         legend boxoff
%         set(gca,'FontSize',16)
%         set(gcf,'color','w');
%         set(gca,'color','None');
%         box on
%         set(leg,'FontSize',14);
%         set(h1,'TickLength',[.02 .1])
%         set(h1,'XMinorTick','on')

        chi_value(1)=mean(mean_chi);
        sd=std(mean_chi);
        SE_all=sd/sqrt(5);
        chi_value(2)=SE_all;
        ts=tinv([0.025 0.975],2);
       
        chi_value(3)=-SE_all*tinv(0.975,2)+mean(mean_chi);
        chi_value(4)=SE_all*tinv(0.975,2)+mean(mean_chi);
        
chi_value

% hold off 
% time_pot=1:5608;
% 
% figure
% plot(time_pot,M1.PotEng)
%     xlabel('1000 femtoseconds');
%     ylabel('Potential Energy mix');
% hold on 
% plot(time_pot,M2.PotEng)
% plot(time_pot,M3.PotEng)
% hold off 
% 
% figure
% time_pol=1:5608;
% plot(time_pol,R1.PotEng)
%     xlabel('1000 femtoseconds');
%     ylabel('Potential Energy polymer')
% hold on 
% 
% plot(time_pot,R2.PotEng)
% plot(time_pot,R3.PotEng)
% hold off 
% 
% figure
% time_bin=1:5608;
% plot(time_bin,B1.PotEng)
%     xlabel('1000 femtoseconds');
%     ylabel('Potential Energy binder')
% hold on 
% plot(time_bin,B2.PotEng)
% plot(time_bin,B3.PotEng)
% hold off
% 
% figure
% plot(time_pot,M1.PotEng)
%     xlabel('1000 femtoseconds');
%     ylabel('Potential Energy mix');
% hold on 
% plot(time_pot,M2.PotEng)
% plot(time_pot,M3.PotEng)
% hold off 


