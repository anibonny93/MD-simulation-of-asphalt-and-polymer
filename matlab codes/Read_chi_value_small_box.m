clear all
close all
path_p=['/mnt/Shared_Data/anic/LONI_output/Small_box/PS_Benzene/PS/PS10_']; %%%Change for polymer
%prefix_p=['Polypropylene_'];%%%Change for polymer
path_b=['/mnt/Shared_Data/anic/LONI_output/Small_box/PS_Cyclohexane/Cyclohexane/Cyclohexane_']; %%%Change for binder
%prefix_b=['Asphaltene_phenol_poly_'];%%%Change for binder
path_mix=['/mnt/Shared_Data/anic/LONI_output/Small_box/PS_Cyclohexane/Mix/PS_CH_']; %%%Change for mix
%prefix_mix=['AP_PP_'];%%%Change for mix

% mass=23530.34*(8^3)/27; %rAAC-1
% mass=27086.022*(8^3)/27; %PMMA
% mass=66786.304*(8^3)/64; %PS
% mass=66786.304*(8^3)/64+23530.34*(8^3)/27; %PS AAC1
% mass=66786.304*(8^3)/64+27086.022*(8^3)/27; %PS PMMA
[mass_p] = ReadMassFromLAMMPSData([path_p '1_replicated.data'])
%mass=(3*310.60999+843.636) %use for mix with components
mass_p=mass_p*4^3;

dcool_p=[300]';
% dheat=[300:50:700]';
Ucool_p=dcool_p;
% Uheat=dheat;
N=5
for i=1:N;
    filename=[path_p num2str(i) '.log'];
    [R] = readlammpslog(filename)
    %tstart=716+250; %use for mix 
    %tstart=716+250-408; % use for pure
    %tstart=700; %for new molecules mix
    tstart=730;%for new molecules pure
    
    for j=1:size(dcool_p,1);
        tend=tstart+100;
        dcool_p(j,i+1)=mean(R.Volume(tstart:tend));
        Ucool_p(j,i+1)=mean(R.PotEng(tstart:tend));
%         [tstart R.Temp(tstart)]
%         [tend R.Temp(tend)]
        %tstart=tstart+600;%use for pure files
        tstart=tstart+400;%use for mix files
    end
%     tstart=9700;
%     for j=1:length(dheat)
%         tend=tstart+1000;
%         dheat(j,i+1)=mean(R.Volume(tstart:tend));
%         Uheat(j,i+1)=mean(R.PotEng(tstart:tend));
%         tstart=tstart+2000;
%     end
end
dcool_p(:,2:end)=((dcool_p(:,2:end)/mass_p).^-1)*1.661; %'_p' signifies polymer
% dheat(:,2:end)=((dheat(:,2:end)/mass).^-1)*1.661;
dall_p=flipud(dcool_p);
% dall(1:end-1,2:end)=dall(1:end-1,2:end)+flipud(dcool(:,2:end));
% dall(1:end-1,2:end)=dall(1:end-1,2:end)/2;
dall_p(:,end+1)=mean(dall_p(:,2:end)')';
dall_p(:,end+1)=std(dall_p(:,2:end-1)')';
SE=dall_p(:,end)/sqrt(N);
dall_p(:,end+1)=-SE*tinv(0.975,N-1)+dall_p(:,N+2);
dall_p(:,end+1)=SE*tinv(0.975,N-1)+dall_p(:,N+2);

Uall_p=flipud(Ucool_p);
% Uall(1:end-1,2:end)=Uall(1:end-1,2:end)+flipud(Ucool(:,2:end));
% Uall(1:end-1,2:end)=Uall(1:end-1,2:end)/2;
Uall_p(:,end+1)=mean(Uall_p(:,2:end)')';
Uall_p(:,end+1)=std(Uall_p(:,2:end-1)')';
SE=Uall_p(:,end)/sqrt(N);
Uall_p(:,end+1)=-SE*tinv(0.975,N-1)+Uall_p(:,N+2);
Uall_p(:,end+1)=SE*tinv(0.975,N-1)+Uall_p(:,N+2);



dall_p
Uall_p

dmean_p=dall_p(:,7);
Umean_p=Uall_p(:,7);
mpm_p=[10 ; 10]; %monomer per molecule for polymer
vmpm = (mass_p/256)./(mpm_p.*dmean_p); %volume per mole of monomer for T 400K
%vmpm2 = (mass_p/512)/(mpm_p.*dmean_p); %volume per mole of monomer for T 500K
vmpm
%vmpm2
rep=512 %replication of molecules
tm_p=rep.*mpm_p;%total number of monomers for polymer for T 400K and 500K



% mass=23530.34*(8^3)/27; %rAAC-1
% mass=27086.022*(8^3)/27; %PMMA
% mass=66786.304*(8^3)/64; %PS
% mass=66786.304*(8^3)/64+23530.34*(8^3)/27; %PS AAC1
% mass=66786.304*(8^3)/64+27086.022*(8^3)/27; %PS PMMA
[mass_b] = ReadMassFromLAMMPSData([path_b '1_replicated.data']) %'_b' denotes binder
%mass=(3*310.60999+843.636) %use for mix with components
mass_b=mass_b*4^3;

dcool_b=[325]';
% dheat=[300:50:700]';
Ucool_b=dcool_b;
% Uheat=dheat;
N=5
for i=1:N;
    filename=[path_b num2str(i) '.log'];
    [R] = readlammpslog(filename)
    %tstart=716+250; %use for mix 
    %tstart=716+250-408; % use for pure
    %tstart=700; %for new molecules mix
    tstart=730;%for new molecules pure
    
    for j=1:size(dcool_b,1);
        tend=tstart+100;
        dcool_b(j,i+1)=mean(R.Volume(tstart:tend));
        Ucool_b(j,i+1)=mean(R.PotEng(tstart:tend));
%         [tstart R.Temp(tstart)]
%         [tend R.Temp(tend)]
        %tstart=tstart+600;%use for pure files
        tstart=tstart+400;%use for mix files
    end
%     tstart=9700;
%     for j=1:length(dheat)
%         tend=tstart+1000;
%         dheat(j,i+1)=mean(R.Volume(tstart:tend));
%         Uheat(j,i+1)=mean(R.PotEng(tstart:tend));
%         tstart=tstart+2000;
%     end
end
dcool_b(:,2:end)=((dcool_b(:,2:end)/mass_b).^-1)*1.661;
% dheat(:,2:end)=((dheat(:,2:end)/mass).^-1)*1.661;
dall_b=flipud(dcool_b);
% dall(1:end-1,2:end)=dall(1:end-1,2:end)+flipud(dcool(:,2:end));
% dall(1:end-1,2:end)=dall(1:end-1,2:end)/2;
dall_b(:,end+1)=mean(dall_b(:,2:end)')';
dall_b(:,end+1)=std(dall_b(:,2:end-1)')';
SE=dall_b(:,end)/sqrt(N);
dall_b(:,end+1)=-SE*tinv(0.975,N-1)+dall_b(:,N+2);
dall_b(:,end+1)=SE*tinv(0.975,N-1)+dall_b(:,N+2);

Uall_b=flipud(Ucool_b);
% Uall(1:end-1,2:end)=Uall(1:end-1,2:end)+flipud(Ucool(:,2:end));
% Uall(1:end-1,2:end)=Uall(1:end-1,2:end)/2;
Uall_b(:,end+1)=mean(Uall_b(:,2:end)')';
Uall_b(:,end+1)=std(Uall_b(:,2:end-1)')';
SE=Uall_b(:,end)/sqrt(N);
Uall_b(:,end+1)=-SE*tinv(0.975,N-1)+Uall_b(:,N+2);
Uall_b(:,end+1)=SE*tinv(0.975,N-1)+Uall_b(:,N+2);



dall_b
Uall_b

dmean_b=dall_b(:,7)
Umean_b=Uall_b(:,7)
mpm_b=(mass_b/256)./(vmpm.*dmean_b);
%mpm_b_2=(mass_b/512)/(vmpm2*dmean_b(2));
%mpm_b= [mpm_b_1, mpm_b_2];
tm_b=rep.*mpm_b;%total number of monomers for binder for T 400K 
%tm_b_2=rep*mpm_b_2;%total number of monomers for binder for T 500K
%tm_b=[tm_b,tm_b_2];


% mass=23530.34*(8^3)/27; %rAAC-1
% mass=27086.022*(8^3)/27; %PMMA
% mass=66786.304*(8^3)/64; %PS
% mass=66786.304*(8^3)/64+23530.34*(8^3)/27; %PS AAC1
% mass=66786.304*(8^3)/64+27086.022*(8^3)/27; %PS PMMA
[mass_mix] = ReadMassFromLAMMPSData([path_mix '1_replicated.data'])
%mass=(3*310.60999+843.636) %use for mix with components
mass_mix=mass_mix*4^3;

dcool_mix=[325]';
% dheat=[300:50:700]';
Ucool_mix=dcool_mix;
% Uheat=dheat;
N=5
for i=1:N;
    filename=[path_mix num2str(i) '.log'];
    [R] = readlammpslog(filename)
    %tstart=716+250; %use for mix 
    %tstart=716+250-408; % use for pure
    %tstart=700; %for new molecules mix
    tstart=730;%for new molecules pure
    
    for j=1:size(dcool_mix,1);
        tend=tstart+100;
        dcool_mix(j,i+1)=mean(R.Volume(tstart:tend));
        Ucool_mix(j,i+1)=mean(R.PotEng(tstart:tend));
%         [tstart R.Temp(tstart)]
%         [tend R.Temp(tend)]
        %tstart=tstart+600;%use for pure files
        tstart=tstart+400;%use for mix files
    end
%     tstart=9700;
%     for j=1:length(dheat)
%         tend=tstart+1000;
%         dheat(j,i+1)=mean(R.Volume(tstart:tend));
%         Uheat(j,i+1)=mean(R.PotEng(tstart:tend));
%         tstart=tstart+2000;
%     end
end
dcool_mix(:,2:end)=((dcool_mix(:,2:end)/mass_mix).^-1)*1.661;
% dheat(:,2:end)=((dheat(:,2:end)/mass).^-1)*1.661;
dall_mix=flipud(dcool_mix);
% dall(1:end-1,2:end)=dall(1:end-1,2:end)+flipud(dcool(:,2:end));
% dall(1:end-1,2:end)=dall(1:end-1,2:end)/2;
dall_mix(:,end+1)=mean(dall_mix(:,2:end)')';
dall_mix(:,end+1)=std(dall_mix(:,2:end-1)')';
SE=dall_mix(:,end)/sqrt(N);
dall_mix(:,end+1)=-SE*tinv(0.975,N-1)+dall_mix(:,N+2);
dall_mix(:,end+1)=SE*tinv(0.975,N-1)+dall_mix(:,N+2);

Uall_mix=flipud(Ucool_mix);
% Uall(1:end-1,2:end)=Uall(1:end-1,2:end)+flipud(Ucool(:,2:end));
% Uall(1:end-1,2:end)=Uall(1:end-1,2:end)/2;
Uall_mix(:,end+1)=mean(Uall_mix(:,2:end)')';
Uall_mix(:,end+1)=std(Uall_mix(:,2:end-1)')';
SE=Uall_mix(:,end)/sqrt(N);
Uall_mix(:,end+1)=-SE*tinv(0.975,N-1)+Uall_mix(:,N+2);
Uall_mix(:,end+1)=SE*tinv(0.975,N-1)+Uall_mix(:,N+2);



dall_mix
Uall_mix

dmean_mix=dall_mix(:,7)
Umean_mix=Uall_mix(:,7)



Umix_pm = Umean_mix./(tm_b+tm_p)
U_p_pm = Umean_p./tm_p
U_b_pm = Umean_b./tm_b
mass_mix=(((mass_p/256)*tm_p)./mpm_p.^2)+(((mass_b/256)*tm_b)./mpm_b.^2)
mass_1=((mass_p/256)*tm_p)./mpm_p.^2
mass_2=((mass_b/256)*tm_b)./mpm_b.^2
PV_m=(1./dmean_mix)*(1/1.661).*mass_mix
PV_1=(1./dmean_p)*(1/1.661).*mass_1
PV_2=(1./dmean_b)*(1/1.661).*mass_2
PV_mix_pm=(0.00001458*PV_m)./(tm_p+tm_b)
PV_1_pm=(0.00001458*PV_1)./(tm_p)
PV_2_pm=(0.00001458*PV_2)./(tm_b)
phi1=(PV_1_pm.*tm_p)./(PV_1_pm.*tm_p+PV_2_pm.*tm_b)
phi2=1-phi1
deltaU=(Umix_pm.*(tm_p+tm_b)-(U_p_pm.*tm_p)-(U_b_pm.*tm_b))./(tm_b+tm_p)
deltaPV=(PV_mix_pm.*(tm_p+tm_b)-PV_1_pm.*tm_p-PV_2_pm.*tm_b)
deltaH=deltaU+deltaPV
chi=deltaH./(phi1.*phi2*325*0.001985875)
temp=[325]
chivalue=[temp chi]