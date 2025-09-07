clear all
close all
path_p=['/mnt/Shared_Data/anic/LONI_output/Small_box/Kspace_ewald/Asphalt_PE30/Pure_runs/PE30/PE30_small_']; %%%Change for polymer

path_b=['/mnt/Shared_Data/anic/LONI_output/Small_box/Kspace_ewald/Asphalt_PE30/Pure_runs/AP/AP_small_']; %%%Change for binder

path_mix=['/mnt/Shared_Data/anic/LONI_output/Small_box/Kspace_ewald/Asphalt_PE30/Mix_runs/AP_PE30/AP_PE30_small_']; %%%Change for mix

%rep=input('Enter replication number:')
%mpm=input('Enter monomer per molecule for polymer:')
%temperature=input('Enter the temperature:')

rep=256;
mpm=30;
temperature=500;

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

%timeframe=input('put the timeframe:')
timeframe=600:10:2200;
N=5
for i=1:N;
    filename_p=[path_p num2str(i) '.log'];
    eval(['[R' num2str(i) '] = readlammpslog(filename_p)'])
end

N=5
for i=1:N;
    filename_b=[path_b num2str(i) '.log'];
    eval(['[B' num2str(i) '] = readlammpslog(filename_b)'])
end

N=5
for i=1:N;
    filename_mix=[path_mix num2str(i) '.log'];
    eval(['[M' num2str(i) '] = readlammpslog(filename_mix)'])
end

twidth=100;
counter=0;

for tstart=timeframe; %loops over time
    tend=tstart+twidth;
    counter=counter+1;
    
    %This section is for the polymer
    
        dcool_p(2)=mean(R1.Volume(tstart:tend));
        Ucool_p(2)=mean(R1.PotEng(tstart:tend));
     
        dcool_p(3)=mean(R2.Volume(tstart:tend));
        Ucool_p(3)=mean(R2.PotEng(tstart:tend));
    
        dcool_p(4)=mean(R3.Volume(tstart:tend));
        Ucool_p(4)=mean(R3.PotEng(tstart:tend));
   
        dcool_p(5)=mean(R4.Volume(tstart:tend));
        Ucool_p(5)=mean(R4.PotEng(tstart:tend));
    
        dcool_p(6)=mean(R5.Volume(tstart:tend));
        Ucool_p(6)=mean(R5.PotEng(tstart:tend));

    
    dcool_p(:,2:end)=((dcool_p(:,2:end)/mass_p).^-1)*1.661; %'_p' signifies polymer
    dall_p=flipud(dcool_p);
    dall_p(:,end+1)=mean(dall_p(:,2:end)')';
    dall_p(:,end+1)=std(dall_p(:,2:end-1)')';
    SE=dall_p(:,end)/sqrt(N);
    dall_p(:,end+1)=-SE*tinv(0.975,N-1)+dall_p(:,N+2);
    dall_p(:,end+1)=SE*tinv(0.975,N-1)+dall_p(:,N+2);

    Uall_p=flipud(Ucool_p);
    Uall_p(:,end+1)=mean(Uall_p(:,2:end)')';
    Uall_p(:,end+1)=std(Uall_p(:,2:end-1)')';
    SE=Uall_p(:,end)/sqrt(N);
    Uall_p(:,end+1)=-SE*tinv(0.975,N-1)+Uall_p(:,N+2);
    Uall_p(:,end+1)=SE*tinv(0.975,N-1)+Uall_p(:,N+2);

    dall_p
    Uall_p

    dmean_p=dall_p(:,7);
    Umean_p=Uall_p(:,7);
    mpm_p=[mpm]; %monomer per molecule for polymer
    vmpm = (mass_p/rep)./(mpm_p.*dmean_p); %volume per mole of monomer for T 400K
    tm_p=rep.*mpm_p;%total number of monomers for polymer for T 400K and 500K

    
    %This section is for the binder

        dcool_b(2)=mean(B1.Volume(tstart:tend));
        Ucool_b(2)=mean(B1.PotEng(tstart:tend));

        dcool_b(3)=mean(B2.Volume(tstart:tend));
        Ucool_b(3)=mean(B2.PotEng(tstart:tend));
        
        dcool_b(4)=mean(B3.Volume(tstart:tend));
        Ucool_b(4)=mean(B3.PotEng(tstart:tend));
    
        dcool_b(5)=mean(B4.Volume(tstart:tend));
        Ucool_b(5)=mean(B4.PotEng(tstart:tend));
 
        dcool_b(6)=mean(B5.Volume(tstart:tend));
        Ucool_b(6)=mean(B5.PotEng(tstart:tend));
    
    dcool_b(:,2:end)=((dcool_b(:,2:end)/mass_b).^-1)*1.661;
    dall_b=flipud(dcool_b);
    dall_b(:,end+1)=mean(dall_b(:,2:end)')';
    dall_b(:,end+1)=std(dall_b(:,2:end-1)')';
    SE=dall_b(:,end)/sqrt(N);
    dall_b(:,end+1)=-SE*tinv(0.975,N-1)+dall_b(:,N+2);
    dall_b(:,end+1)=SE*tinv(0.975,N-1)+dall_b(:,N+2);

    Uall_b=flipud(Ucool_b);
    Uall_b(:,end+1)=mean(Uall_b(:,2:end)')';
    Uall_b(:,end+1)=std(Uall_b(:,2:end-1)')';
    SE=Uall_b(:,end)/sqrt(N);
    Uall_b(:,end+1)=-SE*tinv(0.975,N-1)+Uall_b(:,N+2);
    Uall_b(:,end+1)=SE*tinv(0.975,N-1)+Uall_b(:,N+2);



    dall_b
    Uall_b

    dmean_b=dall_b(:,7)
    Umean_b=Uall_b(:,7)
    mpm_b=(mass_b/rep)./(vmpm.*dmean_b);
    tm_b=rep.*mpm_b;%total number of monomers for binder for T 400K 



     %This section is for the mix
     

        dcool_mix(2)=mean(M1.Volume(tstart:tend));
        Ucool_mix(2)=mean(M1.PotEng(tstart:tend));

        dcool_mix(3)=mean(M2.Volume(tstart:tend));
        Ucool_mix(3)=mean(M2.PotEng(tstart:tend));
    
        dcool_mix(4)=mean(M3.Volume(tstart:tend));
        Ucool_mix(4)=mean(M3.PotEng(tstart:tend));
    
        dcool_mix(5)=mean(M4.Volume(tstart:tend));
        Ucool_mix(5)=mean(M4.PotEng(tstart:tend));
    
        dcool_mix(6)=mean(M5.Volume(tstart:tend));
        Ucool_mix(6)=mean(M5.PotEng(tstart:tend));
  

    dcool_mix(:,2:end)=((dcool_mix(:,2:end)/mass_mix).^-1)*1.661;
    dall_mix=flipud(dcool_mix);
    dall_mix(:,end+1)=mean(dall_mix(:,2:end)')';
    dall_mix(:,end+1)=std(dall_mix(:,2:end-1)')';
    SE=dall_mix(:,end)/sqrt(N);
    dall_mix(:,end+1)=-SE*tinv(0.975,N-1)+dall_mix(:,N+2);
    dall_mix(:,end+1)=SE*tinv(0.975,N-1)+dall_mix(:,N+2);

    Uall_mix=flipud(Ucool_mix);
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
    mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b)./mpm_b.^2)
    mass_1=((mass_p/rep)*tm_p)./mpm_p.^2
    mass_2=((mass_b/rep)*tm_b)./mpm_b.^2
    PV_m=(1./dmean_mix)*(1/1.661).*mass_mix_m
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
    chi=deltaH./(phi1.*phi2*temperature*0.001985875)
    %temp=[400;500]
    %chivalue=[temp chi1]'
    chiall(counter)=chi;
    pdenall_mix(counter)=dall_mix(7);
    den_p(counter)=dmean_p;
    den_b(counter)=dmean_b;
    den_mix(counter)=dmean_mix;
    U_b(counter)=Umean_b;
    U_p(counter)=Umean_p;
    U_mix(counter)=Umean_mix;
    %densites
    %potential energies 
    %for polymer, binder, and mixture
    %
    %xlabel('T (K)');
    %ylabel('\chi');
end

    Umix_pm = Uall_mix(2)./(tm_b+tm_p)
    U_p_pm = Umean_p./tm_p
    U_b_pm = Umean_b./tm_b
    mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b)./mpm_b.^2)
    mass_1=((mass_p/rep)*tm_p)./mpm_p.^2
    mass_2=((mass_b/rep)*tm_b)./mpm_b.^2
    PV_m=(1./dcool_mix(2))*(1/1.661).*mass_mix_m
    PV_1=(1./dcool_p(2))*(1/1.661).*mass_1
    PV_2=(1./dcool_b(2))*(1/1.661).*mass_2
    PV_mix_pm=(0.00001458*PV_m)./(tm_p+tm_b)
    PV_1_pm=(0.00001458*PV_1)./(tm_p)
    PV_2_pm=(0.00001458*PV_2)./(tm_b)
    phi1=(PV_1_pm.*tm_p)./(PV_1_pm.*tm_p+PV_2_pm.*tm_b)
    phi2=1-phi1
    deltaU=(Umix_pm.*(tm_p+tm_b)-(U_p_pm.*tm_p)-(U_b_pm.*tm_b))./(tm_b+tm_p)
    deltaPV=(PV_mix_pm.*(tm_p+tm_b)-PV_1_pm.*tm_p-PV_2_pm.*tm_b)
    deltaH=deltaU+deltaPV
    chi1=deltaH./(phi1.*phi2*temperature*0.001985875)

    Umix_pm = Uall_mix(3)./(tm_b+tm_p)
    U_p_pm = Umean_p./tm_p
    U_b_pm = Umean_b./tm_b
    mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b)./mpm_b.^2)
    mass_1=((mass_p/rep)*tm_p)./mpm_p.^2
    mass_2=((mass_b/rep)*tm_b)./mpm_b.^2
    PV_m=(1./dcool_mix(3))*(1/1.661).*mass_mix_m
    PV_1=(1./dcool_p(3))*(1/1.661).*mass_1
    PV_2=(1./dcool_b(3))*(1/1.661).*mass_2
    PV_mix_pm=(0.00001458*PV_m)./(tm_p+tm_b)
    PV_1_pm=(0.00001458*PV_1)./(tm_p)
    PV_2_pm=(0.00001458*PV_2)./(tm_b)
    phi1=(PV_1_pm.*tm_p)./(PV_1_pm.*tm_p+PV_2_pm.*tm_b)
    phi2=1-phi1
    deltaU=(Umix_pm.*(tm_p+tm_b)-(U_p_pm.*tm_p)-(U_b_pm.*tm_b))./(tm_b+tm_p)
    deltaPV=(PV_mix_pm.*(tm_p+tm_b)-PV_1_pm.*tm_p-PV_2_pm.*tm_b)
    deltaH=deltaU+deltaPV
    chi(2)=deltaH./(phi1.*phi2*temperature*0.001985875)
    
    Umix_pm = Uall_mix(4)./(tm_b+tm_p)
    U_p_pm = Umean_p./tm_p
    U_b_pm = Umean_b./tm_b
    mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b)./mpm_b.^2)
    mass_1=((mass_p/rep)*tm_p)./mpm_p.^2
    mass_2=((mass_b/rep)*tm_b)./mpm_b.^2
    PV_m=(1./dcool_mix(4))*(1/1.661).*mass_mix_m
    PV_1=(1./dcool_p(4))*(1/1.661).*mass_1
    PV_2=(1./dcool_b(4))*(1/1.661).*mass_2
    PV_mix_pm=(0.00001458*PV_m)./(tm_p+tm_b)
    PV_1_pm=(0.00001458*PV_1)./(tm_p)
    PV_2_pm=(0.00001458*PV_2)./(tm_b)
    phi1=(PV_1_pm.*tm_p)./(PV_1_pm.*tm_p+PV_2_pm.*tm_b)
    phi2=1-phi1
    deltaU=(Umix_pm.*(tm_p+tm_b)-(U_p_pm.*tm_p)-(U_b_pm.*tm_b))./(tm_b+tm_p)
    deltaPV=(PV_mix_pm.*(tm_p+tm_b)-PV_1_pm.*tm_p-PV_2_pm.*tm_b)
    deltaH=deltaU+deltaPV
    chi(3)=deltaH./(phi1.*phi2*temperature*0.001985875)
    
    Umix_pm = Uall_mix(5)./(tm_b+tm_p)
    U_p_pm = Umean_p./tm_p
    U_b_pm = Umean_b./tm_b
    mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b)./mpm_b.^2)
    mass_1=((mass_p/rep)*tm_p)./mpm_p.^2
    mass_2=((mass_b/rep)*tm_b)./mpm_b.^2
    PV_m=(1./dcool_mix(5))*(1/1.661).*mass_mix_m
    PV_1=(1./dcool_p(5))*(1/1.661).*mass_1
    PV_2=(1./dcool_b(5))*(1/1.661).*mass_2
    PV_mix_pm=(0.00001458*PV_m)./(tm_p+tm_b)
    PV_1_pm=(0.00001458*PV_1)./(tm_p)
    PV_2_pm=(0.00001458*PV_2)./(tm_b)
    phi1=(PV_1_pm.*tm_p)./(PV_1_pm.*tm_p+PV_2_pm.*tm_b)
    phi2=1-phi1
    deltaU=(Umix_pm.*(tm_p+tm_b)-(U_p_pm.*tm_p)-(U_b_pm.*tm_b))./(tm_b+tm_p)
    deltaPV=(PV_mix_pm.*(tm_p+tm_b)-PV_1_pm.*tm_p-PV_2_pm.*tm_b)
    deltaH=deltaU+deltaPV
    chi(4)=deltaH./(phi1.*phi2*temperature*0.001985875)
    
    Umix_pm = Uall_mix(6)./(tm_b+tm_p)
    U_p_pm = Umean_p./tm_p
    U_b_pm = Umean_b./tm_b
    mass_mix_m=(((mass_p/rep)*tm_p)./mpm_p.^2)+(((mass_b/rep)*tm_b)./mpm_b.^2)
    mass_1=((mass_p/rep)*tm_p)./mpm_p.^2
    mass_2=((mass_b/rep)*tm_b)./mpm_b.^2
    PV_m=(1./dcool_mix(6))*(1/1.661).*mass_mix_m
    PV_1=(1./dcool_p(6))*(1/1.661).*mass_1
    PV_2=(1./dcool_b(6))*(1/1.661).*mass_2
    PV_mix_pm=(0.00001458*PV_m)./(tm_p+tm_b)
    PV_1_pm=(0.00001458*PV_1)./(tm_p)
    PV_2_pm=(0.00001458*PV_2)./(tm_b)
    phi1=(PV_1_pm.*tm_p)./(PV_1_pm.*tm_p+PV_2_pm.*tm_b)
    phi2=1-phi1
    deltaU=(Umix_pm.*(tm_p+tm_b)-(U_p_pm.*tm_p)-(U_b_pm.*tm_b))./(tm_b+tm_p)
    deltaPV=(PV_mix_pm.*(tm_p+tm_b)-PV_1_pm.*tm_p-PV_2_pm.*tm_b)
    deltaH=deltaU+deltaPV
    chi(5)=deltaH./(phi1.*phi2*temperature*0.001985875)
    
    chi(:,end+1)=mean(chi(:,2:end)')';
    chi(:,end+1)=std(chi(:,2:end-1)')';
    SE=chi(:,end)/sqrt(N);eh
    chi(:,end+1)=-SE*tinv(0.975,N-1)+chi(:,N+1);
    chi(:,end+1)=SE*tinv(0.975,N-1)+chi(:,N+1);
    
    %chi
    
    %chiplot=[chi(1) chi(7) chi(8) chi(9); chi(2) chi(7) chi(8) chi(9); chi(3) chi(7) chi(8) chi(9); chi(4) chi(7) chi(8) chi(9); chi(5) chi(7) chi(8) chi(9)]
    chiplot=[chi(6) chi(7) chi(8) chi(9)]
    
