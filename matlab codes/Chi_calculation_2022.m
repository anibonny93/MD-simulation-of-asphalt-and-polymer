close all
clear all

path_p=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2023/PET6_2023/PET6_']; %%%Change for polymer

path_b=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2023/Pure_asphalts_2023/AP/AP_']; %%%Change for binder

path_mix=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2023/Asphalt_PET_2023/AP_PET/AP_PET6_2023_']; %%%Change for mix

rep=256;%input('Enter replication number:')
mpm=6;%input('Enter monomer per molecule for polymer:')
temperature=500;%input('Enter the temperature:')
in=3;%input('Enter the run number')

[mass_p] = ReadMassFromLAMMPSData([path_p '3_replicated.data'])
mass_p=mass_p*rep;
dcool_p=[temperature]';
Ucool_p=dcool_p;

[mass_b] = ReadMassFromLAMMPSData([path_b '3_replicated.data']) %'_b' denotes binder
mass_b=mass_b*rep;
dcool_b=[temperature]';
Ucool_b=dcool_b;

[mass_mix] = ReadMassFromLAMMPSData([path_mix '3_replicated.data'])
mass_mix=mass_mix*rep;
dcool_mix=[temperature];
Ucool_mix=dcool_mix;

timeframe=2000:10:5500;%input('put the timeframe:')


for i=in;
    filename_p=[path_p num2str(i) '.log'];
    eval(['[R' num2str(i) '] = readlammpslog(filename_p)'])
    filename_b=[path_b num2str(i) '.log'];
    eval(['[B' num2str(i) '] = readlammpslog(filename_b)'])
    filename_mix=[path_mix num2str(i) '.log'];
    eval(['[M' num2str(i) '] = readlammpslog(filename_mix)'])
end

twidth=100;
counter=0;

for tstart=timeframe; %loops over time
    tend=tstart+twidth;
    counter=counter+1;
    
    %This section is for the polymer
    
        dcool_p(2)=mean(R2.Volume(tstart:tend));
        Ucool_p(2)=mean(R2.PotEng(tstart:tend));
    
        dcool_p
        Ucool_p
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

    dmean_p=dcool_p(2);
    Umean_p=Ucool_p(2);
    mpm_p=[mpm]; %monomer per molecule for polymer
    vmpm = (mass_p/rep)./(mpm_p.*dmean_p); %volume per mole of monomer for T 400K
    tm_p=rep.*mpm_p;%total number of monomers for polymer for T 400K and 500K

    
    %This section is for the binder

        dcool_b(2)=mean(B2.Volume(tstart:tend));
        Ucool_b(2)=mean(B2.PotEng(tstart:tend));

    
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



    dcool_b
    Ucool_b

    dmean_b=dcool_b(2)
    Umean_b=Ucool_b(2)
    mpm_b=(mass_b/rep)./(vmpm.*dmean_b);
    tm_b=rep.*mpm_b;%total number of monomers for binder for T 400K 



     %This section is for the mix
     

        dcool_mix(2)=mean(M2.Volume(tstart:tend));
        Ucool_mix(2)=mean(M2.PotEng(tstart:tend))
  

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


    dcool_mix
    Ucool_mix

    dmean_mix=dcool_mix(2)
    Umean_mix=Ucool_mix(2)


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
    pdenall_mix(counter)=dcool_mix(2);
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

chi_value(1)=500;
chi_value(2)=mean(chiall);
%chi_value(:,end+1)=std(chiall);
SE=chi_value(:,end)/sqrt(counter);
chi_value(:,end+1)=SE;
chi_value(:,end+1)=-SE*tinv(0.975,counter-1)+mean(chiall);
chi_value(:,end+1)=SE*tinv(0.975,counter-1)+mean(chiall);
mean_chi=mean(chiall);

errorbar(1,chi_value(2),SE)

chi_value 

plot(timeframe,chiall)
