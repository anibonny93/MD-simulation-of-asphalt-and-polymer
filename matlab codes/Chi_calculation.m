close all
clear all

path_p=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/Asphalt_PVC/Asphalt_PVC10/PVC/PVC_']; %%%Change for polymer

path_b=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/Pure_asphaltene_runs_pppm_9/AP/AP_']; %%%Change for binder

path_mix=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/Asphalt_PVC/Asphalt_PVC10/AP_PVC/AP_PVC_']; %%%Change for mix

rep=256; %input('Enter replication number:')
mpm=10; %input('Enter monomer per molecule for polymer:')
temperature=500; %input('Enter the temperature:')

[mass_p] = ReadMassFromLAMMPSData([path_p '1_replicated.data']) %'_p' denotes polymer
mass_p=mass_p*rep;
dcool_p=[temperature]';
Ucool_p=dcool_p;

[mass_b] = ReadMassFromLAMMPSData([path_b '1_replicated.data']) %'_b' denotes binder
mass_b=mass_b*rep;
dcool_b=[temperature]';
Ucool_b=dcool_b;
[mass_mix] = ReadMassFromLAMMPSData([path_mix '1_replicated.data']) %'_m' denotes mix
mass_mix=mass_mix*rep;
dcool_mix=[temperature];
Ucool_mix=dcool_mix;

timeframe=input('put the timeframe:')

N=3
for i=1:N;
    filename_p=[path_p num2str(i) '.log'];
    eval(['[R' num2str(i) '] = readlammpslog(filename_p)'])
end

N=3
for i=1:N;
    filename_b=[path_b num2str(i) '.log'];
    eval(['[B' num2str(i) '] = readlammpslog(filename_b)'])
end

N=3
for i=1:N;
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
    phi1=(PV_1_pm)./(PV_1_pm+PV_2_pm)
    phi2=1-phi1
    deltaU=(Umix_pm.*(tm_p+tm_b)-(U_p_pm.*tm_p)-(U_b_pm.*tm_b))./(tm_b+tm_p)
    deltaPV=(PV_mix_pm.*(tm_p+tm_b)-PV_1_pm.*tm_p-PV_2_pm.*tm_b)
    deltaH=deltaU+deltaPV
    chi_1=deltaH./(phi1.*phi2*temperature*0.001985875)
    %temp=[400;500]
    %chivalue=[temp chi1]'
    chiall_1(counter)=chi_1;
    pdenall_mix(counter)=dall_mix(7);
    den_p(counter)=dmean_p;
    den_b(counter)=dmean_b;
    den_mix(counter)=dmean_mix;
    U_b(counter)=Umean_b;
    U_p(counter)=Umean_p;
    U_mix(counter)=Umean_mix;
    mean_chi=mean(chiall_1);
    %densites
    %potential energies 
    %for polymer, binder, and mixture
    %
    %xlabel('T (K)');
    %ylabel('\chi');
end

plot(timeframe,chiall_1);
hold on;
plot(timeframe,ones(length(timeframe),1)*mean_chi, 'k')
xlabel('1000 femtoseconds');
ylabel('\chi');
export_fig 'Chi_AP_PE30' -png -r800 -a1
hold off;

plot(timeframe,den_p);
xlabel('1000 femtoseconds');
ylabel('Density of polymer');
export_fig 'Density_PE30' -png -r800 -a1

plot(timeframe,den_b);
xlabel('1000 femtoseconds')
ylabel('Density of binder')
export_fig 'Density_AP' -png -r800 -a1

plot(timeframe,den_mix);
xlabel('1000 femtoseconds')
ylabel('Density mixture')
export_fig 'Density_AP_PE30' -png -r800 -a1

plot(timeframe,U_p);
xlabel('1000 femtoseconds');
ylabel('PotEng polymer');
export_fig 'PotEng_PE30' -png -r800 -a1

plot(timeframe,U_b);
xlabel('1000 femtoseconds')
ylabel('PotEng binder')
export_fig 'PotEng_AP' -png -r800 -a1

plot(timeframe,U_mix);
xlabel('1000 femtoseconds')
ylabel('PotEng mixture')
export_fig 'PotEng_AP_PE30' -png -r800 -a1