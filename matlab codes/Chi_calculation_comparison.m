path_p_2=['/mnt/Shared_Data/anic/LONI_output/Small_box/Kspace_ewald/Asphalt_PE30/Pure_runs/PE30/PE30_small_']; %%%Change for polymer

path_b_2=['/mnt/Shared_Data/anic/LONI_output/Small_box/Kspace_ewald/Asphalt_PE30/Pure_runs/TIR/TIR_small_']; %%%Change for binder

path_mix_2=['/mnt/Shared_Data/anic/LONI_output/Small_box/Kspace_ewald/Asphalt_PE30/Mix_runs/TIR_PE30/TIR_PE30_small_']; %%%Change for mix


rep_2=input('Enter replication number:')
mpm_2=input('Enter monomer per molecule for polymer:')
temperature_2=input('Enter the temperature:')

[mass_p_2] = ReadMassFromLAMMPSData([path_p_2 '1_replicated.data'])
mass_p_2=mass_p_2*rep_2;
dcool_p_2=[temperature_2]';
Ucool_p_2=dcool_p_2;

[mass_b_2] = ReadMassFromLAMMPSData([path_b_2 '1_replicated.data']) %'_b' denotes binder
mass_b_2=mass_b_2*rep_2;
dcool_b_2=[temperature_2]';
Ucool_b_2=dcool_b_2;

[mass_mix_2] = ReadMassFromLAMMPSData([path_mix_2 '1_replicated.data'])
mass_mix_2=mass_mix_2*rep_2;
dcool_mix_2=[temperature_2];
Ucool_mix_2=dcool_mix_2;

timeframe_2=input('put the timeframe:')

N_2=5
for i=1:N_2;
    filename_p_2=[path_p_2 num2str(i) '.log'];
    eval(['[R' num2str(i) '] = readlammpslog(filename_p_2)'])
end

N_2=5
for i=1:N_2;
    filename_b_2=[path_b_2 num2str(i) '.log'];
    eval(['[B' num2str(i) '] = readlammpslog(filename_b_2)'])
end

N_2=5
for i=1:N_2;
    filename_mix_2=[path_mix_2 num2str(i) '.log'];
    eval(['[M' num2str(i) '] = readlammpslog(filename_mix_2)'])
end

twidth_2=100;
counter_2=0;

for tstart_2=timeframe_2; %loops over time
    tend_2=tstart_2+twidth_2;
    counter_2=counter_2+1;
    
    %This section is for the polymer
    
        dcool_p_2(2)=mean(R1.Volume(tstart_2:tend_2));
        Ucool_p_2(2)=mean(R1.PotEng(tstart_2:tend_2));
     
        dcool_p_2(3)=mean(R2.Volume(tstart_2:tend_2));
        Ucool_p_2(3)=mean(R2.PotEng(tstart_2:tend_2));
    
        dcool_p_2(4)=mean(R3.Volume(tstart_2:tend_2));
        Ucool_p_2(4)=mean(R3.PotEng(tstart_2:tend_2));
   
        dcool_p_2(5)=mean(R4.Volume(tstart_2:tend_2));
        Ucool_p_2(5)=mean(R4.PotEng(tstart_2:tend_2));
    
        dcool_p_2(6)=mean(R5.Volume(tstart_2:tend_2));
        Ucool_p_2(6)=mean(R5.PotEng(tstart_2:tend_2));

    
    dcool_p_2(:,2:end)=((dcool_p_2(:,2:end)/mass_p_2).^-1)*1.661; %'_p' signifies polymer
    dall_p_2=flipud(dcool_p_2);
    dall_p_2(:,end+1)=mean(dall_p_2(:,2:end)')';
    dall_p_2(:,end+1)=std(dall_p_2(:,2:end-1)')';
    SE_2=dall_p_2(:,end)/sqrt(N_2);
    dall_p_2(:,end+1)=-SE_2*tinv(0.975,N_2-1)+dall_p_2(:,N_2+2);
    dall_p_2(:,end+1)=SE_2*tinv(0.975,N_2-1)+dall_p_2(:,N_2+2);

    Uall_p_2=flipud(Ucool_p_2);
    Uall_p_2(:,end+1)=mean(Uall_p_2(:,2:end)')';
    Uall_p_2(:,end+1)=std(Uall_p_2(:,2:end-1)')';
    SE_2=Uall_p_2(:,end)/sqrt(N_2);
    Uall_p_2(:,end+1)=-SE_2*tinv(0.975,N_2-1)+Uall_p_2(:,N_2+2);
    Uall_p_2(:,end+1)=SE_2*tinv(0.975,N_2-1)+Uall_p_2(:,N_2+2);

    dall_p_2
    Uall_p_2

    dmean_p_2=dall_p_2(:,7);
    Umean_p_2=Uall_p_2(:,7);
    mpm_p_2=[mpm_2]; %monomer per molecule for polymer
    vmpm_2 = (mass_p_2/rep_2)./(mpm_p_2.*dmean_p_2); %volume per mole of monomer for T 400K
    tm_p_2=rep_2.*mpm_p_2;%total number of monomers for polymer for T 400K and 500K

    
    %This section is for the binder

        dcool_b_2(2)=mean(B1.Volume(tstart_2:tend_2));
        Ucool_b_2(2)=mean(B1.PotEng(tstart_2:tend_2));

        dcool_b_2(3)=mean(B2.Volume(tstart_2:tend_2));
        Ucool_b_2(3)=mean(B2.PotEng(tstart_2:tend_2));
        
        dcool_b_2(4)=mean(B3.Volume(tstart_2:tend_2));
        Ucool_b_2(4)=mean(B3.PotEng(tstart_2:tend_2));
    
        dcool_b_2(5)=mean(B4.Volume(tstart_2:tend_2));
        Ucool_b_2(5)=mean(B4.PotEng(tstart_2:tend_2));
 
        dcool_b_2(6)=mean(B5.Volume(tstart_2:tend_2));
        Ucool_b_2(6)=mean(B5.PotEng(tstart_2:tend_2));
    
    dcool_b_2(:,2:end)=((dcool_b_2(:,2:end)/mass_b_2).^-1)*1.661;
    dall_b_2=flipud(dcool_b_2);
    dall_b_2(:,end+1)=mean(dall_b_2(:,2:end)')';
    dall_b_2(:,end+1)=std(dall_b_2(:,2:end-1)')';
    SE_2=dall_b_2(:,end)/sqrt(N_2);
    dall_b_2(:,end+1)=-SE_2*tinv(0.975,N_2-1)+dall_b_2(:,N_2+2);
    dall_b_2(:,end+1)=SE_2*tinv(0.975,N_2-1)+dall_b_2(:,N_2+2);

    Uall_b_2=flipud(Ucool_b_2);
    Uall_b_2(:,end+1)=mean(Uall_b_2(:,2:end)')';
    Uall_b_2(:,end+1)=std(Uall_b_2(:,2:end-1)')';
    SE_2=Uall_b_2(:,end)/sqrt(N_2);
    Uall_b_2(:,end+1)=-SE_2*tinv(0.975,N_2-1)+Uall_b_2(:,N_2+2);
    Uall_b_2(:,end+1)=SE_2*tinv(0.975,N_2-1)+Uall_b_2(:,N_2+2);



    dall_b_2
    Uall_b_2

    dmean_b_2=dall_b_2(:,7)
    Umean_b_2=Uall_b_2(:,7)
    mpm_b_2=(mass_b_2/rep_2)./(vmpm_2.*dmean_b_2);
    tm_b_2=rep_2.*mpm_b_2;%total number of monomers for binder for T 400K 



     %This section is for the mix
     

        dcool_mix_2(2)=mean(M1.Volume(tstart_2:tend_2));
        Ucool_mix_2(2)=mean(M1.PotEng(tstart_2:tend_2));

        dcool_mix_2(3)=mean(M2.Volume(tstart_2:tend_2));
        Ucool_mix_2(3)=mean(M2.PotEng(tstart_2:tend_2));
    
        dcool_mix_2(4)=mean(M3.Volume(tstart_2:tend_2));
        Ucool_mix_2(4)=mean(M3.PotEng(tstart_2:tend_2));
    
        dcool_mix_2(5)=mean(M4.Volume(tstart_2:tend_2));
        Ucool_mix_2(5)=mean(M4.PotEng(tstart_2:tend_2));
    
        dcool_mix_2(6)=mean(M5.Volume(tstart_2:tend_2));
        Ucool_mix_2(6)=mean(M5.PotEng(tstart_2:tend_2));
  

    dcool_mix_2(:,2:end)=((dcool_mix_2(:,2:end)/mass_mix_2).^-1)*1.661;
    dall_mix_2=flipud(dcool_mix_2);
    dall_mix_2(:,end+1)=mean(dall_mix_2(:,2:end)')';
    dall_mix_2(:,end+1)=std(dall_mix_2(:,2:end-1)')';
    SE_2=dall_mix_2(:,end)/sqrt(N_2);
    dall_mix_2(:,end+1)=-SE_2*tinv(0.975,N_2-1)+dall_mix_2(:,N_2+2);
    dall_mix_2(:,end+1)=SE_2*tinv(0.975,N_2-1)+dall_mix_2(:,N_2+2);

    Uall_mix_2=flipud(Ucool_mix_2);
    Uall_mix_2(:,end+1)=mean(Uall_mix_2(:,2:end)')';
    Uall_mix_2(:,end+1)=std(Uall_mix_2(:,2:end-1)')';
    SE_2=Uall_mix_2(:,end)/sqrt(N_2);
    Uall_mix_2(:,end+1)=-SE_2*tinv(0.975,N_2-1)+Uall_mix_2(:,N_2+2);
    Uall_mix_2(:,end+1)=SE_2*tinv(0.975,N_2-1)+Uall_mix_2(:,N_2+2);
    
    


    dall_mix_2
    Uall_mix_2

    dmean_mix_2=dall_mix_2(:,7)
    Umean_mix_2=Uall_mix_2(:,7)


    Umix_pm_2 = Umean_mix_2./(tm_b_2+tm_p_2)
    U_p_pm_2 = Umean_p_2./tm_p_2
    U_b_pm_2 = Umean_b_2./tm_b_2
    mass_mix_m_2=(((mass_p_2/rep_2)*tm_p_2)./mpm_p_2.^2)+(((mass_b_2/rep_2)*tm_b_2)./mpm_b_2.^2)
    mass_1_2=((mass_p_2/rep_2)*tm_p_2)./mpm_p_2.^2
    mass_2_2=((mass_b_2/rep_2)*tm_b_2)./mpm_b_2.^2
    PV_m_2=(1./dmean_mix_2)*(1/1.661).*mass_mix_m_2
    PV_1_2=(1./dmean_p_2)*(1/1.661).*mass_1_2
    PV_2_2=(1./dmean_b_2)*(1/1.661).*mass_2_2
    PV_mix_pm_2=(0.00001458*PV_m_2)./(tm_p_2+tm_b_2)
    PV_1_pm_2=(0.00001458*PV_1_2)./(tm_p_2)
    PV_2_pm_2=(0.00001458*PV_2_2)./(tm_b_2)
    phi1_2=(PV_1_pm_2.*tm_p_2)./(PV_1_pm_2.*tm_p_2+PV_2_pm_2.*tm_b_2)
    phi2_2=1-phi1_2
    deltaU_2=(Umix_pm_2.*(tm_p_2+tm_b_2)-(U_p_pm_2.*tm_p_2)-(U_b_pm_2.*tm_b_2))./(tm_b_2+tm_p_2)
    deltaPV_2=(PV_mix_pm_2.*(tm_p_2+tm_b_2)-PV_1_pm_2.*tm_p_2-PV_2_pm_2.*tm_b_2)
    deltaH_2=deltaU_2+deltaPV_2
    chi1_2=deltaH_2./(phi1_2.*phi2_2*temperature_2*0.001985875)
    %temp=[400;500]
    %chivalue=[temp chi1]'
    chiall_2(counter_2)=chi1_2;
    pdenall_mix(counter_2)=dall_mix_2(7);
    den_p_2(counter_2)=dmean_p_2;
    den_b_2(counter_2)=dmean_b_2;
    den_mix_2(counter_2)=dmean_mix_2;
    U_b_2(counter_2)=Umean_b_2;
    U_p_2(counter_2)=Umean_p_2;
    U_mix_2(counter_2)=Umean_mix_2;
    %densites
    %potential energies 
    %for polymer, binder, and mixture
    %
    %xlabel('T (K)');
    %ylabel('\chi');
end

%plot(timeframe_2,chiall_2);
%xlabel('1000 femtoseconds');
%ylabel('\chi');
%export_fig 'Chi_BZB_PE30' -png -r800 -a1

%plot(timeframe_2,den_p_2);
%xlabel('1000 femtoseconds');
%ylabel('Density of polymer');
%export_fig 'Density_PE30' -png -r800 -a1

%plot(timeframe_2,den_b_2);
%xlabel('1000 femtoseconds')
%ylabel('Density of binder')
%export_fig 'Density_BZB' -png -r800 -a1

%plot(timeframe_2,den_mix_2);
%xlabel('1000 femtoseconds')
%ylabel('Density mixture')
%export_fig 'Density_BZB_PE30' -png -r800 -a1

%plot(timeframe_2,U_p_2);
%xlabel('1000 femtoseconds');
%ylabel('PotEng polymer');
%export_fig 'PotEng_PE30' -png -r800 -a1

%plot(timeframe_2,U_b_2);
%xlabel('1000 femtoseconds')
%ylabel('PotEng binder')
%export_fig 'PotEng_BZB' -png -r800 -a1

%plot(timeframe_2,U_mix_2);
%xlabel('1000 femtoseconds')
%ylabel('PotEng mixture')
%export_fig 'PotEng_BZB_PE30' -png -r800 -a1

%plot(timeframe_2,chiall,chiall_2);
%xlabel('1000 femtoseconds');
%ylabel('Density of polymer');
%export_fig 'Density_PE30' -png -r800 -a1