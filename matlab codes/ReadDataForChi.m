clear all
close all
path=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/2023/Asphalt_PE30_2nd_version/BZB_PE30/']; %%%Change
prefix=['BZB_PE30_'];%%%Change
% mass=23530.34*(8^3)/27; %rAAC-1
% mass=27086.022*(8^3)/27; %PMMA
% mass=66786.304*(8^3)/64; %PS
% mass=66786.304*(8^3)/64+23530.34*(8^3)/27; %PS AAC1
% mass=66786.304*(8^3)/64+27086.022*(8^3)/27; %PS PMMA
[mass] = ReadMassFromLAMMPSData([path prefix '1_replicated.data'])
%mass=(3*310.60999+843.636) %use for mix with components
mass=mass*256;
for i=1:5;
    filename=[path prefix num2str(i) '.log'];
    [R] = readlammpslog(filename)
    dcool=[500]';
% dheat=[300:50:700]';
    Ucool=dcool;
% Uheat=dheat;
    tstart=800;
    while tstart<=5500;
         ite_start=20:20:220;
         ite=100:20:300;
         i=1:11;
         a(i)=ite(i)
    for j=1:size(dcool,1);
        tend=tstart+ite(i);
        dcool(j,i+1)=mean(R.Volume(tstart:tend));
        Ucool(j,i+1)=mean(R.PotEng(tstart:tend));
%         [tstart R.Temp(tstart)]
%         [tend R.Temp(tend)]
        %tstart=tstart+600;%use for pure files
        tstart=tstart+ite_start(i+1);%use for mix files
    end 
%     tstart=9700;
%     for j=1:length(dheat)
%         tend=tstart+1000;
%         dheat(j,i+1)=mean(R.Volume(tstart:tend));
%         Uheat(j,i+1)=mean(R.PotEng(tstart:tend));
%         tstart=tstart+2000;
%     end
  end
    %tstart=716+250; %use for mix 
    %tstart=716+250-408; % use for pure
    %tstart=700; %for new molecules mix
    %tstart=730;%for new molecules pure
    %tstart=730;
   
end
dcool(:,2:end)=((dcool(:,2:end)/mass).^-1)*1.661;
% dheat(:,2:end)=((dheat(:,2:end)/mass).^-1)*1.661;
dall=flipud(dcool);
% dall(1:end-1,2:end)=dall(1:end-1,2:end)+flipud(dcool(:,2:end));
% dall(1:end-1,2:end)=dall(1:end-1,2:end)/2;
dall(:,end+1)=mean(dall(:,2:end)')';
dall(:,end+1)=std(dall(:,2:end-1)')';
SE=dall(:,end)/sqrt(N);
dall(:,end+1)=-SE*tinv(0.975,N-1)+dall(:,N+2);
dall(:,end+1)=SE*tinv(0.975,N-1)+dall(:,N+2);

Uall=flipud(Ucool);
% Uall(1:end-1,2:end)=Uall(1:end-1,2:end)+flipud(Ucool(:,2:end));
% Uall(1:end-1,2:end)=Uall(1:end-1,2:end)/2;
Uall(:,end+1)=mean(Uall(:,2:end)')';
Uall(:,end+1)=std(Uall(:,2:end-1)')';
SE=Uall(:,end)/sqrt(N);
Uall(:,end+1)=-SE*tinv(0.975,N-1)+Uall(:,N+2);
Uall(:,end+1)=SE*tinv(0.975,N-1)+Uall(:,N+2);


dall
Uall

dmean=dall(:,7)
Umean=Uall(:,7)
