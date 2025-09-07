function [Bonddis,Anglek,Dihedralk]=get_atomistic_data_asphaltCG(prefix,flags)
% prefix='AP';
path=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Oxidized/' prefix '/'];
%CGpath=[path 'CG_1/']
run([path prefix '_config.m'])
bonds_flag=flags(1);
angles_flag=flags(2);
dihedrals_flag=flags(3);
Buildflag=flags(4);
ExtractDataFlag=flags(5); % 1 = read data, 2 = load matlab file
SaveFlag=flags(6);
COMflag=flags(7); %1=center of mass; 2=center of geometry
%% Extract Data
if ExtractDataFlag==1
    for i=simnums
        path=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Oxidized/' prefix '/'];
        datafile=[prefix '_' num2str(i) '_replicated.data'];
        dumpfile=[prefix '_' num2str(i) '.dump'];

        eval(['[system' num2str(i) ']=Read_LAMMPS_Data([path datafile],[path dumpfile])']);
    end
    for j=1:5
        if eval(['~exist(''system' num2str(j) ''',''var'')'])
            eval(['system' num2str(j) '=[]']);
        end
    end
    save(['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Oxidized/' prefix '/' prefix '_systemdata'],'system1','system2','system3','system4','system5');
    clear dis system system1 system2 system3 system4 system5 h dis1 dis 2 dis 3 dis4 dis 5
    close all
    load(['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Oxidized/' prefix '/' prefix '_systemdata']);
elseif ExtractDataFlag==2
    clear dis system system1 system2 system3 system4 system5 h dis1 dis 2 dis 3 dis4 dis 5
    close all
    load(['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Oxidized/' prefix '/' prefix '_systemdata']);
end
%% Build
clear systemwrite
if Buildflag==1
    system=system1;
    atomcounter=1;
    k=1; %k is molecule number
    for k=1:256
        for i=1:size(Beadlist,2)
            systemwrite.attype(atomcounter)=Beadtype{i};
            atoms1=Beadlist{i};
            atoms1_temp=atoms1+(length(system.mass)/256)*(k-1);
            [com1,~]=COM(system.pos(atoms1_temp,:,end),system.lammpsdim1D(end,:),system.mass(atoms1));
            systemwrite.mass(atomcounter)=sum(system.mass(atoms1));
            if systemwrite.attype(atomcounter)=='A' || systemwrite.attype(atomcounter)=='B'
                systemwrite.mass(atomcounter)=sum(system.mass(atoms1))/2;
            end
            systemwrite.pos(atomcounter,:)=(com1-system.lammpsdim(1,:,end));
            for b=find(Bonds(:,1)==i)'
                if ~isfield(systemwrite,'bonds')
                    systemwrite.bonds=[atomcounter Bonds(b,2)-i+atomcounter];
                    systemwrite.bondnames(1,:)=[num2str(i,'%02.0f') '-' num2str(Bonds(b,2),'%02.0f')];
                else
                    systemwrite.bonds(end+1,:)=[atomcounter Bonds(b,2)-i+atomcounter];
                    systemwrite.bondnames(end+1,:)=[num2str(i,'%02.0f') '-' num2str(Bonds(b,2),'%02.0f')];
                end
            end
            for b=find(Angles(:,1)==i)'
                if ~isfield(systemwrite,'angles')
                    systemwrite.angles=[atomcounter Angles(b,2)-i+atomcounter Angles(b,3)-i+atomcounter];
                    systemwrite.anglenames(1,:)=[num2str(i,'%02.0f') '-' num2str(Angles(b,2),'%02.0f') '-' num2str(Angles(b,3),'%02.0f')];
                else
                    systemwrite.angles(end+1,:)=[atomcounter Angles(b,2)-i+atomcounter Angles(b,3)-i+atomcounter];
                    systemwrite.anglenames(end+1,:)=[num2str(i,'%02.0f') '-' num2str(Angles(b,2),'%02.0f') '-' num2str(Angles(b,3),'%02.0f')];
                end
            end
            for b=find(Dihedrals(:,1)==i)'
                if ~isfield(systemwrite,'impropers')
                    systemwrite.impropers=[atomcounter Dihedrals(b,2)-i+atomcounter Dihedrals(b,3)-i+atomcounter Dihedrals(b,4)-i+atomcounter];
                    systemwrite.impropernames(1,:)=[num2str(i,'%02.0f') '-' num2str(Dihedrals(b,2),'%02.0f') '-' num2str(Dihedrals(b,3),'%02.0f') '-' num2str(Dihedrals(b,4),'%02.0f')];
                else
                    systemwrite.impropers(end+1,:)=[atomcounter Dihedrals(b,2)-i+atomcounter Dihedrals(b,3)-i+atomcounter Dihedrals(b,4)-i+atomcounter];
                    systemwrite.impropernames(end+1,:)=[num2str(i,'%02.0f') '-' num2str(Dihedrals(b,2),'%02.0f') '-' num2str(Dihedrals(b,3),'%02.0f') '-' num2str(Dihedrals(b,4),'%02.0f')];
                end
            end
            atomcounter=atomcounter+1;
        end




    end


    systemwrite.bonds=systemwrite.bonds-1;
    systemwrite.bondnames=systemwrite.bondnames;
    systemwrite.angles=systemwrite.angles-1;
    systemwrite.impropers=systemwrite.impropers-1;
    systemwrite.dim=system.lammpsdim1D(end,:);
    systemwrite.natoms=length(systemwrite.attype);
    systemwrite.pos=systemwrite.pos/10;
    systemwrite.dim=systemwrite.dim/10;

    systemwrite.pos=systemwrite.pos+(systemwrite.pos<0).*repmat(systemwrite.dim,systemwrite.natoms,1);
    systemwrite.pos=systemwrite.pos-repmat(systemwrite.dim,systemwrite.natoms,1)/2;

    Write_GSD(systemwrite,[CGpath 'CG_start.gsd'])


    %% Write py file
    BondWrite=xlsread(['/mnt/Shared_Data/PlasticRoads/CoarseGraining/' prefix '/Interactions.xlsx'],'bonds','B6:H300');
    AngleWrite=xlsread(['/mnt/Shared_Data/PlasticRoads/CoarseGraining/' prefix '/Interactions.xlsx'],'angles','B6:I300');
    DihedralWrite=xlsread(['/mnt/Shared_Data/PlasticRoads/CoarseGraining/' prefix '/Interactions.xlsx'],'dihedrals','C7:J301');
    Params=[BondWrite(:,[7 6]);AngleWrite(:,[8 7]);DihedralWrite(:,[7 8])];
    csvwrite([CGpath 'params.csv'],Params);


    WriteOpeningAsphaltCG(prefix,Params)



end
%% Bonds
if bonds_flag==1

    % bondnum=17;
    % atoms1=[31	32	33	34	38	40
    %     ];
    % atoms2=[22	23	38	39	40	41
    %     ];


    for b=1:size(Bonds,1) % loops bonds
        atoms1=Beadlist{Bonds(b,1)};
        atoms2=Beadlist{Bonds(b,2)};
        bondnum=b;
        close all
        clear dis dis1 dis2 dis3 dis4 dis5
        for j=simnums %loop simulations
            system=eval(['system' num2str(j)]);
            %     system.mass=repmat(system.mass,256,1);
            clear dis;
            counter=1;
            for t=tstart:tend %loop timesteps
                for k=1:256 % loop molecules
                    atoms1_temp=atoms1+(length(system.mass)/256)*(k-1);
                    atoms2_temp=atoms2+(length(system.mass)/256)*(k-1);
                    %         atoms1_all=find(ismember(system.attypenum,atoms1));
                    %         atoms2_all=find(ismember(system.attypenum,atoms2));
                    if COMflag==1
                        [com1,~]=COM(system.pos(atoms1_temp,:,t),system.lammpsdim1D(t,:),system.mass(atoms1));
                        [com2,~]=COM(system.pos(atoms2_temp,:,t),system.lammpsdim1D(t,:),system.mass(atoms2));
                    elseif COMflag==2
                        [com1,~]=COM(system.pos(atoms1_temp,:,t),system.lammpsdim1D(t,:));
                        [com2,~]=COM(system.pos(atoms2_temp,:,t),system.lammpsdim1D(t,:));
                    else
                        error('com flag set wrong')
                    end
                    d=abs(com1-com2);

                    for i=1:3
                        if d(i)>system.lammpsdim1D(t,i)/2;
                            d(i)=system.lammpsdim1D(t,i)-d(i);
                        end
                    end
                    dis(counter)=sqrt(sum(d.^2));
                    counter=counter+1;
                end
            end
            %         eval(['dis' num2str(j) '=dis;']);
        end
        %     dis=[dis1];
        h=histogram(dis,'NumBins',80)




        set(gca,'YTickLabel',[]);
        h1=gca;
        % set(gcf, 'units', 'inches', 'pos', [3.1895 3.9368 17.505/2 10.916/2])
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        % set(h1,'TickLength',[.02 .1])
        % set(h1,'XMinorTick','on')
        if SaveFlag==1
            savefig([path 'Bond_' num2str(bondnum) '_hist'])

            eval(['export_fig ' path 'Bond_' num2str(bondnum) '_hist -png -r800 -a1'])
        end

        Bonddis(b)=mean(dis)
    end
end
%% Angles
if angles_flag==1

    for b=1:size(Angles,1) % loops bonds
        atoms1=Beadlist{Angles(b,1)};
        atoms2=Beadlist{Angles(b,2)};
        atoms3=Beadlist{Angles(b,3)};
        anglenum=b;
        close all
        clear dis dis1 dis2 dis3 dis4 dis5
        for j=simnums %loop simulations
            system=eval(['system' num2str(j)]);
            %     system.mass=repmat(system.mass,256,1);
            clear dis;
            counter=1;
            for t=tstart:tend %loop timesteps
                for k=1:256 % loop molecules
                    atoms1_temp=atoms1+(length(system.mass)/256)*(k-1);
                    atoms2_temp=atoms2+(length(system.mass)/256)*(k-1);
                    atoms3_temp=atoms3+(length(system.mass)/256)*(k-1);
                    %         atoms1_all=find(ismember(system.attypenum,atoms1));
                    %         atoms2_all=find(ismember(system.attypenum,atoms2));
                    if COMflag==1
                        [com1,~]=COM(system.pos(atoms1_temp,:,t),system.lammpsdim1D(t,:),system.mass(atoms1));
                        [com2,~]=COM(system.pos(atoms2_temp,:,t),system.lammpsdim1D(t,:),system.mass(atoms2));
                        [com3,~]=COM(system.pos(atoms3_temp,:,t),system.lammpsdim1D(t,:),system.mass(atoms3));
                    elseif COMflag==2
                        [com1,~]=COM(system.pos(atoms1_temp,:,t),system.lammpsdim1D(t,:));
                        [com2,~]=COM(system.pos(atoms2_temp,:,t),system.lammpsdim1D(t,:));
                        [com3,~]=COM(system.pos(atoms3_temp,:,t),system.lammpsdim1D(t,:));
                    else
                        error('com flag set wrong')
                    end
                    for i=1:3
                        if com2(i)-com1(i)>system.lammpsdim1D(t,i)/2;
                            com2(i)=com2(i)-system.lammpsdim1D(t,i);
                        elseif com2(i)-com1(i)<-system.lammpsdim1D(t,i)/2;
                            com2(i)=com2(i)+system.lammpsdim1D(t,i);
                        end
                    end
                    for i=1:3
                        if com3(i)-com2(i)>system.lammpsdim1D(t,i)/2;
                            com3(i)=com3(i)-system.lammpsdim1D(t,i);
                        elseif com3(i)-com2(i)<-system.lammpsdim1D(t,i)/2;
                            com3(i)=com3(i)+system.lammpsdim1D(t,i);
                        end
                    end
                    P0 = com2;
                    P1 = com1;
                    P2 = com3;
                    n1 = (P2 - P0) / norm(P2 - P0);  % Normalized vectors
                    n2 = (P1 - P0) / norm(P1 - P0);
                    angle1 = acosd(dot(n1, n2));
                    angle(counter)=angle1;
                    counter=counter+1;
                end
            end
            %         eval(['dis' num2str(j) '=dis;']);
        end
        %     dis=[dis1];
        h=histogram(angle,'NumBins',80)




        set(gca,'YTickLabel',[]);
        h1=gca;
        % set(gcf, 'units', 'inches', 'pos', [3.1895 3.9368 17.505/2 10.916/2])
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        % set(h1,'TickLength',[.02 .1])
        % set(h1,'XMinorTick','on')
        if SaveFlag==1
            savefig([path 'Angle_' num2str(anglenum) '_hist'])

            eval(['export_fig ' path 'Angle_' num2str(anglenum) '_hist -png -r800 -a1'])
        end
        Anglek(b)=mean(angle)
    end
end

%% Improper Dihedrals
if dihedrals_flag==1

    for b=1:size(Dihedrals,1) % loops bonds
        atoms1=Beadlist{Dihedrals(b,1)};
        atoms2=Beadlist{Dihedrals(b,2)};
        atoms3=Beadlist{Dihedrals(b,3)};
        atoms4=Beadlist{Dihedrals(b,4)};
        dihedralnum=b;
        close all
        clear dis dis1 dis2 dis3 dis4 dis5
        for j=simnums %loop simulations
            system=eval(['system' num2str(j)]);
            %     system.mass=repmat(system.mass,256,1);
            clear dis;
            counter=1;
            for t=tstart:tend %loop timesteps
                for k=1:256 % loop molecules
                    atoms1_temp=atoms1+(length(system.mass)/256)*(k-1);
                    atoms2_temp=atoms2+(length(system.mass)/256)*(k-1);
                    atoms3_temp=atoms3+(length(system.mass)/256)*(k-1);
                    atoms4_temp=atoms4+(length(system.mass)/256)*(k-1);
                    %         atoms1_all=find(ismember(system.attypenum,atoms1));
                    %         atoms2_all=find(ismember(system.attypenum,atoms2));
                    if COMflag==1
                        [com1,~]=COM(system.pos(atoms1_temp,:,t),system.lammpsdim1D(t,:),system.mass(atoms1));
                        [com2,~]=COM(system.pos(atoms2_temp,:,t),system.lammpsdim1D(t,:),system.mass(atoms2));
                        [com3,~]=COM(system.pos(atoms3_temp,:,t),system.lammpsdim1D(t,:),system.mass(atoms3));
                        [com4,~]=COM(system.pos(atoms4_temp,:,t),system.lammpsdim1D(t,:),system.mass(atoms4));
                    elseif COMflag==2
                        [com1,~]=COM(system.pos(atoms1_temp,:,t),system.lammpsdim1D(t,:));
                        [com2,~]=COM(system.pos(atoms2_temp,:,t),system.lammpsdim1D(t,:));
                        [com3,~]=COM(system.pos(atoms3_temp,:,t),system.lammpsdim1D(t,:));
                        [com4,~]=COM(system.pos(atoms4_temp,:,t),system.lammpsdim1D(t,:));
                    else
                        error('com flag set wrong')
                    end
                    for i=1:3
                        if com2(i)-com1(i)>system.lammpsdim1D(t,i)/2;
                            com2(i)=com2(i)-system.lammpsdim1D(t,i);
                        elseif com2(i)-com1(i)<-system.lammpsdim1D(t,i)/2;
                            com2(i)=com2(i)+system.lammpsdim1D(t,i);
                        end
                    end
                    for i=1:3
                        if com3(i)-com1(i)>system.lammpsdim1D(t,i)/2;
                            com3(i)=com1(i)-system.lammpsdim1D(t,i);
                        elseif com3(i)-com1(i)<-system.lammpsdim1D(t,i)/2;
                            com3(i)=com3(i)+system.lammpsdim1D(t,i);
                        end
                    end
                    for i=1:3
                        if com4(i)-com2(i)>system.lammpsdim1D(t,i)/2;
                            com4(i)=com4(i)-system.lammpsdim1D(t,i);
                        elseif com4(i)-com2(i)<-system.lammpsdim1D(t,i)/2;
                            com4(i)=com4(i)+system.lammpsdim1D(t,i);
                        end
                    end
                    for i=1:3
                        if com4(i)-com3(i)>system.lammpsdim1D(t,i)/2;
                            com4(i)=com4(i)-system.lammpsdim1D(t,i);
                        elseif com4(i)-com3(i)<-system.lammpsdim1D(t,i)/2;
                            com4(i)=com4(i)+system.lammpsdim1D(t,i);
                        end
                    end
                    Q=com1;
                    R=com2;
                    S=com3;
                    T=com4;
                    %Plane 1 is QRS; Plane 2 is RST
                    QR=Q-R;
                    QS=Q-S;
                    TR=T-R;
                    TS=T-S;

                    P1=cross(QR*2,QS*2);
                    P1=P1/norm(P1);

                    P2=cross(TR*2,TS*2);
                    P2=P2/norm(P2);

                    dihedral1=acosd((dot(P1,P2)));
                    dihedral(counter)=dihedral1;
                    counter=counter+1;
                end
            end
            %         eval(['dis' num2str(j) '=dis;']);
        end
        %     dis=[dis1];
        h=histogram(dihedral,'NumBins',80)




        set(gca,'YTickLabel',[]);
        h1=gca;
        % set(gcf, 'units', 'inches', 'pos', [3.1895 3.9368 17.505/2 10.916/2])
        set(gca,'FontSize',12)
        set(gcf,'color','w');
        set(gca,'color','None');
        box on
        % set(h1,'TickLength',[.02 .1])
        % set(h1,'XMinorTick','on')
        if SaveFlag==1
            savefig([path 'Dihedral_' num2str(dihedralnum) '_hist'])

            eval(['export_fig ' path 'Dihedral_' num2str(dihedralnum) '_hist -png -r800 -a1'])
        end
        Dihedralk(b)=mean(dihedral)
    end
end