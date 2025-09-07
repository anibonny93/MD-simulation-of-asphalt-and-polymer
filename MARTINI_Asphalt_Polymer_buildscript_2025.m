TypeMAT=[9   9    189    0   180    0    90    0    90    0    36    0    90    0    90   0   0   0   0   9   9   9   0   0   0   0    0   0   0]
RestrictPolymerFlag=0;
%TypeMAT=[0     0     0     0     1     1     0    0    0    0    0     0     0     0     0    0    0     0      0     0     0     0     0     0     0     0];
%TypeMAT(18)=41%[32    32   104    88    32    32   120    40    32    24    24    16     0     4     0     0     0   0]%8*[4 4 13 11 4 4 15 5 4 3 3 2 0 0 0 0 0]%[8    8    26    22    8    8    30    10    8    6    6     4     0     6     0     0     0]%4*[4 4 13 11 4 4 15 5 4 3 3 2 0 1 0 0 0]; %[0 0 0 0 0 0 0 0 0 24 24 24 900] %10*[4 4 13 11 4 4 15 5 4 3 3 2 0];  
%TypeMAT=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 300 0];
%1 = squalane 2=hopane 3=dochn 4=phpn 5= quinolinohopane 6=pyrodinohopane
%7=BZB 8 = TBO 9=TIR 10=Asphaltene Thiophene (AT) 11=AP 12= APy 13 =
%Heptane/Octane 14 = PE 15 = PP 16 = PB 17 = PS 18 = SBS 19 = SBSr
BoxSize=[100 100 100]*5
DegreeofPolymerization=[400 100 100 200]; %PE PP PB PS 
SBSDegreeOfPolymerization=[10 80 10]; %SBS
SBSRadialDegreeOfPolymerization=[5 100 99 5] %SBS radial 1 = first S 2 =B horizontal 3=B vertical (probably want this to be 1 less since the midpoint is on the first B 4= the other 3 S
clear system


%write out list of ranges for molecules
Molecules={'Squalane','Hopane','DOCHN','DOCHN_oxi','PHPN','PHPN_oxi','Quinolinohopane','QH_oxi','Pyridinohopane','PH_oxi','Benzobisbenzothiophene','BZB_oxi','Trimethylbenzene-oxane','TBO_oxi','Thio-Isorenieratane','TIR_oxi','AP_oxi_d30','APy_oxi_d30', 'AT_oxi_d30','Asphaltene Thiophene','Asphaltene Phenol','Asphaltene Pyrrole','Heptane','PE','PP','PB','PS','SBS','SBSr'};
f_MOL_ID=fopen('MoleculeList.txt','w');
for i=1:29
    if i==1
        LimStart=1;
        SatStart=1;
    else
        LimStart=LimStart+TypeMAT(i-1);
    end
    LimEnd=LimStart+TypeMAT(i);
    if i==2
        SatEnd=LimEnd;
    elseif i==3
        AroStart=LimStart;
    elseif i==6
        AroEnd=LimEnd;
    elseif i==7
        ResStart=LimStart;
    elseif i==16
        ResEnd=LimEnd;
    elseif i==17
        AsphStart=LimStart;
    elseif i==22
        AsphEnd=LimEnd;
    elseif i==23
        PolyStart=LimStart;
    elseif i==29
        PolyEnd=LimEnd;
    end
    str=Molecules{i};
    for j=length(str)+1:29
        str=[str ' '];
    end

    fprintf(f_MOL_ID,[str '\t' num2str(LimStart-1) '\t' num2str(LimEnd) '\n'])
end

% Write SARA
fprintf(f_MOL_ID,['\n\n'])
str='Saturates';
for j=length(str)+1:29
    str=[str ' '];
end
fprintf(f_MOL_ID,[str '\t' num2str(SatStart-1) '\t' num2str(SatEnd) '\n'])
str='Aromatics';
for j=length(str)+1:29
    str=[str ' '];
end
fprintf(f_MOL_ID,[str '\t' num2str(AroStart-1) '\t' num2str(AroEnd) '\n'])
str='Resins';
for j=length(str)+1:29
    str=[str ' '];
end
fprintf(f_MOL_ID,[str '\t' num2str(ResStart-1) '\t' num2str(ResEnd) '\n']);
str='Asphaltenes';
for j=length(str)+1:29
    str=[str ' '];
end
fprintf(f_MOL_ID,[str '\t' num2str(AsphStart-1) '\t' num2str(AsphEnd) '\n']);

str='Polymer';
for j=length(str)+1:29
    str=[str ' '];
end
fprintf(f_MOL_ID,[str '\t' num2str(PolyStart-1) '\t' num2str(PolyEnd) '\n']);



MOLPOS=[];
counter=0;
molcounter=0;

system.natoms=0;
system.dim=BoxSize;
system.lammpsdim=[0 BoxSize(1);0 BoxSize(2);0 BoxSize(3)];
system.attype=[];
system.mass=[];
system.charge=[];
system.pos=[];
system.molid=[];

system.bonds=[];
system.bondnames=[];
system.BondPotential=[];

system.angles=[];
system.anglenames=[];
system.AnglePotential=[];

system.impropers=[];
system.impropernames=[];
system.ImproperPotential=[];

system.dihedrals=[];

TotalMASS=0;

Namelist={'SQ','HP','DCN','DOCHN_oxi','PHPN','PHPN_oxi','QH','QH_oxi','PH','PH_oxi','BZB','BZB_oxi','TBO','TBO_oxi','TIR','TIR_oxi','AP_oxi_d30','APy_oxi_d30','AT_oxi_d30','AT','AP','APy','HEP','PE','PP','PB','PS','SBS','SBSr'};

for j=1:29 %1 = squalane 2=hopane 3=dochn 4=phpn 5= quinolinohopane 6=pyrodinohopane 7=BZB 8 = TBO 9=TIR 10=Asphaltene Thiophene (AT) 11=AP 12= APy 13 = Heptane 14 = PE 15 = PP 16 = PB 17 = PS 18 = SBS 19 = SBSr

    prefix=Namelist{j};
    for k=1:TypeMAT(j)
        if j<23
            [system2] = asphaltCG_make_LAMMPS_NonAsphaltene(prefix);
%         elseif j==22
%             [system2] = asphaltCG_make_AT();
%             system2.mass=707.117;
%         elseif j==25
%             [system2] = asphaltCG_make_AP();
%             system2.mass=574.893;
%         elseif j==26
%             [system2] = asphaltCG_make_APy();
%             system2.mass=888.381;
        elseif j==23 %Heptane
            system2.natoms=2;
            system2.attype={'RC1','RC1'};
            system2.mass=[50 50];
            system2.charge=[0 0 0];
            system2.pos=[0 0 0;0.51 0 0]*10;
            system2.dim=[0.55 0.55 0.55]*10;
            system2.bonds=[1 2];
            system2.bondnames=[1];
            system2.BondPotential=[1 300 5.1];
            system2.angles=[];
            system2.anglenames=[];
            system2.AnglePotential=[];
            system2.dihedrals=[];
            system2.dihedralnames=[];
            system2.DihedralPotential=[];
        elseif j==24 %PE
            [system2]=MakePE(DegreeofPolymerization);
        elseif j==25 %PP
            [system2]=MakePP(DegreeofPolymerization);
        elseif j==26 %PB
            [system2]=MakePB(DegreeofPolymerization);
        elseif j==27 %PS
            [system2]=MakePS(DegreeofPolymerization);
        elseif j==28 %SBS
            [system2]=MakeSBS(SBSDegreeOfPolymerization);
        elseif j==29 %SBSr
            [system2]=MakeSBSradial(SBSRadialDegreeOfPolymerization);
        else 
            error('invalid size')
        end
        TotalMASS=TotalMASS+sum(system2.mass);
        molcounter=molcounter+1;
        counter=counter+1;
        system3=system2;
        maxoffset=BoxSize-system2.dim;
        minboxoffset=[0 0 0];
        
        if j>=14 && RestrictPolymerFlag==1
            masoffset=[0 30 30];
            minboxoffset =[(BoxSize(1)-system2.dim(1))/2 480 480];
        end
        testvar=1;
        while testvar==1
            randoffset=rand(1,3).*maxoffset+minboxoffset;
            [R]=CreateRotationMatrix(180*rand,1-2*rand(1,3));
            system3.pos=system2.pos*R+randoffset
            system3.pos=system2.pos+randoffset;
            if molcounter>1
                D=pdist2(system.pos,system3.pos);
                if min(D)<1.5
                else
                    testvar=0;
                end
            else
                testvar=0;
            end
        end
        MOLPOS(counter,:)=randoffset;

        %%%POS
        %     system.pos=[system.pos;system3.pos];
        system.molid=[system.molid; molcounter*ones(size(system3.pos,1),1)];

        oldnatoms=system.natoms;
        system.natoms=system.natoms+system3.natoms;
        system.attype=[system.attype system3.attype];
        %     system.mass=[system.mass system3.mass];
        system.charge=[system.charge system3.charge];
        system.pos=[system.pos;system3.pos];

        %BONDS
        SYSTEMBONDS=max(system.bondnames);
        if isempty(SYSTEMBONDS)
            SYSTEMBONDS=0;
        end
        system.bonds=[system.bonds; system3.bonds+oldnatoms];
        system.bondnames=[system.bondnames,system3.bondnames+SYSTEMBONDS];
        system3.BondPotential(:,1)=system3.BondPotential(:,1)+SYSTEMBONDS;
        system.BondPotential=[ system.BondPotential; system3.BondPotential];

        %Angles
        SYSTEMANGLES=max(system.anglenames);
        if isempty(SYSTEMANGLES)
            SYSTEMANGLES=0;
        end
        system.angles=[system.angles; system3.angles+oldnatoms];
        system.anglenames=[system.anglenames,system3.anglenames+SYSTEMANGLES];
        if ~isempty(system3.AnglePotential)
            system3.AnglePotential(:,1)=system3.AnglePotential(:,1)+SYSTEMANGLES;
            system.AnglePotential=[ system.AnglePotential; system3.AnglePotential];
        end

        %Dihedrals/Impropers
        SYSTEMIMPROPERS=max(system.impropernames);
        if isempty(SYSTEMIMPROPERS)
            SYSTEMIMPROPERS=0;
        end
        system.impropers=[system.impropers; system3.dihedrals+oldnatoms];
        system.impropernames=[system.impropernames,system3.dihedralnames+SYSTEMIMPROPERS];
        if ~isempty(system3.DihedralPotential)
            system3.DihedralPotential(:,1)=system3.DihedralPotential(:,1)+SYSTEMIMPROPERS;
            system.ImproperPotential=[ system.ImproperPotential; system3.DihedralPotential];
        end
    end

end


% system.impropers=[];
system.bondcoeffs=system.BondPotential;
system.anglecoeffs=system.AnglePotential;
system.impropercoeffs=system.ImproperPotential;
system.dihedralcoeffs=[];
system.paircoeffs=[];
% system.mass=system.mass';


%convert
system.bondcoeffs(:,2)=system.bondcoeffs(:,2)*0.239;
if ~isempty(system.anglecoeffs)
    system.anglecoeffs(:,2)=system.anglecoeffs(:,2)*0.239;
end
if ~isempty(system.impropercoeffs)
    system.impropercoeffs(:,2)=system.impropercoeffs(:,2)*0.239;
end


attypelistcell={'RC1','RC2','RC3','RC4','RC5','SC1','SC2','SC3','SC4','SC5','SC6','TC1','TC2','TC3','TC4','TC5','TC5e','TC6','SN3a','SN4a','TN6','TN6a','TN6d','TP6','SP6','ATA','ATB','ATS','APA','APB','APO','ApA','ApB','ApN'};
system.attypelistcell=attypelistcell;
masslist=[72 72 72 72 72 54 54 54 54 54 54 36 36 36 36 36 36 36 54 54 36 36 36 36 54 36 36 36 36 36 36 36 36 36];
for i=1:length(system.attype)
    try
    system.attypenum(i)=find(strcmp(attypelistcell,system.attype{i}));
    catch
        1==1
    end
end
system.mass=masslist';
if size(system.pos,1)~=system.natoms
    error('check number of atoms vs number of positions')
end
Write_LAMMPS_Data(system,'final_molecule.data')

PairCoeffMARTINI


%paircoeffs for asphaltenes %Lorentz-Berthelot rules
AsphaltSigmas=[7.1 7.1 4 6 6 4 8.4 8.4 4.55]; %ATA ATB ATS APA APB APO
AsphaltEpsilons=[0.33 0.33 0.5 .31 .31 2 .26 .26 .97]; %ATA ATB ATS APA APB APO
for i=1:25
    for j=26:34
        size1=attypelistcell{i};
        % size2=attypelistcell{j};
        size1=size1(1);
        if size1=='R'
            size1=4.7;
        elseif size1=='S'
            size1=4.1;
        elseif size1=='T'
            size1=3.4;
        end
        % size2=size2(1);
        size2=AsphaltSigmas(j-25);
%         eps1=% this will take some work %go self with fit parameter
        eps1=geteps1(i);
        eps2=AsphaltEpsilons(j-25);


        sig=(size1+size2)/2;
        eps=sqrt(eps1*eps2);
        sig=sig;
        eps=eps;
        %Get Eps1
        fprintf(f1,['pair_coeff\t%d\t%d\t%.3f\t%.3f\n'],i,j,eps,sig)
        
    end
end

for i=26:34
    for j=26:34
        if j>=i
            size1=AsphaltSigmas(i-25);
            size2=AsphaltSigmas(j-25);
            sig=(size1+size2)/2;
            eps1=AsphaltEpsilons(i-25);
            eps2=AsphaltEpsilons(j-25);
            eps=sqrt(eps1*eps2);
            sig=sig;
            eps=eps;
            fprintf(f1,['pair_coeff\t%d\t%d\t%.3f\t%.3f\n'],i,j,eps,sig)
        end
    end
end




TotalMASS;
f_MOL_ID=fopen('MoleculeList.txt','a');
str='Total Mass';
for j=length(str)+1:24
    str=[str ' '];
end
fprintf(f_MOL_ID,['\n\n' str '\t' num2str(TotalMASS)  '\n']);
fclose(f_MOL_ID);
fclose(f1);
fclose all;
function [eps]=geteps1(i)

Table1=[  3  4  5  6  6  7   8 10 11 12 13 14 14 15 16 17 17 19
  4  5  6  7  7  8   9 10 10 12 13 14 14 14 15 16 16 17
  5  6  6  7  7  8   9 10 10 11 12 13 13 14 15 16 16 17
  6  7   7 7   7 8   9 10 10 11 11 12 12 13 14 16 16 17
  6  7   7 7   7 8   8  9 10 11 11 11 12 13 14 15 15 16
  7  8   8 8  8  8   8  9  9 11 11 11 12 13 14 15 15 16
  8  9   9 9  8  8   8  8  8 10 11 11 11 12 13 14 14 15
  10 10 10 10 9  9   8  9  9 10 10 11 11 12 13 13 14 15
  11 10 10 10 10 9   8  9  9 10 10 11 11 12 12 13 14 15
  12 12 11 11 11 11 10 10 10 10 10 11 11 12 12 13 13 14
  13 13 12 11 11 11 11 10 10 10 10 11 11 12 12 12 13 14
  14 14 13 12 11 11 11 11 11 11 11 10 10 11 11 12 13 14
  14 14 13 12 12 12 11 11 11 11 11 10 10 11 11 12 13 14
  15 14 14 13 13 13 12 12 12 12 12 11 11 11 11 12 13 13
  16 15 15 14 14 14 13 13 12 12 12 11 11 11 11 11 12 13
  17 16 16 16 15 15 14 13 13 13 12 12 12 12 11 11 11 11
  17 16 16 16 15 15 14 14 14 13 13 13 13 13 12 11 11 10
  19 17 17 17 16 16 15 15 15 14 14 14 14 13 13 11 10 11
];

Table2=[  0.470    5.490     0.430     5.210    0.410     4.840    0.395     4.870    0.365     4.510     0.340    4.020
  0.470    5.240     0.430    4.920     0.410     4.550    0.395     4.540    0.365     4.220     0.340    3.750
  0.470    4.990     0.430    4.680     0.410     4.290    0.395     4.230    0.365     3.930     0.340    3.470
  0.470    4.730     0.430     4.400    0.410     4.040    0.395     3.910    0.365     3.640     0.340    3.220
  0.470    4.480     0.430     4.170    0.410     3.820    0.395     3.670    0.365     3.370     0.340    2.990
  0.470    4.250     0.430     3.970    0.410     3.550    0.395     3.430    0.365     3.100     0.340    2.740
  0.470    4.060     0.430     3.770    0.410     3.310    0.395     3.190    0.365     2.830     0.340    2.500
  0.470    3.880     0.430     3.590    0.410     3.070    0.395     2.940    0.365     2.570     0.340    2.280
  0.470    3.690     0.430     3.380    0.410     2.840    0.395     2.700    0.365     2.340     0.340    2.020
  0.470    3.520     0.430     3.160    0.410     2.600    0.395     2.500    0.365     2.110     0.340    1.770
  0.470    3.390     0.430     2.920    0.410     2.350    0.395     2.310    0.365     1.910     0.340    1.510
  0.470    3.240     0.430     2.770    0.410     2.200    0.395     2.110    0.365     1.750     0.340    1.230
  0.470    3.070     0.430     2.530    0.410     2.040    0.395     1.850    0.366     1.600     0.352    1.050
  0.470    2.960     0.430     2.320    0.410     1.910    0.395     1.670    0.378     1.450     0.366    0.935
  0.470    2.790     0.430     2.160    0.410     1.750    0.398     1.450    0.394     1.290     0.385    0.855
  0.470    2.620     0.430     1.930    0.414     1.560    0.399     1.180    0.404     1.130     0.411    0.801
  0.470    2.440     0.430     1.730    0.418     1.400    0.401     0.983    0.409     1.080     0.424    0.776
  0.485    2.240     0.443     1.480    0.421     1.250    0.404     0.830    0.423     0.981     0.438    0.731
  0.520    2.120     0.473     1.280    0.453     1.150    0.465     0.703    0.484     0.890     0.505    0.680
  0.570    2.000     0.514     1.110    0.496     1.050    0.511     0.659    0.537     0.808     0.561    0.630
  0.620     1.960    0.552     1.010    0.535     0.986    0.524     0.629    0.588     0.763     0.599    0.620
];

% attypelistcell={'RC1','RC2','RC3','RC4','RC5','SC1','SC2','SC3','SC4','SC5','SC6','TC1','TC2','TC3','TC4','TC5','TC5e','TC6','SN3a','SN4a','TN6a','TP6','SP6'};
attypelistcell={'RC1','RC2','RC3','RC4','RC5','SC1','SC2','SC3','SC4','SC5','SC6','TC1','TC2','TC3','TC4','TC5','TC5e','TC6','SN3a','SN4a','TN6','TN6a','TN6d','TP6','SP6'}

Table1Axis={'P6','P5','P4','P3','P2','P1','N6','N5','N4','N3','N2','N1','C6','C5','C4','C3','C2','C1'};
        x=attypelistcell{i};
        x=x(2:3);
        y=attypelistcell{i};
        y=y(2:3);
        X=find(strcmp(Table1Axis,x));
        Y=find(strcmp(Table1Axis,y));
        level=Table1(X,Y);
        %check labels
        size1=attypelistcell{i};
        size2=attypelistcell{i};
        size1=size1(1);
        size2=size2(1);
        sizes=sort([size1 size2]);
        L={};
        x=attypelistcell{i};
        y=attypelistcell{i};
        if isletter(x(end))
            L(end+1)={x(end)};
        else
            L(end+1)={'o'};
        end
        if isletter(y(end))
            L(end+1)={y(end)};
        else
            L(end+1)={'o'};
        end
        L=[L{1} L{2}];
        L=sort(L);
        if strcmp(L,'aa')
            if strcmp(sizes,'RR') || strcmp(sizes,'RS') 
                level=level+2;
                if level>12
                    level=12;
                end
            elseif strcmp(sizes,'SS') || strcmp(sizes,'ST') || strcmp(sizes,'RT') 
                level=level+1;
                if level>12
                    level=12;
                end
            elseif strcmp(sizes,'TT')
                level=level+1;
                if level>12
                    level=12;
                end
            end
        elseif strcmp(L,'dd')
            if strcmp(sizes,'RR') || strcmp(sizes,'RS') 
                level=level+2;
                if level>12
                    level=12;
                end
            elseif strcmp(sizes,'SS') || strcmp(sizes,'ST') || strcmp(sizes,'RT') 
                level=level+1;
                if level>12
                    level=12;
                end
            elseif strcmp(sizes,'TT')
                level=level+1;
                if level>12
                    level=12;
                end
            end
        elseif strcmp(L,'ad')
            if strcmp(sizes,'RR') || strcmp(sizes,'RS') 
                level=level+0;
                if level>12
                    level=12;
                end
            elseif strcmp(sizes,'SS') || strcmp(sizes,'ST') || strcmp(sizes,'RT') 
                level=level-1;
                if level>12
                    level=12;
                end
            elseif strcmp(sizes,'TT')
                level=level-3;
                if level>12
                    level=12;
                end
            end
        elseif strcmp(L,'de')
            if strcmp(sizes,'RR') || strcmp(sizes,'RS') || strcmp(sizes,'SS') || strcmp(sizes,'ST') || strcmp(sizes,'RT') || strcmp(sizes,'TT') 
                level=level-1;
            end
        elseif strcmp(L,'ee')
            level=level+1;
                if level>14
                    level=14;
                end
        elseif strcmp(L,'ae')
            if strcmp(sizes,'RR') || strcmp(sizes,'RS') 
                level=level+2;
            elseif strcmp(sizes,'SS') || strcmp(sizes,'ST') || strcmp(sizes,'RT') || strcmp(sizes,'TT') 
                level=level+1;
            end
        elseif strcmp(L,'eo')
            level=level-1;
        elseif strcmp(L,'oo') || strcmp(L,'ao')|| strcmp(L,'do')
        else
            lerror('check labels')
        end



        size1=attypelistcell{i};
        size2=attypelistcell{i};
        size1=size1(1);
        size2=size2(1);
        
        C=-1;
        if strcmp(size1,size2) %if they are the same
            if strcmp(size1,'R')
                C=1;
            elseif strcmp(size1,'S')
                C=3;
            elseif strcmp(size1,'T')
                C=6;
            end
        else % if they are different
            if strcmp(size1,'R')
                if strcmp(size2,'S')
                    C=2;
                elseif strcmp(size2,'T')
                    C=4;
                end
            elseif strcmp(size1,'S')
                if strcmp(size2,'R')
                    C=2;
                elseif strcmp(size2,'T')
                    C=5;
                end
            elseif strcmp(size1,'T')
                if strcmp(size2,'R')
                    C=4;
                elseif strcmp(size2,'S')
                    C=5;
                end
            end
        end

        eps= Table2(level,C*2)*0.239;
        sig= Table2(level,C*2-1)*10;

end


function [system2]=MakePE(DegreeofPolymerization)
if mod(DegreeofPolymerization(1),2)~=0
    error('Degree of Polymerization for PE must be mulitple of 2');
end
system2.natoms=DegreeofPolymerization(1)/2;
system2.attype=repmat({'RC1'},1,system2.natoms);
system2.mass=56*ones(1,system2.natoms);
system2.charge=zeros(1,system2.natoms);
system2.pos=[[0:5.1:5.1*(system2.natoms-1)]', zeros(system2.natoms,1), zeros(system2.natoms,1)];
system2.dim=[5.1*system2.natoms 7 7];
for i=1:system2.natoms-1
    system2.bonds(i,:)=[i i+1];
    system2.bondnames(i)=[i];
    system2.BondPotential(i,:)=[i 300 5.1];
end

system2.angles=[];
system2.anglenames=[];
system2.AnglePotential=[];
system2.dihedrals=[];
system2.dihedralnames=[];
system2.DihedralPotential=[];

end


function [system2]=MakePP(DegreeofPolymerization)
system2.natoms=DegreeofPolymerization(2);
system2.attype=repmat({'SC3'},1,system2.natoms);
system2.mass=42*ones(1,system2.natoms);
system2.charge=zeros(1,system2.natoms);
system2.pos=[[0:2.96:2.96*(system2.natoms-1)]', zeros(system2.natoms,1), zeros(system2.natoms,1)];
system2.dim=[3*system2.natoms 10 10];
for i=1:system2.natoms-1
    system2.bonds(i,:)=[i i+1];
    system2.bondnames(i)=[i];
    system2.BondPotential(i,:)=[i 300 2.96];
end

system2.angles=[];
system2.anglenames=[];
system2.AnglePotential=[];
system2.dihedrals=[];
system2.dihedralnames=[];
system2.DihedralPotential=[];
end


function [system2]=MakePB(DegreeofPolymerization)
system2.natoms=DegreeofPolymerization(3);
system2.attype=repmat({'RC4'},1,system2.natoms);
system2.mass=54*ones(1,system2.natoms);
system2.charge=zeros(1,system2.natoms);
system2.pos=[[0:4.85:4.85*(system2.natoms-1)]', zeros(system2.natoms,1), zeros(system2.natoms,1)];
system2.dim=[5*system2.natoms 7 7];
for i=1:system2.natoms-1
    system2.bonds(i,:)=[i i+1];
    system2.bondnames(i)=[i];
    system2.BondPotential(i,:)=[i 300 4.85];
end

system2.angles=[];
system2.anglenames=[];
system2.AnglePotential=[];
system2.dihedrals=[];
system2.dihedralnames=[];
system2.DihedralPotential=[];

end


function [system2]=MakePS(DegreeofPolymerization)
system2.natoms=DegreeofPolymerization(4)*4;
system2.attype=repmat({'TC5','TC5','TC5','TC3'},1,system2.natoms/4);
system2.mass=104*ones(1,system2.natoms)/4;
system2.charge=zeros(1,system2.natoms);
i=1;
system2.pos(1:4,:)=[(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']
for i=2:DegreeofPolymerization(4)
    system2.pos([end+1:end+4],:)=[(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']
end
system2.pos=system2.pos*10;
system2.dim=[2.6*DegreeofPolymerization(4) 7 7];

system2.bonds=[1 2;1 3; 2 3; 2 4; 3 4];
system2.bondnames=[1 2 3 4 5];
system2.BondPotential=[1 6000 2.91;2 6000 2.91;3 6000 2.91;4 300 3.2; 5 300 3.2];
system2.dihedrals=[1 2 3 4];
system2.dihedralnames(1)=1;
system2.DihedralPotential(1,:)=[1 100 180];

for i=2:DegreeofPolymerization(4)
    system2.bonds(end+1:end+6,:)=4*(i-1)+[0 4; 1 2;1 3; 2 3; 2 4; 3 4];
    
    system2.BondPotential(end+1:end+6,:)=[1+length(system2.bondnames) 300 2.55; 2+length(system2.bondnames) 6000 2.91;3+length(system2.bondnames) 6000 2.91;4+length(system2.bondnames) 6000 2.91;5+length(system2.bondnames) 300 3.2; 6+length(system2.bondnames) 300 3.2];
    system2.bondnames(end+1:end+6)=[1 2 3 4 5 6]+length(system2.bondnames);
    system2.dihedrals(i,:)=[1 2 3 4]+4*(i-1);
    system2.dihedralnames(i)=i;
    system2.DihedralPotential(i,:)=[i 100 180];
end

system2.angles=[];
system2.anglenames=[];
system2.AnglePotential=[];
end


function [system2]=MakeSBS(SBSDegreeOfPolymerization)
    % PS
    system2.natoms=SBSDegreeOfPolymerization(1)*4;
    natomsS1=system2.natoms;
    system2.attype=repmat({'TC5','TC5','TC5','TC3'},1,system2.natoms/4);
    system2.mass=104*ones(1,system2.natoms)/4;
    system2.charge=zeros(1,system2.natoms);
    i=1;
    system2.pos(1:4,:)=[(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']
    for i=2:SBSDegreeOfPolymerization(1)
        system2.pos([end+1:end+4],:)=[(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']
    end
    system2.pos=system2.pos*10;
    system2.dim=[2.6*SBSDegreeOfPolymerization(1) 7 7];
    
    system2.bonds=[1 2;1 3; 2 3; 2 4; 3 4];
    system2.bondnames=[1 2 3 4 5];
    system2.BondPotential=[1 6000 2.91;2 6000 2.91;3 6000 2.91;4 300 3.2; 5 300 3.2];
    system2.dihedrals=[1 2 3 4];
    system2.dihedralnames(1)=1;
    system2.DihedralPotential(1,:)=[1 100 180];
    
    for i=2:SBSDegreeOfPolymerization(1)
        system2.bonds(end+1:end+6,:)=4*(i-1)+[0 4; 1 2;1 3; 2 3; 2 4; 3 4];
        
        system2.BondPotential(end+1:end+6,:)=[1+length(system2.bondnames) 300 2.55; 2+length(system2.bondnames) 6000 2.91;3+length(system2.bondnames) 6000 2.91;4+length(system2.bondnames) 6000 2.91;5+length(system2.bondnames) 300 3.2; 6+length(system2.bondnames) 300 3.2];
        system2.bondnames(end+1:end+6)=[1 2 3 4 5 6]+length(system2.bondnames);
        system2.dihedrals(i,:)=[1 2 3 4]+4*(i-1);
        system2.dihedralnames(i)=i;
        system2.DihedralPotential(i,:)=[i 100 180];
    end
    
    system2.angles=[];
    system2.anglenames=[];
    system2.AnglePotential=[];


    %PB
    natomsB=SBSDegreeOfPolymerization(2);
    system2.natoms=system2.natoms+natomsB;
    system2.attype=[system2.attype repmat({'RC4'},1,natomsB)];
    system2.mass=[system2.mass 54*ones(1,natomsB)];
    system2.charge=[system2.charge zeros(1,natomsB)];
    posend=max(system2.pos(:,1))+4.85;
    system2.pos=[system2.pos; [posend+[0:4.85:4.85*(natomsB-1)]', zeros(natomsB,1), zeros(natomsB,1)]];
%     system2.dim=[5*system2.natoms 7 7];
%     system2.bonds(end+1,:)=[natomsS1 natomsS1+1] %S to B bond
%     system2.bondnames(end+1)=[natomsS1+1];
%     system2.BondPotential(end+1,:)=[natomsS1+1 300 4.85];
    nbondsS1=size(system2.bonds,1);
    counter=1;
    for i=nbondsS1+1:nbondsS1+natomsB
        system2.bonds(i,:)=[natomsS1+counter-1 natomsS1+counter];
        system2.bondnames(i)=[nbondsS1+counter];
        system2.BondPotential(i,:)=[nbondsS1+counter 300 4.85];
        counter=counter+1;
    end
    
%     system2.angles=[];
%     system2.anglenames=[];
%     system2.AnglePotential=[];
%     system2.dihedrals=[];
%     system2.dihedralnames=[];
%     system2.DihedralPotential=[];

    %PS2

    natomsS2=SBSDegreeOfPolymerization(3)*4;
    system2.natoms=system2.natoms+natomsS2;
    system2.attype=[system2.attype repmat({'TC5','TC5','TC5','TC3'},1,natomsS2/4)];
    system2.mass=[system2.mass 104*ones(1,natomsS2)/4];
    system2.charge=[ system2.charge zeros(1,natomsS2)];
    i=1;
    
    posend=max(system2.pos(:,1))+4.85;
    system2.pos(end+1:end+4,:)=[posend/10+(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']*10
    
    for i=2:SBSDegreeOfPolymerization(3)
        system2.pos([end+1:end+4],:)=[posend/10+2.55/10+(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']*10
    end
%     system2.pos=system2.pos*10;
    system2.dim=[2.6*SBSDegreeOfPolymerization(1)+5*SBSDegreeOfPolymerization(2)+2.6*SBSDegreeOfPolymerization(3) 7 7];
    
%     system2.bonds=[1 2;1 3; 2 3; 2 4; 3 4];
%     system2.bondnames=[1 2 3 4 5];
%     system2.BondPotential=[1 6000 2.91;2 6000 2.91;3 6000 2.91;4 300 3.2; 5 300 3.2];
%     system2.dihedrals=[1 2 3 4];
%     system2.dihedralnames(1)=1;
%     system2.DihedralPotential(1,:)=[1 100 180];
    DihedralCounter=size(system2.dihedrals,1)+1;
    counter=1;
    for i=1+natomsS1+natomsB:natomsS1+natomsB+SBSDegreeOfPolymerization(3)
        system2.bonds(end+1:end+6,:)=natomsS1+natomsB+4*(i-natomsS1-natomsB-1)+[0 4; 1 2;1 3; 2 3; 2 4; 3 4];
        
        system2.BondPotential(end+1:end+6,:)=[1+length(system2.bondnames) 300 2.55; 2+length(system2.bondnames) 6000 2.91;3+length(system2.bondnames) 6000 2.91;4+length(system2.bondnames) 6000 2.91;5+length(system2.bondnames) 300 3.2; 6+length(system2.bondnames) 300 3.2];
        system2.bondnames(end+1:end+6)=[1 2 3 4 5 6]+length(system2.bondnames);
        system2.dihedrals(DihedralCounter,:)=[1 2 3 4]+natomsS1+natomsB+4*(counter-1);
        system2.dihedralnames(DihedralCounter)=DihedralCounter;
        system2.DihedralPotential(DihedralCounter,:)=[DihedralCounter 100 180];
        counter=counter+1;
        DihedralCounter=DihedralCounter+1;
    end

%     system2.bonds=system2.bonds-1;




end


function [system2]=MakeSBSradial(SBSRadialDegreeOfPolymerization)
    SBSDegreeOfPolymerization=SBSRadialDegreeOfPolymerization;
    % PS
    system2.natoms=SBSDegreeOfPolymerization(1)*4;
    natomsS1=system2.natoms;
    system2.attype=repmat({'TC5','TC5','TC5','TC3'},1,system2.natoms/4);
    system2.mass=104*ones(1,system2.natoms);
    system2.charge=zeros(1,system2.natoms);
    i=1;
    system2.pos(1:4,:)=[(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']
    for i=2:SBSDegreeOfPolymerization(1)
        system2.pos([end+1:end+4],:)=[(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']
    end
    system2.pos=system2.pos*10;
    system2.dim=[2.6*SBSDegreeOfPolymerization(1) 7 7];
    
    system2.bonds=[1 2;1 3; 2 3; 2 4; 3 4];
    system2.bondnames=[1 2 3 4 5];
    system2.BondPotential=[1 6000 2.91;2 6000 2.91;3 6000 2.91;4 300 3.2; 5 300 3.2];
    system2.dihedrals=[1 2 3 4];
    system2.dihedralnames(1)=1;
    system2.DihedralPotential(1,:)=[1 100 180];
    
    for i=2:SBSDegreeOfPolymerization(1)
        system2.bonds(end+1:end+6,:)=4*(i-1)+[0 4; 1 2;1 3; 2 3; 2 4; 3 4];
        
        system2.BondPotential(end+1:end+6,:)=[1+length(system2.bondnames) 300 2.55; 2+length(system2.bondnames) 6000 2.91;3+length(system2.bondnames) 6000 2.91;4+length(system2.bondnames) 6000 2.91;5+length(system2.bondnames) 300 3.2; 6+length(system2.bondnames) 300 3.2];
        system2.bondnames(end+1:end+6)=[1 2 3 4 5 6]+length(system2.bondnames);
        system2.dihedrals(i,:)=[1 2 3 4]+4*(i-1);
        system2.dihedralnames(i)=i;
        system2.DihedralPotential(i,:)=[i 100 180];
    end
    
    system2.angles=[];
    system2.anglenames=[];
    system2.AnglePotential=[];


    %PB1
    natomsB=SBSDegreeOfPolymerization(2);
    system2.natoms=system2.natoms+natomsB;
    system2.attype=[system2.attype repmat({'RC4'},1,natomsB)];
    system2.mass=[system2.mass 54*ones(1,natomsB)];
    system2.charge=[system2.charge zeros(1,natomsB)];
    posend=max(system2.pos(:,1))+4.85;
    system2.pos=[system2.pos; [posend+[0:4.85:4.85*(natomsB-1)]', zeros(natomsB,1), zeros(natomsB,1)]];
%     system2.dim=[5*system2.natoms 7 7];
%     system2.bonds(end+1,:)=[natomsS1 natomsS1+1] %S to B bond
%     system2.bondnames(end+1)=[natomsS1+1];
%     system2.BondPotential(end+1,:)=[natomsS1+1 300 4.85];
    nbondsS1=size(system2.bonds,1);
    counter=1;
    for i=nbondsS1+1:nbondsS1+natomsB
        system2.bonds(i,:)=[natomsS1+counter-1 natomsS1+counter];
        system2.bondnames(i)=[nbondsS1+counter];
        system2.BondPotential(i,:)=[nbondsS1+counter 300 4.85];
        counter=counter+1;
    end
    
%     system2.angles=[];
%     system2.anglenames=[];
%     system2.AnglePotential=[];
%     system2.dihedrals=[];
%     system2.dihedralnames=[];
%     system2.DihedralPotential=[];


    %PB2
    MidBead=ceil(SBSDegreeOfPolymerization(2)/2);
    MidBead2=floor(SBSDegreeOfPolymerization(3)/2)+1;
    natomsB1=natomsB;
    natomsB=SBSDegreeOfPolymerization(3);
    natomsS1old=natomsS1;
    natomsS1=system2.natoms;
    system2.natoms=system2.natoms+natomsB;
    system2.attype=[system2.attype repmat({'RC4'},1,natomsB)];
    system2.mass=[system2.mass 54*ones(1,natomsB)];
    system2.charge=[system2.charge zeros(1,natomsB)];
    posendx=system2.pos(natomsS1old+MidBead,1);
    posendy=system2.pos(natomsS1old+MidBead,2);
    system2.pos=[system2.pos; posendx*ones(natomsB,1), [4.85:4.85:4.85*(MidBead2-1) -4.85:-4.85:-4.85*(natomsB-MidBead2+1)]',zeros(natomsB,1)]; 
%     system2.pos=[system2.pos; [posend+[0:4.85:4.85*(natomsB-1)]', zeros(natomsB,1), zeros(natomsB,1)]];
%     system2.dim=[5*system2.natoms 7 7];
%     system2.bonds(end+1,:)=[natomsS1 natomsS1+1] %S to B bond
%     system2.bondnames(end+1)=[natomsS1+1];
%     system2.BondPotential(end+1,:)=[natomsS1+1 300 4.85];
    nbondsS1old=nbondsS1;
    nbondsS1=size(system2.bonds,1);
    TOPBBEAD=system2.natoms;
    %Middle Bead to Bond to
    
    
%     system2.pos(end-SBSDegreeOfPolymerization(3)+1:end,2)=system2.pos(end-SBSDegreeOfPolymerization(3)+1:end,2)-5;
%     system2.pos(end-SBSDegreeOfPolymerization(3)+1:end-SBSDegreeOfPolymerization(3)+MidBead2-1,2)=system2.pos(end-SBSDegreeOfPolymerization(3)+1:end-SBSDegreeOfPolymerization(3)+MidBead2-1,2)+10;

    counter=1;
    for i=nbondsS1+1:nbondsS1+natomsB
        system2.bonds(i,:)=[natomsS1+counter-1 natomsS1+counter];
        if i==nbondsS1+1;
           system2.bonds(i,:)= [system2.natoms-natomsB-natomsB1+MidBead natomsS1+counter]
        elseif i==nbondsS1+MidBead2
            system2.bonds(i,:)= [system2.natoms-natomsB-natomsB1+MidBead natomsS1+counter]
        end
%         system2.pos(natomsS1+counter,:)
        system2.bondnames(i)=[nbondsS1+counter];
        system2.BondPotential(i,:)=[nbondsS1+counter 300 4.85];
        counter=counter+1;
    end


    %PS2
    
    natomsS2=SBSDegreeOfPolymerization(4)*4;
    system2.natoms=system2.natoms+natomsS2;
    system2.attype=[system2.attype repmat({'TC5','TC5','TC5','TC3'},1,natomsS2/4)];
    system2.mass=[system2.mass 104*ones(1,natomsS2)];
    system2.charge=[ system2.charge zeros(1,natomsS2)];
    i=1;
    
    posend=max(system2.pos(:,1))+4.85;
    system2.pos(end+1:end+4,:)=[posend/10+(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']*10
    
    for i=2:SBSDegreeOfPolymerization(4)
        system2.pos([end+1:end+4],:)=[posend/10+2.55/10+(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]',[.1454 0 .1454*2 .1454]']*10
    end
%     system2.pos=system2.pos*10;
    
    
%     system2.bonds=[1 2;1 3; 2 3; 2 4; 3 4];
%     system2.bondnames=[1 2 3 4 5];
%     system2.BondPotential=[1 6000 2.91;2 6000 2.91;3 6000 2.91;4 300 3.2; 5 300 3.2];
%     system2.dihedrals=[1 2 3 4];
%     system2.dihedralnames(1)=1;
%     system2.DihedralPotential(1,:)=[1 100 180];
    DihedralCounter=size(system2.dihedrals,1)+1;
    counter=1;
    for i=1+natomsS1+natomsB:natomsS1+natomsB+SBSDegreeOfPolymerization(4)
        system2.bonds(end+1:end+6,:)=natomsS1+natomsB+4*(i-natomsS1-natomsB-1)+[0 4; 1 2;1 3; 2 3; 2 4; 3 4];
        if i==1+natomsS1+natomsB
            system2.bonds(end-5,:)=[system2.natoms-natomsS2-natomsB natomsS1+natomsB+4*(i-natomsS1-natomsB-1)+4]
        end
        
        system2.BondPotential(end+1:end+6,:)=[1+length(system2.bondnames) 300 2.55; 2+length(system2.bondnames) 6000 2.91;3+length(system2.bondnames) 6000 2.91;4+length(system2.bondnames) 6000 2.91;5+length(system2.bondnames) 300 3.2; 6+length(system2.bondnames) 300 3.2];
        system2.bondnames(end+1:end+6)=[1 2 3 4 5 6]+length(system2.bondnames);
        system2.dihedrals(DihedralCounter,:)=[1 2 3 4]+natomsS1+natomsB+4*(counter-1);
        system2.dihedralnames(DihedralCounter)=DihedralCounter;
        system2.DihedralPotential(DihedralCounter,:)=[DihedralCounter 100 180];
        counter=counter+1;
        DihedralCounter=DihedralCounter+1;
    end

%%%%%%%%%%%%%%%%%%%%%%PS3_1
    natomssum=system2.natoms;
    natomsS2=SBSDegreeOfPolymerization(4)*4;
    system2.natoms=system2.natoms+natomsS2;
    system2.attype=[system2.attype repmat({'TC5','TC5','TC5','TC3'},1,natomsS2/4)];
    system2.mass=[system2.mass 104*ones(1,natomsS2)];
    system2.charge=[ system2.charge zeros(1,natomsS2)];
    i=1;
    
%     posend=max(system2.pos(:,1))+4.85;
%     posend(1)=posendx;
%     posend(2)=posendy;
%     posend(3)=0;
    posend=posendx;
    posendy=(SBSDegreeOfPolymerization(3)/2+1)*0.485;
    system2.pos(end+1:end+4,:)=[posend/10+(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]'+posendy,[.1454 0 .1454*2 .1454]']*10;
    R = rotx(180);
        system2.pos(end-3:end,:)=system2.pos(end-3:end,:)*R;
    
    for i=2:SBSDegreeOfPolymerization(4)
        system2.pos([end+1:end+4],:)=[posend/10+2.55/10+(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]'+posendy,[.1454 0 .1454*2 .1454]']*10
        R = rotx(180);
        system2.pos(end-3:end,:)=system2.pos(end-3:end,:)*R;
    end
%     system2.pos=system2.pos*10;
    
    
%     system2.bonds=[1 2;1 3; 2 3; 2 4; 3 4];
%     system2.bondnames=[1 2 3 4 5];
%     system2.BondPotential=[1 6000 2.91;2 6000 2.91;3 6000 2.91;4 300 3.2; 5 300 3.2];
%     system2.dihedrals=[1 2 3 4];
%     system2.dihedralnames(1)=1;
%     system2.DihedralPotential(1,:)=[1 100 180];
    DihedralCounter=size(system2.dihedrals,1)+1;
    counter=1;
    
    for i=natomssum:natomssum+SBSDegreeOfPolymerization(4)-1
        system2.bonds(end+1:end+6,:)=natomssum+4*(i-natomssum)+[0 4; 1 2;1 3; 2 3; 2 4; 3 4];
        if i==natomssum
%             system2.bonds(end-5,:)=[natomssum natomssum-1]
            system2.bonds(end-5,:)=[natomssum+4 TOPBBEAD]
        end
        
        system2.BondPotential(end+1:end+6,:)=[1+length(system2.bondnames) 300 2.55; 2+length(system2.bondnames) 6000 2.91;3+length(system2.bondnames) 6000 2.91;4+length(system2.bondnames) 6000 2.91;5+length(system2.bondnames) 300 3.2; 6+length(system2.bondnames) 300 3.2];
        system2.bondnames(end+1:end+6)=[1 2 3 4 5 6]+length(system2.bondnames);
        system2.dihedrals(DihedralCounter,:)=[1 2 3 4]+natomsS1+natomsB+4*(counter-1);
        system2.dihedralnames(DihedralCounter)=DihedralCounter;
        system2.DihedralPotential(DihedralCounter,:)=[DihedralCounter 100 180];
        counter=counter+1;
        DihedralCounter=DihedralCounter+1;
    end


    %%%%%%%%%%%%PS_4
    %%%%%%%%%%%%%%%%%%%%%%PS3_1
    natomssum=system2.natoms;
    natomsS2=SBSDegreeOfPolymerization(4)*4;
    system2.natoms=system2.natoms+natomsS2;
    system2.attype=[system2.attype repmat({'TC5','TC5','TC5','TC3'},1,natomsS2/4)];
    system2.mass=[system2.mass 104*ones(1,natomsS2)];
    system2.charge=[ system2.charge zeros(1,natomsS2)];
    i=1;
    
%     posend=max(system2.pos(:,1))+4.85;
%     posend(1)=posendx;
%     posend(2)=posendy;
%     posend(3)=0;
    posend=posendx;
    posendy=(SBSDegreeOfPolymerization(3)/2+1)*0.485;
    system2.pos(end+1:end+4,:)=[posend/10+(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]'+posendy,[.1454 0 .1454*2 .1454]']*10;
    R = rotx(0);
        system2.pos(end-3:end,:)=system2.pos(end-3:end,:)*R;
    
    for i=2:SBSDegreeOfPolymerization(4)
        system2.pos([end+1:end+4],:)=[posend/10+2.55/10+(0.255*(i-1))+[0 0 0 0]', [.285+.252 .285 .285 0]'+posendy,[.1454 0 .1454*2 .1454]']*10
        R = rotx(0);
        system2.pos(end-3:end,:)=system2.pos(end-3:end,:)*R;
    end
%     system2.pos=system2.pos*10;
    
    
%     system2.bonds=[1 2;1 3; 2 3; 2 4; 3 4];
%     system2.bondnames=[1 2 3 4 5];
%     system2.BondPotential=[1 6000 2.91;2 6000 2.91;3 6000 2.91;4 300 3.2; 5 300 3.2];
%     system2.dihedrals=[1 2 3 4];
%     system2.dihedralnames(1)=1;
%     system2.DihedralPotential(1,:)=[1 100 180];
    DihedralCounter=size(system2.dihedrals,1)+1;
    counter=1;
    
    for i=natomssum:natomssum+SBSDegreeOfPolymerization(4)-1
        system2.bonds(end+1:end+6,:)=natomssum+4*(i-natomssum)+[0 4; 1 2;1 3; 2 3; 2 4; 3 4];
        if i==natomssum
%             system2.bonds(end-5,:)=[natomssum natomssum-1]
            system2.bonds(end-5,:)=[natomssum+4 TOPBBEAD-SBSDegreeOfPolymerization(3)+1]
        end
        
        system2.BondPotential(end+1:end+6,:)=[1+length(system2.bondnames) 300 2.55; 2+length(system2.bondnames) 6000 2.91;3+length(system2.bondnames) 6000 2.91;4+length(system2.bondnames) 6000 2.91;5+length(system2.bondnames) 300 3.2; 6+length(system2.bondnames) 300 3.2];
        system2.bondnames(end+1:end+6)=[1 2 3 4 5 6]+length(system2.bondnames);
        system2.dihedrals(DihedralCounter,:)=[1 2 3 4]+natomsS1+natomsB+4*(counter-1);
        system2.dihedralnames(DihedralCounter)=DihedralCounter;
        system2.DihedralPotential(DihedralCounter,:)=[DihedralCounter 100 180];
        counter=counter+1;
        DihedralCounter=DihedralCounter+1;
    end

    system2.pos(:,1)=system2.pos(:,1)/3;
    system2.pos(:,2)=system2.pos(:,2)/3;
    system2.pos(:,2)=system2.pos(:,2)-min(system2.pos(:,2));
    system2.pos(:,3)=system2.pos(:,3)-min(system2.pos(:,3));
%     system2.dim=[2.6*SBSDegreeOfPolymerization(1)+5*SBSDegreeOfPolymerization(2)+2.6*SBSDegreeOfPolymerization(4) 5*SBSDegreeOfPolymerization(3) 7];
    system2.dim=max(system2.pos)+ [2 2 2];
%     system2.bonds=system2.bonds-1;




end

