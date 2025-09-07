function [system2] = asphaltCG_make_LAMMPS_NonAsphaltene(prefix);
if isunix
    path=['/mnt/Shared_Data/anic/LONI_output/New_files/Kspace_pppm_long/CGing/Oxidized/' prefix '/'];
else
    path=['S:/PlasticRoads/CoarseGraining/' prefix '/'];
end
run([path prefix '_config.m'])


system.attype=Beadtype;
system.bonds=Bonds;
system.masses=masses;

system.Angles=Angles;
system.Dihedrals=Dihedrals;
system.bondcoeffs=bondcoeffs;
system.anglecoeffs=anglecoeffs;
system.dihedralcoeffs=dihedralcoeffs;



clear system2
% charge=0.6
% epsilon=0.35
% sigma=3.9
run([path prefix '_config.m'])
system2.natoms=length(Beadtype);
system2.attype=Beadtype;

system2.mass=masses;
system2.charge=zeros(1,system2.natoms);
system2.pos=pos;
system2.dim=MAXDIM;

%%% Bonds
system2.bonds=Bonds;
for i=1:size(Bonds,1)
    system2.bondnames(i)=i;
end
for i=1:size(Bonds,1)
    system2.BondPotential(i,:)=[i bondcoeffs(i,:)];
end

%%% Angles
system2.angles=Angles;
for i=1:size(Angles,1)
    system2.anglenames(i)=i;
end
for i=1:size(Angles,1)
    system2.AnglePotential(i,:)=[i anglecoeffs(i,:)];
end
if size(Angles,1)==0
    system2.anglenames=[];
    system2.AnglePotential=[];
end

%%% Dihedrals
system2.dihedrals=Dihedrals;
for i=1:size(Dihedrals,1)
    system2.dihedralnames(i)=i;
end
for i=1:size(Dihedrals,1)
    system2.DihedralPotential(i,:)=[i dihedralcoeffs(i,:)];
end
if size(Dihedrals,1)==0
    system2.dihedralnames=[];
    system2.DihedralPotential=[];
end



end


function  [R]=rotx(a)
R=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
end
function [R]=roty(a)
R=[cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
end
function [R]=rotz(a)
R=[cos(a) -sin(a) 0;sin(a) cos(a) 0; 0 0 1];
end
