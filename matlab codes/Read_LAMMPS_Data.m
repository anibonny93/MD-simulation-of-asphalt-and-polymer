function [system]=Read_LAMMPS_Data(file_name,dumpname)

fid=fopen(file_name);

y='NO';
typeflag=0;
stepcount=0;
readxflag=1;
while strcmp(y,'NO')
    if readxflag==1
        x=fgetl(fid);
    else
        readxflag=1;
    end
    if isnumeric(x)
        if x==-1
            y='YES';
            break
        end
    end
    
    %readnatoms
    if endsWith(x,' atoms')
        x=x(1:find(x=='a',1)-1);
        natoms=str2num(x);
    end
    
    %boxsize
    if endsWith(x,'xhi')
        x=x(1:find(x=='x',1)-1);
        x2=str2num(x);
        x=fgetl(fid);
        x=x(1:find(x=='y',1)-1);
        y2=str2num(x);
        ystore=y2;
        x=fgetl(fid);
        x=x(1:find(x=='z',1)-1);
        z2=str2num(x);
        
        dim=[x2(2)-x2(1),y2(2)-y2(1),z2(2)-z2(1)];
        lammpsdim=[x2(1),y2(1),z2(1);x2(2),y2(2),z2(2)];
    end
    
    
    %atomtypes and positions, only works with full style for now
    if startsWith(x,'Atoms')
        if ~startsWith(x,'Atoms # full')
            warning('Only works with full style. It may error or produce nonsense if another style is used')
        end
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                pos(xvec(1),:)=xvec(5:7);
                attype(xvec(1))=xvec(3);
                charge(xvec(1))=xvec(4);
                atid(xvec(1))=xvec(1);
                molid(xvec(1))=xvec(2);
            end
        end
    elseif startsWith(x,'Atoms') %Libpargen style
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                pos(xvec(1),:)=xvec(5:7);
                attype(xvec(1))=xvec(3);
                charge(xvec(1))=xvec(4);
            end
        end
    end
    
    %pair coeffs
    if startsWith(x,'Pair Coeffs')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if ~exist('paircoeffs','var')
                    paircoeffs=xvec(2:3);
                else
                    paircoeffs(end+1,:)=xvec(2:3);
                end
            end
        end
    end
        
    %bonds
    if startsWith(x,'Bonds')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if ~exist('bonds','var')
                    bonds=xvec(3:4);
                else
                    bonds(end+1,:)=xvec(3:4);
                end
                if ~exist('bondtypes','var')
                    bondtypes=xvec(2);
                else
                    bondtypes(end+1,:)=xvec(2);
                end
            end
        end
    end
    
    %bond coeffs
    if startsWith(x,'Bond Coeffs')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if ~exist('bondcoeffs','var')
                    bondcoeffs=xvec(2:3);
                else
                    bondcoeffs(end+1,:)=xvec(2:3);
                end
            end
        end
    end
    
    %angles
    if startsWith(x,'Angles')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isnumeric(x)
                    y2='YES';
                    break
            elseif isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if ~exist('angles','var')
                    angles=xvec(3:5);
                else
                    angles(end+1,:)=xvec(3:5);
                end
                if ~exist('angletypes','var')
                    angletypes=xvec(2);
                else
                    angletypes(end+1,:)=xvec(2);
                end
            end
        end
    end
    
    if isnumeric(x)
        if x==-1
            break
        end
    end
    
    %angle coeffs
    if startsWith(x,'Angle Coeffs')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if ~exist('anglecoeffs','var')
                    anglecoeffs=xvec(2:3);
                else
                    anglecoeffs(end+1,:)=xvec(2:3);
                end
            end
        end
    end
            
    %dihedrals
    if startsWith(x,'Dihedrals')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if ~exist('dihedrals','var')
                    dihedrals=xvec(3:6);
                else
                    dihedrals(end+1,:)=xvec(3:6);
                end
                if ~exist('dihedraltypes','var')
                    dihedraltypes=xvec(2);
                else
                    dihedraltypes(end+1,:)=xvec(2);
                end
            end
        end
    end
    
    %dihedral coeffs
    if startsWith(x,'Dihedral Coeffs')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if ~exist('dihedralcoeffs','var')
                    dihedralcoeffs=xvec(2:5);
                else
                    dihedralcoeffs(end+1,:)=xvec(2:5);
                end
            end
        end
    end
    
    %impropers
    if startsWith(x,'Impropers')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isnumeric(x)
                if x==-1
                    y='YES';
                    break
                end
            end
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if ~exist('impropers','var')
                    impropers=xvec(3:6);
                else
                    impropers(end+1,:)=xvec(3:6);
                end
            end
        end
    end
    if isnumeric(x)
        if x==-1
            y='YES';
            break
        end
    end
    %improper coeffs
    if startsWith(x,'Improper Coeffs')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if ~exist('impropercoeffs','var')
                    impropercoeffs=xvec(2:end);
                else
                    impropercoeffs(end+1,:)=xvec(2:end);
                end
            end
        end
    end
    
    %masses to get atom type
    if startsWith(x,'Masses')
        x=fgetl(fid);
        y2='NO';
        while strcmp(y2,'NO')
            x=fgetl(fid);
            if contains(x,'#')
                x=split(x,'#');
                x=x{1};
            end
            if isempty(str2num(x))
                y2='YES';
                break
            else
                xvec=str2num(x);
                if round(xvec(2))==14
                    typematch(xvec(1))='N';
                elseif round(xvec(2))==12
                    typematch(xvec(1))='C';
                elseif round(xvec(2))==1
                    typematch(xvec(1))='H';
                elseif round(xvec(2))==16
                    typematch(xvec(1))='O';
                elseif round(xvec(2))==32
                    typematch(xvec(1))='S';
                elseif round(xvec(2))==72
                    typematch=['A','B','C'];
                else
%                     error('atom type not yet defined')
                    typematch(xvec(1))='Z';
                end
                mass(xvec(1))=xvec(2);
            end
        end
    end
    
        
    
    
    
    
%     %Timestep
%     if startsWith(x,'ITEM: TIMESTEP')
%         x=fgetl(fid);
%         if exist('timesteps','var')
%             timesteps(end+1)=str2num(x);
%         else
%             timesteps=str2num(x);
%         end
%         stepcount=stepcount+1
%     end
%     
%     % natoms
%     if strcmp(x,'ITEM: NUMBER OF ATOMS')
%         x=fgetl(fid);
%         if exist('natoms','var')
%             if natoms~=str2num(x)
%                 error('number of atoms changed')
%             end
%         else
%             natoms=str2num(x);
%         end
%     end
%     
%     %boxsize
%     if strcmp(x,'ITEM: BOX BOUNDS pp pp pp')
%         x=fgetl(fid);
%         xdim=str2num(x);
%         xdim=xdim(2)-xdim(1);
%         x=fgetl(fid);
%         ydim=str2num(x);
%         ydim=ydim(2)-ydim(1);
%         x=fgetl(fid);
%         zdim=str2num(x);
%         zdim=zdim(2)-zdim(1);
%         if exist('dim','var')
%             if sum(dim==[xdim ydim zdim])~=3
%                 warning('box size changes, using the last box size')
%             end
%         end
%         dim=[xdim ydim zdim];
%     end
%     
%     %atom types and positions
%     if strcmp(x,'ITEM: ATOMS id type xs ys zs')
%         y2='NO';
%         atomcount=0;
%         while strcmp(y2,'NO')
%             
%             x=fgetl(fid);
%             if isnumeric(x)
%                 if x==-1
%                     break
%                 end
%             end
%             if startsWith(x,'ITEM') || isempty(x)
%                 y2='YES';
%                 readxflag=0;
%                 break
%             end
% 
% 
%             xvec=str2num(x);
%             atomcount=atomcount+1;
%             if typeflag==0
%                 if ~exist('attype','var')
%                     attype(1)=xvec(2);
%                 else
%                     attype(end+1)=xvec(2);
%                 end
%             end
%             if ~exist('pos','var')
%                 pos(1,:,1)=xvec(3:5).*dim;
%             else
%                 pos(atomcount,:,stepcount)=xvec(3:5).*dim;
%             end
%         end
%     end
%     
    
    
    
end
%move to be centered on zero
% pos(:,1)=(pos(:,1)-x2(1))-dim(1)/2;
% pos(:,2)=(pos(:,2)-ystore(1))-dim(2)/2;
% pos(:,3)=(pos(:,3)-z2(1))-dim(3)/2;

% pos(:,1,:)=pos(:,1,:)-dim(1)/2;
% pos(:,2,:)=pos(:,2,:)-dim(2)/2;
% pos(:,3,:)=pos(:,3,:)-dim(3)/2;



system.natoms=natoms;
system.dim=dim;
system.lammpsdim=lammpsdim;

%convert types to letters
attypenum=attype;
attype=attype';
for k=1:length(typematch)
    system.attype(attype==k)=typematch(k);
end
% system.attype(attype==1)='N';
% system.attype(attype==2)='C';
% system.attype(attype==3)='O';
% system.attype(attype==4)='O';
% system.attype(attype==5)='H';
% system.attype(attype==6)='H';
% system.attype(attype==7)='C';
% system.attype(attype==8)='C';
% system.attype(attype==9)='H';
% system.attype(attype==10)='C';



system.attype=system.attype';
system.attypenum=attypenum';
system.atid=atid';
system.molid=molid';
system.mass=mass';

%Check Mass matches with atoms, Ligpargen is different than moltemplate. 
for k=1:size(system.attype,1)
    if strcmp(typematch,'ABC')
        system.mass(k)=72;
    else
        if system.attype(k)=='N'
           system.mass(k)=14.0070;
        elseif system.attype(k)=='C'
           system.mass(k)=12.0110;
           elseif system.attype(k)=='H'
           system.mass(k)=1.008;
           elseif system.attype(k)=='O'
           system.mass(k)=15.999;
           elseif system.attype(k)=='S'
           system.mass(k)=32.065;
        else
            error('atom type mismatched in selecting mass')
        end
    end
end

if exist('paircoeffs','var')
    system.paircoeffs=paircoeffs;
end
% system.timesteps=timesteps;
system.pos=pos;
system.charge=charge;
system.bonds=bonds;
system.bondtypes=bondtypes;
if exist('bondcoeffs','var')
    system.bondcoeffs=bondcoeffs;
end
system.angles=angles;
if exist('anglecoeffs','var')
    system.anglecoeffs=anglecoeffs;
    system.angletypes=angletypes;
end
if exist('dihedrals','var')
    system.dihedrals=dihedrals;
    
end
if exist('dihedralcoeffs','var')
    system.dihedralcoeffs=dihedralcoeffs;
    system.dihedraltypes=dihedraltypes;
end
if exist('impropers','var')
    system.impropers=impropers;
end
if exist('impropercoeffs','var')
    system.impropercoeffs=impropercoeffs;
end

if exist('dumpname','var')
    fid2=fopen(dumpname);
    TEST=1;
    timestep=0;
    while TEST==1;
        x=fgetl(fid2);
        if isnumeric(x)
            if x==-1
                TEST=2;
                return
            end
        end
        if startsWith(x,'ITEM: TIMESTEP')
            x=fgetl(fid2);
            xvec=str2num(x);
            timestep=timestep+1;
            if mod(timestep,10)==0
                timestep
            end
            system.timestep(timestep)=xvec;
        end
        if startsWith(x,'ITEM: NUMBER OF ATOMS')
            x=fgetl(fid2);
            xvec=str2num(x);
            numatom=xvec;
        end
        if startsWith(x,'ITEM: BOX BOUNDS')
            x=fgetl(fid2);
            xvec=str2num(x);
            system.lammpsdim(:,1,timestep)=xvec';
            x=fgetl(fid2);
            xvec=str2num(x);
            system.lammpsdim(:,2,timestep)=xvec';
            x=fgetl(fid2);
            xvec=str2num(x);
            system.lammpsdim(:,3,timestep)=xvec';
            system.lammpsdim1D(timestep,:)=abs(system.lammpsdim(2,:,timestep)-system.lammpsdim(1,:,timestep));
        end
        if startsWith(x,'ITEM: ATOMS id type xs ys zs')
            for i=1:numatom
                x=fgetl(fid2);
                xvec=str2num(x);
                system.pos(xvec(1),:,timestep)=xvec(3:5).*system.lammpsdim1D(timestep,:);
            end
        end
    end
end
        
        
        