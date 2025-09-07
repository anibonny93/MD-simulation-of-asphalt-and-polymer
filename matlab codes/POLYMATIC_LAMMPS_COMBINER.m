 function [] = POLYMATIC_LAMMPS_COMBINER(filelist)
%This function takes lammps files in a cell array and creates new lammps
%files for use with polymatic. One file is created for each element in
%filelist. p
if ~iscell(filelist)
    error('input must be a cell array of file names including full paths')
end

atomtypes=0;
bondtypes=0;
angletypes=0;
dihedraltypes=0;
impropertypes=0;

OLDatomtypes=0;
OLDbondtypes=0;
OLDangletypes=0;
OLDdihedraltypes=0;
OLDimpropertypes=0;

for i=1:numel(filelist)
    filein=filelist{i};
    %     if endsWith(filein,'.lmp')
    %         fileout=[filein(1:end-4) '_mixed.lmp'];
    %         fileouttemp=[filein(1:end-4) '_mixed.lmptemp'];
    %     elseif endsWith(filein,'.data')
    %         fileout=[filein(1:end-4) '_mixed.data'];
    %         fileouttemp=[filein(1:end-4) '_mixed.datatemp'];
    %     else
    %         error('files must end with .lmp or .data')
    %     end
    fid_in=fopen(filein,'r');
    %     fid_outtemp=fopen(fileouttemp,'w');
    
    x=fgetl(fid_in);
    %     fprintf(fid_outtemp,'temporary lammps file for polymatic mixing. Can delete after code is run\n')
    
    testflag=0;
    while testflag==0
        x=fgetl(fid_in);
        if isnumeric(x)
            if x==-1
                testflag=0;
                OLDatomtypes=OLDatomtypes+atomtypes;
                OLDbondtypes=OLDbondtypes+bondtypes;
                OLDangletypes=OLDangletypes+angletypes;
                OLDdihedraltypes=OLDdihedraltypes+dihedraltypes;
                OLDimpropertypes=OLDimpropertypes+impropertypes;
                break
            end
        end
        %Don't need all this because I just need to renumber everything
        %firse time through
        %         if isempty(x)
        %             fprintf(fid_outtemp,'\n');
        %         end
        %         if endsWith(x,'atoms') || endsWith(x,'bonds') || endsWith(x,'angles') || endsWith(x,'dihedrals') || endsWith(x,'impropers')
        %             fprintf(fid_outtemp,[x '\n']);
        %         end
        %
        if endsWith(x,'atom types')
            xsplit=split(x);
            atomtypes=str2num(xsplit{2});
            %                     fprintf(fid_outtemp,['    ' num2str(atomtypes) ' atom types\n']);
            
        end
        if endsWith(x,'bond types')
            xsplit=split(x);
            bondtypes=str2num(xsplit{2});
            %                     fprintf(fid_outtemp,['    ' num2str(atomtypes) ' bond types\n']);
        end
        if endsWith(x,'angle types')
            xsplit=split(x);
            angletypes=str2num(xsplit{2});
            %                     fprintf(fid_outtemp,['    ' num2str(atomtypes) ' angle types\n']);
        end
        if endsWith(x,'dihedral types')
            xsplit=split(x);
            dihedraltypes=str2num(xsplit{2});
            %                     fprintf(fid_outtemp,['    ' num2str(atomtypes) ' dihedral types\n']);
        end
        if endsWith(x,'improper types')
            xsplit=split(x);
            impropertypes=str2num(xsplit{2});
            %                     fprintf(fid_outtemp,['    ' num2str(atomtypes) ' improper types\n']);
        end
        %
        %         if endsWith(x,'hi')
        %             fprintf(fid_outtemp,[x '\n']);
        %         end
        
        if startsWith(x,'Masses') %Write base masses file
            %             fprintf(fid_outtemp,['Masses \n\n']);
            x=fgetl(fid_in);
            x=fgetl(fid_in);
            xsplit=split(x);
            if exist('fMASSES','var')~=1
                fMASSES=fopen('MASSES','w');
                fprintf(fMASSES,'Masses\n\n');
                fprintf(fMASSES,['    ' num2str(str2num(xsplit{2})+OLDatomtypes) ' ' num2str(xsplit{3}) '  \n'])
            else
                fprintf(fMASSES,['    ' num2str(str2num(xsplit{2})+OLDatomtypes) ' ' num2str(xsplit{3}) '  \n'])
            end
            testflag2=0;
            while testflag2==0
                x=fgetl(fid_in);
                xsplit=split(x);
                if ~isempty(x)
                    fprintf(fMASSES,['    ' num2str(str2num(xsplit{2})+OLDatomtypes) ' ' num2str(xsplit{3}) '  \n'])
                else
                    testflag2=1;
                    %fprintf(fMASSES,'\n');
                end
            end
        end
        if startsWith(x,'Pair Coeffs') %Write base paircoeffs
            x=fgetl(fid_in);
            x=fgetl(fid_in);
            xsplit=split(x);
            if exist('fPCOEFF','var')~=1
                fPCOEFF=fopen('PCOEFF','w');
                fprintf(fPCOEFF,'Pair Coeffs\n\n');
                fprintf(fPCOEFF,['    ' num2str(str2num(xsplit{2})+OLDatomtypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  \n'])
            else
                fprintf(fPCOEFF,['    ' num2str(str2num(xsplit{2})+OLDatomtypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  \n'])
            end
            testflag2=0;
            while testflag2==0
                x=fgetl(fid_in);
                xsplit=split(x);
                if ~isempty(x)
                    fprintf(fPCOEFF,['    ' num2str(str2num(xsplit{2})+OLDatomtypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  \n'])
                else
                    testflag2=1;
                end
            end
        end
        
        
        if startsWith(x,'Bond Coeffs') %Write base bonds coeffs
            x=fgetl(fid_in);
            x=fgetl(fid_in);
            xsplit=split(x);
            if exist('fBCOEFF','var')~=1
                fBCOEFF=fopen('BCOEFF','w');
                fprintf(fBCOEFF,'Bond Coeffs\n\n');
                fprintf(fBCOEFF,['    ' num2str(str2num(xsplit{2})+OLDbondtypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  \n'])
            else
                fprintf(fBCOEFF,['    ' num2str(str2num(xsplit{2})+OLDbondtypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  \n'])
            end
            testflag2=0;
            while testflag2==0
                x=fgetl(fid_in);
                xsplit=split(x);
                if ~isempty(x)
                    fprintf(fBCOEFF,['    ' num2str(str2num(xsplit{2})+OLDbondtypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  \n'])
                else
                    testflag2=1;
                end
            end
        end
        
        
        if startsWith(x,'Angle Coeffs') %Write base angle coeffs
            x=fgetl(fid_in);
            x=fgetl(fid_in);
            xsplit=split(x);
            if exist('fACOEFF','var')~=1
                fACOEFF=fopen('ACOEFF','w');
                fprintf(fACOEFF,'Angle Coeffs\n\n');
                fprintf(fACOEFF,['    ' num2str(str2num(xsplit{2})+OLDangletypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  \n'])
            else
                fprintf(fACOEFF,['    ' num2str(str2num(xsplit{2})+OLDangletypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  \n'])
            end
            testflag2=0;
            while testflag2==0
                x=fgetl(fid_in);
                xsplit=split(x);
                if ~isempty(x)
                    fprintf(fACOEFF,['    ' num2str(str2num(xsplit{2})+OLDangletypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  \n'])
                else
                    testflag2=1;
                end
            end
        end
        
        
        if startsWith(x,'Dihedral Coeffs') %Write base dihedral coeffs
            x=fgetl(fid_in);
            x=fgetl(fid_in);
            xsplit=split(x);
            if exist('fDCOEFF','var')~=1
                fDCOEFF=fopen('DCOEFF','w');
                fprintf(fDCOEFF,'Dihedral Coeffs\n\n');
                fprintf(fDCOEFF,['    ' num2str(str2num(xsplit{2})+OLDdihedraltypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) ' ' num2str(xsplit{5}) ' ' num2str(xsplit{6}) '  \n'])
            else
                fprintf(fDCOEFF,['    ' num2str(str2num(xsplit{2})+OLDdihedraltypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) ' ' num2str(xsplit{5}) ' ' num2str(xsplit{6}) '  \n'])
            end
            testflag2=0;
            while testflag2==0
                x=fgetl(fid_in);
                xsplit=split(x);
                if ~isempty(x)
                    fprintf(fDCOEFF,['    ' num2str(str2num(xsplit{2})+OLDdihedraltypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) ' ' num2str(xsplit{5}) ' ' num2str(xsplit{6}) '  \n'])
                else
                    testflag2=1;
                end
            end
        end
        
        
        if startsWith(x,'Improper Coeffs') %Write base improper coeffs
            x=fgetl(fid_in);
            x=fgetl(fid_in);
            xsplit=split(x);
            if exist('fICOEFF','var')~=1
                fICOEFF=fopen('ICOEFF','w');
                fprintf(fICOEFF,'Improper Coeffs\n\n');
                fprintf(fICOEFF,['    ' num2str(str2num(xsplit{2})+OLDimpropertypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  ' num2str(xsplit{5}) '  \n'])
            else
                fprintf(fICOEFF,['    ' num2str(str2num(xsplit{2})+OLDimpropertypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  ' num2str(xsplit{5}) '  \n'])
            end
            testflag2=0;
            while testflag2==0
                x=fgetl(fid_in);
                xsplit=split(x);
                if ~isempty(x)
                    fprintf(fICOEFF,['    ' num2str(str2num(xsplit{2})+OLDimpropertypes) ' ' num2str(xsplit{3}) '  ' num2str(xsplit{4}) '  ' num2str(xsplit{5}) '  \n'])
                else
                    testflag2=1;
                end
            end
        end
        
        
        
        
        
        
    end
end


TotalAtomTypes=OLDatomtypes;
TotalBondTypes=OLDbondtypes;
TotalAngleTypes=OLDangletypes;
TotalDihedralTypes=OLDdihedraltypes;
TotalImproperTypes=OLDimpropertypes;

atomtypes=0;
bondtypes=0;
angletypes=0;
dihedraltypes=0;
impropertypes=0;

OLDatomtypes=0;
OLDbondtypes=0;
OLDangletypes=0;
OLDdihedraltypes=0;
OLDimpropertypes=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% through
for i=1:numel(filelist)
    filein=filelist{i};
    if endsWith(filein,'.lmp')
        fileout=[filein(1:end-4) '_mixed.lmp'];
        %         fileouttemp=[filein(1:end-4) '_mixed.lmptemp'];
    elseif endsWith(filein,'.data')
        fileout=[filein(1:end-4) '_mixed.data'];
        %         fileouttemp=[filein(1:end-4) '_mixed.datatemp'];
    else
        error('files must end with .lmp or .data')
    end
    fid_in=fopen(filein,'r');
    fid_out=fopen(fileout,'w');
    fprintf(fid_out,['LAMMPS files created by LAMMPS Compiler - Peters research group \n\n'])
    %     fid_outtemp=fopen(fileouttemp,'w');
    
    x=fgetl(fid_in);
    %     fprintf(fid_outtemp,'temporary lammps file for polymatic mixing. Can delete after code is run\n')
    
    testflag=0;
    while testflag==0
        x=fgetl(fid_in);
        if isnumeric(x)
            if x==-1
                testflag=0;
                OLDatomtypes=OLDatomtypes+atomtypes;
                OLDbondtypes=OLDbondtypes+bondtypes;
                OLDangletypes=OLDangletypes+angletypes;
                OLDdihedraltypes=OLDdihedraltypes+dihedraltypes;
                OLDimpropertypes=OLDimpropertypes+impropertypes;
                break
            end
        end
        if isempty(x) || strcmp(x,' ')
            fprintf(fid_out,[x '\n']);
        end
        if endsWith(x,'atoms') || endsWith(x,'bonds') || endsWith(x,'angles') || endsWith(x,'dihedrals') || endsWith(x,'impropers')
            fprintf(fid_out,[x '\n']);
        end
        
        if endsWith(x,'atom types')
            xsplit=split(x);
            atomtypes=str2num(xsplit{2});
            fprintf(fid_out,['    ' num2str(TotalAtomTypes) ' atom types\n']);
            
        end
        if endsWith(x,'bond types')
            xsplit=split(x);
            bondtypes=str2num(xsplit{2});
            fprintf(fid_out,['    ' num2str(TotalBondTypes) ' bond types\n']);
        end
        if endsWith(x,'angle types')
            xsplit=split(x);
            angletypes=str2num(xsplit{2});
            fprintf(fid_out,['    ' num2str(TotalAngleTypes) ' angle types\n']);
        end
        if endsWith(x,'dihedral types')
            xsplit=split(x);
            dihedraltypes=str2num(xsplit{2});
            fprintf(fid_out,['    ' num2str(TotalDihedralTypes) ' dihedral types\n']);
        end
        if endsWith(x,'improper types')
            xsplit=split(x);
            impropertypes=str2num(xsplit{2});
            fprintf(fid_out,['    ' num2str(TotalImproperTypes) ' improper types\n']);
        end
        
        if endsWith(x,'hi')
            fprintf(fid_out,[x '\n']);
        end
        
        if startsWith(x,'Masses')
            fMASSES=fopen('MASSES','r');
            testMASS=0;
            while testMASS==0
                x=fgetl(fMASSES);
                if isnumeric(x)
                    if x==-1
                        break
                    end
                end
                fprintf(fid_out,[x '\n']);
                
            end
            x='0';
        end
        
        if startsWith(x,'Pair Coeffs')
            fPCOEFF=fopen('PCOEFF','r');
            testPCOEFF=0;
            while testPCOEFF==0
                x=fgetl(fPCOEFF);
                if isnumeric(x)
                    if x==-1
                        break
                    end
                end
                fprintf(fid_out,[x '\n']);
                
            end
            x='0';
        end
        
        if startsWith(x,'Bond Coeffs')
            fBCOEFF=fopen('BCOEFF','r');
            testBCOEFF=0;
            while testBCOEFF==0
                x=fgetl(fBCOEFF);
                if isnumeric(x)
                    if x==-1
                        break
                    end
                end
                fprintf(fid_out,[x '\n']);
                
            end
            x='0';
        end
        
        if startsWith(x,'Angle Coeffs')
            fACOEFF=fopen('ACOEFF','r');
            testACOEFF=0;
            while testACOEFF==0
                x=fgetl(fACOEFF);
                if isnumeric(x)
                    if x==-1
                        break
                    end
                end
                fprintf(fid_out,[x '\n']);
                
            end
            x='0';
        end
        
        if startsWith(x,'Dihedral Coeffs')
            fDCOEFF=fopen('DCOEFF','r');
            testDCOEFF=0;
            while testDCOEFF==0
                x=fgetl(fDCOEFF);
                if isnumeric(x)
                    if x==-1
                        break
                    end
                end
                fprintf(fid_out,[x '\n']);
                
            end
            x='0';
        end
        
        if startsWith(x,'Improper Coeffs')
            fICOEFF=fopen('ICOEFF','r');
            testICOEFF=0;
            while testICOEFF==0
                x=fgetl(fICOEFF);
                if isnumeric(x)
                    if x==-1
                        break
                    end
                end
                fprintf(fid_out,[x '\n']);
                
            end
            x='0';
        end
        
        if startsWith(x,'Atoms')
            fprintf(fid_out,'Atoms \n\n');
            x=fgetl(fid_in);
            testATOMS=0;
            while testATOMS==0
                x=fgetl(fid_in);
                if isempty(x)
                    break
                end
                xsplit=split(x);
                if length(xsplit{1})==0
                    fprintf(fid_out,['    ' xsplit{2} '    ' xsplit{3} '    ' num2str(str2num(xsplit{4})+OLDatomtypes) '   ' xsplit{5} '   ' xsplit{6} '   ' xsplit{7} '   ' xsplit{8} '\n'])
                elseif length(xsplit{1})>0
                    fprintf(fid_out,['    ' xsplit{1} '    ' xsplit{2} '    ' num2str(str2num(xsplit{3})+OLDatomtypes) '   ' xsplit{4} '   ' xsplit{5} '   ' xsplit{6} '   ' xsplit{7} '\n'])
                else
                    error('atom line read not handled correctly')
                end
            end
            fprintf(fid_out,'\n');
            x='0';
        end
        
        
        if startsWith(x,'Bonds')
            fprintf(fid_out,'Bonds \n\n');
            x=fgetl(fid_in);
            testATOMS=0;
            while testATOMS==0
                x=fgetl(fid_in);
                if isempty(x)
                    break
                end
                xsplit=split(x);
                fprintf(fid_out,['    ' xsplit{2} '    ' num2str(str2num(xsplit{3})+OLDbondtypes) '   ' num2str(str2num(xsplit{4})) '   ' num2str(str2num(xsplit{5})) '\n'])
            end
            fprintf(fid_out,'\n');
            x='0';
        end
            
        if startsWith(x,'Angles')
            fprintf(fid_out,'Angles \n\n');
            x=fgetl(fid_in);
            testATOMS=0;
            while testATOMS==0
                x=fgetl(fid_in);
                if isempty(x)
                    break
                end
                xsplit=split(x);
                fprintf(fid_out,['    ' xsplit{2} '    ' num2str(str2num(xsplit{3})+OLDangletypes) '   ' num2str(str2num(xsplit{4})) '   ' num2str(str2num(xsplit{5})) '   ' num2str(str2num(xsplit{6})) '\n'])
            end
            fprintf(fid_out,'\n');
            x='0';
        end
        
        if startsWith(x,'Dihedrals')
            fprintf(fid_out,'Dihedrals \n\n');
            x=fgetl(fid_in);
            testATOMS=0;
            while testATOMS==0
                x=fgetl(fid_in);
                if isempty(x)
                    break
                end
                xsplit=split(x);
                fprintf(fid_out,['    ' xsplit{2} '    ' num2str(str2num(xsplit{3})+OLDdihedraltypes) '   ' num2str(str2num(xsplit{4})) '   ' num2str(str2num(xsplit{5})) '   ' num2str(str2num(xsplit{6})) '   ' num2str(str2num(xsplit{7})) '\n'])
            end
            fprintf(fid_out,'\n');
            x='0';
        end
        
        if startsWith(x,'Impropers')
            fprintf(fid_out,'Impropers \n\n');
            x=fgetl(fid_in);
            testATOMS=0;
            while testATOMS==0
                x=fgetl(fid_in);
                if isempty(x)
                    break
                end
                if isnumeric(x)
                    if x==-1
                        break
                    end
                end
                xsplit=split(x);
                fprintf(fid_out,['    ' xsplit{2} '    ' num2str(str2num(xsplit{3})+OLDimpropertypes) '   ' num2str(str2num(xsplit{4})) '   ' num2str(str2num(xsplit{5})) '   ' num2str(str2num(xsplit{6})) '   ' num2str(str2num(xsplit{7})) '\n'])
            end
            fprintf(fid_out,'\n');
            x='0';
        end
        
        %Wtite rest of files, MASSES, COEFFS, atoms, bonds, etc.
        
        
        
        
        
    end
end
