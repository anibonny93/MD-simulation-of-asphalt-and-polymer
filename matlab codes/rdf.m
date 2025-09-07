function [r,gr,numtimesteps] = rdf({'varargin)
%Calculates the radial pair distibution function using vmd implementation
 

p=inputParser;
defaultInputnames='test.gsd';
defaultSelection1='name A';
defaultSelection2='name A';
defaultGPU=1;
defaultDelta=0.1;
defaultRmax=30;
defaultUsepbc=1;
defaultSelupdate=0;
defaultFirstStep=0;
defaultLastStep=-1;
defaultStepSize=1;
defaultStepSizeTime=0;
defaultOutFile='gr.dat';

lammpstestflag=0;

addParameter(p,'Inputnames',defaultInputnames)
addParameter(p,'Selection1',defaultSelection1)PP_oxi2023_1.
addParameter(p,'Selection2',defaultSelection2)
addParameter(p,'GPU',defaultGPU,@isnumeric)
addParameter(p,'Delta',defaultDelta,@isnumeric)
addParameter(p,'Rmax',defaultRmax,@isnumeric)
addParameter(p,'Usepbc',defaultUsepbc,@isnumeric)
addParameter(p,'Selupdate',defaultSelupdate,@isnumeric)
addParameter(p,'FirstStep',defaultFirstStep,@isnumeric)
addParameter(p,'LastStep',defaultLastStep,@isnumeric)
addParameter(p,'StepSize',defaultStepSize,@isnumeric)
addParameter(p,'StepSizeTime',defaultStepSizeTime,@isnumeric)
addParameter(p,'OutFile',defaultOutFile)

parse(p,varargin{:});
fields=fieldnames(p.Results);
for i=1:numel(fields)
    eval([fields{i} '=' 'p.Results.' fields{i} ';']);
end

%open tcl file to write
tclscript=fopen('rdf.tcl','w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Writing
%write out lines to load mols
if class(Inputnames)=='char'
    if Inputnames(end-2:end)=='gsd'
        fprintf(tclscript,['mol new "' Inputnames '" waitfor -1\n']);
    elseif Inputnames(end-2:end)=='xml'
        fprintf(tclscript,['mol new "' Inputnames '" type {hoomd} waitfor -1\n']);
    else
        error('Invalid File Type');
    end
elseif class(Inputnames)=='cell'
    if Inputnames{1}(end-2:end)=='gsd'
        fprintf(tclscript,['mol new "' Inputnames{1} '"\n']);
    elseif Inputnames{1}(end-2:end)=='xml'
        fprintf(tclscript,['mol new "' Inputnames{1} '" type {hoomd}\n']);
    elseif Inputnames{1}(end-3:end)=='data' %Lammps data and dump
        %check to see if there are exactly two inputs
        if size(Inputnames)~=[1 2];
            error('size of lammps input names is not 2 cells')
        end
        lammpstestflag=1;
        fprintf(tclscript,['package require topotools \n']);
        fprintf(tclscript,['topo readlammpsdata ' Inputnames{1} '\n']);
        fprintf(tclscript,['mol addfile {' Inputnames{2} '} type {lammpstrj} first 0 last -1 step 1 waitfor -1 0\n']);

    else
        error('Invalid File Type');
    end
    if lammpstestflag ~=1
        for i=2:length(Inputnames)
            fprintf(tclscript,['mol addfile "' Inputnames{i} '" waitfor all\n']);
        end
    end
end

% set up atom selections
if iscell(Selection1)
    %check that cell arrays are set up correctly
    if length(Selection1)~=length(Selection2)
        error('Selection cell arrays not the same size');
    end
    if iscell(OutFile)
        if length(OutFile)~=length(Selection1)
            error('Selection cell array size not equal to OutFile cell array size');
        end
    elseif ~strcmp(OutFile,defaultOutFile)
        error('OutFile not cell array when Selection is cell array');
    end
    %writeout set selects to script
    for i=1:length(Selection1)
        fprintf(tclscript,['set sel1 [atomselect top "' Selection1{i} '"]\n']);
        fprintf(tclscript,['set sel2 [atomselect top "' Selection2{i} '"]\n']);
        [maxi,numtimesteps]=calcrg(StepSizeTime,tclscript,Delta,Rmax,Usepbc,Selupdate,FirstStep,LastStep,StepSize,OutFile{i});
    end
else
    fprintf(tclscript,['set sel1 [atomselect top "' Selection1 '"]\n']);
    fprintf(tclscript,['set sel2 [atomselect top "' Selection2 '"]\n']);
    [maxi,numtimesteps]=calcrg(StepSizeTime,tclscript,Delta,Rmax,Usepbc,Selupdate,FirstStep,LastStep,StepSize,OutFile);
end
fprintf(tclscript,['exit\n']);
disp('Running in VMD');
if isunix
    system('nice -n 1 vmd -dispdev text -eofexit <rdf.tcl> output.log')
else
    system('vmd -dispdev text -e rdf.tcl > output.log')
end

%Collect data from written files
if iscell(OutFile)
    for j=1:length(OutFile)
        if StepSizeTime==0
            file_A_A = fopen([OutFile{j}]);
            A_A = textscan(file_A_A, '%f %f');
            fclose(file_A_A);
            delete([OutFile{j}]);
            r_A_A = A_A{1};
            gr_A_A = A_A{2};
            indx_A_A = find(gr_A_A, 1, 'last');
            r = r_A_A(1:indx_A_A);
            gr = gr_A_A(1:indx_A_A);
        else
            for i=1:maxi
                OutFile2=[OutFile{j} num2str(i)];
                file_A_A = fopen([OutFile2]);
                A_A = textscan(file_A_A, '%f %f');
                fclose(file_A_A);
                delete(OutFile2);
                r_A_A = A_A{1};
                gr_A_A = A_A{2};
                indx_A_A = find(gr_A_A, 1, 'last');
                if i==1 && j==1
                    r = r_A_A(1:indx_A_A);
                    gr = gr_A_A(1:indx_A_A);
                    rsize=length(r);
                else
                    r = [r, r_A_A(1:size(r,1))];
                    gr = [gr, gr_A_A(1:size(r,1))];

                end
            end
        end
    end
else
    if StepSizeTime==0
        file_A_A = fopen([OutFile]);
        A_A = textscan(file_A_A, '%f %f');
        fclose(file_A_A);
        delete([OutFile]);
        r_A_A = A_A{1};
        gr_A_A = A_A{2};
        indx_A_A = find(gr_A_A, 1, 'last');
        r = r_A_A(1:indx_A_A);
        gr = gr_A_A(1:indx_A_A);
    else
        for i=1:maxi
            OutFile2=[OutFile num2str(i)];
            file_A_A = fopen([OutFile2]);
            A_A = textscan(file_A_A, '%f %f');
            fclose(file_A_A);
            delete([OutFile2]);
            r_A_A = A_A{1};
            gr_A_A = A_A{2};
            indx_A_A = find(gr_A_A, 1, 'last');
            if i==1
                r = r_A_A(1:indx_A_A);
                gr = gr_A_A(1:indx_A_A);
                rsize=length(r);
            else
                r = [r, r_A_A(1:indx_A_A)];
                gr = [gr, gr_A_A(1:indx_A_A)];
            end
        end
    end
end
end

function [maxi,numtimesteps]=calcrg(StepSizeTime,tclscript,Delta,Rmax,Usepbc,Selupdate,FirstStep,LastStep,StepSize,OutFile)
%calculateg(r) Write the rest of the current script to calculate rg and put
%in the outfile
if StepSizeTime==0
    fprintf(tclscript,['set gr [measure rdf $sel1 $sel2 delta ' num2str(Delta) ' rmax ' num2str(Rmax) ' usepbc ' num2str(Usepbc) ' selupdate ' num2str(Selupdate) ' first ' num2str(FirstStep) ' last ' num2str(LastStep) ' step ' num2str(StepSize) ']\n']);
    %setup up outfile
    fprintf(tclscript,['set outfile [open "' OutFile '" w]\n']);
    fprintf(tclscript,['set r [lindex $gr 0]\n']);
    fprintf(tclscript,['set gr2 [lindex $gr 1]\n']);
    fprintf(tclscript,['foreach p $r x $gr2 {puts $outfile "$p $x"}\n']);
    fprintf(tclscript,['close $outfile\n']);
    numtimesteps=1;
    maxi=1;
%     fprintf(tclscript,['exit\n']);
else
    FStep=FirstStep;
    LStep=FStep+StepSizeTime*StepSize-1;
    for i=1:((LastStep-FirstStep)/StepSize+1)/StepSizeTime
        fprintf(tclscript,['set gr' num2str(i) ' [measure rdf $sel1 $sel2 delta ' num2str(Delta) ' rmax ' num2str(Rmax) ' usepbc ' num2str(Usepbc) ' selupdate ' num2str(Selupdate) ' first ' num2str(FStep) ' last ' num2str(LStep) ' step ' num2str(StepSize) ']\n']);
        %setup up outfile
        OutFile2=[OutFile num2str(i)];
        fprintf(tclscript,['set outfile [open "' OutFile2 '" w]\n']);
        fprintf(tclscript,['set r [lindex $gr' num2str(i) ' 0]\n']);
        fprintf(tclscript,['set gr2 [lindex $gr' num2str(i) ' 1]\n']);
        fprintf(tclscript,['foreach p $r x $gr2 {puts $outfile "$p $x"}\n']);
        fprintf(tclscript,['close $outfile\n']);
        maxi=i;
        FStep=LStep+1;
        LStep=FStep+StepSizeTime*StepSize-1;
    end
    numtimesteps=floor(((LastStep-FirstStep)/StepSize+1)/StepSizeTime);
%     fprintf(tclscript,['exit\n']);
end
end
