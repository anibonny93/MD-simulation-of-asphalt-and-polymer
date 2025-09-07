function [r,gr,q,Sq]=SingleSq({'AP_PP_oxi2023_1_replicated.data','AP_PP_oxi2023_1.dump'}),('Selection1','resid/2.0 != floor(resid/2.0)','Selection2','resid/2.0 = floor(resid/2.0)'))

%Calculates S(q) for single file and plots. Possible parameters are
%FirstStep, Last Step, StepSize, Selections1, and Selections2. Will set
%Selections2 to Selections1 unless otherwise specified. Defaults to all
%timesteps. 








p=inputParser;
defaultfilename='test.gsd';
defaultsinglepath=0;
defaultplotGR=1;
defaultsaveGR=1;
defaultplotSQ=1;
defaultsaveSQ=1;
defaultlabelqstar=0;
defaultplotQmax=1;
defaultsaveQmax=1;
defaultplotD=1;
defaultsaveD=1;
defaultcloseplots=1;
defaultSelections1='name A';
defaultSelections2='default';
defaultStepSizeTime=1;
defaultcustomaxis=[];
defaultXmlStartidIntifierLength=9;
defaultDelta=0.1;
defaultFirstStep=0; %if <=0 then relative to last step
defaultLastStep=-1; %if <=0 then relative to last step
defaultStepSize=1;
defaultOverwriteFiles=0;
defaultShowPlots=0;

% path='/mnt/Shared_Data/Nanoparticles/PCND_Nano_SLJDPD/Refine_PCND_Nano/TestXi0result/testbestsetsagain/';
% filename='PCND_Nano_XiNano0.02_tNano100000_XiPoly0_tPoly100000_longhalf_start.xml';
singlepath=0;%0=singlefile 1=singlepath 2=allpaths
plotGR=1;
saveGR=1;
plotSQ=1;
saveSQ=1;
labelqstar=0;
plotQmax=1;
saveQmax=1;
plotD=1;
saveD=1;
closeplots=1;
Selections1={'name D'};
Selections2=Selections1;
StepSizeTime=1;
customaxis=[];
XmlStartidIntifierLength=9;
Delta=0.1;
FirstStep=1; %if <=0 then relative to last step
LastStep=-1; %if <=0 then relative to last step
StepSize=1;
OverwriteFiles=0;
ShowPlots=0;

% addParameter(p,'filename',defaultfilename)
% addParameter(p,'singlepath',defaultsinglepath,@isnumeric)
% addParameter(p,'plotGR',defaultplotGR,@isnumeric)
% addParameter(p,'saveGR',defaultsaveGR,@isnumeric)
% addParameter(p,'plotSQ',defaultplotSQ,@isnumeric)
% addParameter(p,'saveSQ',defaultsaveSQ,@isnumeric)
% addParameter(p,'labelqstar',defaultlabelqstar,@isnumeric)
% addParameter(p,'plotQmax',defaultplotQmax,@isnumeric)
% addParameter(p,'saveQmax',defaultsaveQmax,@isnumeric)
% addParameter(p,'plotD',defaultplotD,@isnumeric)
% addParameter(p,'saveD',defaultsaveD,@isnumeric)
% addParameter(p,'closeplots',defaultcloseplots,@isnumeric)
addParameter(p,'Selections1',defaultSelections1)
addParameter(p,'Selections2',defaultSelections2)
% addParameter(p,'StepSizeTime',defaultStepSizeTime,@isnumeric)
% addParameter(p,'customaxis',defaultcustomaxis,@isnumeric)
% addParameter(p,'XmlStartidIntifierLength',defaultXmlStartidIntifierLength,@isnumeric)
% addParameter(p,'Delta',defaultDelta,@isnumeric)
addParameter(p,'FirstStep',defaultFirstStep,@isnumeric)
addParameter(p,'LastStep',defaultLastStep,@isnumeric)
addParameter(p,'StepSize',defaultStepSize,@isnumeric)
% addParameter(p,'OverwriteFiles',defaultOverwriteFiles,@isnumeric)
% addParameter(p,'ShowPlots',defaultShowPlots,@isnumeric)

parse(p,varargin{:});
fields=fieldnames(p.Results);
for i=1:numel(fields)
    eval([fields{i} '=' 'p.Results.' fields{i} ';']);
    if strcmp(fields{i},'Selections2')
        if strcmp(Selections2,defaultSelections2)
            Selections2=Selections1;
        end
    end
end



% filename='PCND_Nano_XiNano0.02_tNano100000_XiPoly0_tPoly100000_longhalf_start.xml';
% singlepath=0;%0=singlefile 1=singlepath 2=allpaths
% plotGR=1;
% saveGR=1;
% plotSQ=1;
% saveSQ=1;
% labelqstar=0;
% plotQmax=1;
% saveQmax=1;
% plotD=1;
% saveD=1;
% closeplots=1;
% Selections1={'name D','name A'};
% Selections2=Selections1;
% StepSizeTime=1;
% customaxis=[];
% XmlStartidIntifierLength=9;
% Delta=0.1;
% FirstStep=1; %if <=0 then relative to last step
% LastStep=0; %if <=0 then relative to last step
% StepSize=1;
% OverwriteFiles=1;


%% Find files and run vmd
if length(Selections1)~=length(Selections2)
    error('Selections not of equal size');
end

close all

[r,gr] = rdf('StepSizeTime',0,'Inputnames',filename,'Selection1',Selections1,'Selection2',Selections2,'Delta',Delta,'Rmax',100,'FirstStep',FirstStep,'LastStep',LastStep,'StepSize',StepSize,'OutFile','temp.csv');              

[q,Sq]=localplotting(r,gr)

end

function [q,Sq]=localplotting(r,gr)%,numtimesteps,Outfiles,xmlname,dcdname,nsteps)
%% bring in all variables
ALLVAR = evalin('caller','whos');
for ii = 1:length(ALLVAR)
   C_ =  evalin('caller',[ALLVAR(ii).name ';']);
   eval([ALLVAR(ii).name,'=C_;']);
end
clear ALLVAR C_ ii
Outfiles=1;
numtimesteps=1;

%% Plotting of GR

if plotGR==1
    for i=1:1
        for k=1:1
            if k==1
                figure
                plot(r, gr, 'LineWidth', 2);
                if iscell(Selections1)
                    title({'Radial Pair Distribution Function g(r)', [Selections1{i} ':' Selections2{i}]});
                else
                    title({'Radial Pair Distribution Function g(r)', [Selections1 ':' Selections2]});
                end
                axis([0 max(r) 0 max(gr)+0.5]);
                xlabel('r');
                ylabel('g(r)');
                grid on;
                hold on
            end
%             counter=counter+1;
        end
    end
end

%% Calculation of S(q)  and q* and D          
counter=1;
if plotSQ==1;
    for i=1:1
        for k=1:1
            %%%%%%%%%%%%%%%%%
            rAA=r;
            grAA=gr;
            rmax_AA = rAA(end);
            ru0_AA = 1;
            n1 = 1;
            for m1 = 1:length(rAA)
                if m1 == 1
                    del_r_AA(1,1) = rAA(1,1) - 0;
                else
                    del_r_AA(m1,1) = rAA(m1,1) - rAA(m1-1,1);
                end
            end
            qs_AA = 0.001:0.001:10;
            for q1 = qs_AA
                S_AA = 4.*pi.*ru0_AA.*rAA.*((sin(q1.*rAA))./q1).*(grAA-1).*(sin((pi.*rAA)./rmax_AA)./((pi.*rAA)/rmax_AA));
                Sq_AA(n1,1) = (sum(S_AA(1:end).*del_r_AA))-1;
                n1 = n1 + 1;
            end
            Sq(:,counter)=Sq_AA;
            q(:,counter)=qs_AA;
            z_AA = find(Sq_AA == max(Sq_AA));
            q_star_AA = qs_AA(z_AA);
            Lm_AA=2*pi()/q_star_AA;
            E_AA=0.142*(Lm_AA/rmax_AA)^(1.92);
            D_AA=Lm_AA/(E_AA+1);
            qstar(counter)=q_star_AA;
            D(counter)=D_AA;
            qmax(counter)=max(Sq_AA);
            %%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%plot
            if k==1
                figure
                plot(q(:,counter), Sq(:,counter), 'LineWidth', 2);
                if labelqstar==1
                    text(q_star_AA, Sq_AA(z_AA), ['\leftarrow q* = ' num2str(q_star_AA) '    D = ' num2str(D_AA)]);
                end
                if iscell(Selections1)
                    title({'Structure Factor S(q)', [Selections1{i} ':' Selections2{i}]});
                else
                    title({'Structure Factor S(q)', [Selections1 ':' Selections2]});
                end
                xlabel('q');
                ylabel('S(q)');
                if customaxis~=[]
                    axis(customaxis);
                end
                grid on;
                hold on
            else
                plot(q(:,counter), Sq(:,counter), 'LineWidth', 2);
                if labelqstar==1
                    text(q_star_AA, Sq_AA(z_AA), ['\leftarrow q* = ' num2str(q_star_AA) '    D = ' num2str(D_AA)]);
                end
            end
            counter=counter+1;
        end
    end
end





fclose all
end
