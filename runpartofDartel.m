
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXECUTE PREDARTEL AND/OR DARTEL MATLABBATCHES FROM WRAPPER SCRIPT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Directories where batch files are saved.
batch_dir = '/u/project/sanscn/data/SAFETY/_testBatch/_batches/PreProcessing/DARTEL';
preDARbatch_dir = ''; %pre-dartel batch folder (within batch_dir)
DARbatch_dir    = [batch_dir '/DARTEL_20180923_0423']; %dartel batch folder (within batch_dir)


% What should be executed?
execPREDARTEL = 0;      %run Pre-Dartel? 1=execute; 0=no/ignore
execDARTEL = 1;         %run DARTEL?     1=execute; 0=no/ignore


%%%%%%%%%%%%%%%%%%%%%%%%%%% End User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% File names
preDARTELfile = 'preDARTEL*.mat';
workspaceDARTELfile = 'forDartel_workspace*.mat';
DARTELfile = 'DARTEL*.mat';

cd(batch_dir);

% Find files.
if execPREDARTEL==1
    preDARbatch_dirFull=[batch_dir preDARbatch_dir];
    if ~strcmp(preDARbatch_dirFull,batch_dir)
        cd(preDARbatch_dirFull); preDARTELfiles = dir(preDARTELfile); cd ..
    end
%    cd(preDARbatch_dir); preDARTELfiles = dir(preDARTELfile); cd ..
end
if execDARTEL==1
    DARbatch_dirFull=[batch_dir DARbatch_dir];
    if ~strcmp(DARbatch_dirFull,batch_dir)
        cd(DARbatch_dirFull); workspaceDARTELfiles = dir(workspaceDARTELfile);
        DARTELfiles = dir(DARTELfile); cd ..
    end
     %cd(DARbatch_dir); workspaceDARTELfiles = dir(workspaceDARTELfile);
   % DARTELfiles = dir(DARTELfile); cd ..
end

%Check that necessary files are there
if execPREDARTEL==1
    if isempty(preDARTELfiles)
        sprintf('\nERROR! \nNo preDARTEL*.mat file(s) found. \nUpdate ''preDARTELfile''.')
        return
    end
end
if execDARTEL==1
    if length(workspaceDARTELfiles)>1
        sprintf('\nERROR! \nMore than one workspaceDARTEL*.mat file found. \nUpdate ''workspaceDARTELfile''.')
        return
    elseif isempty(workspaceDARTELfiles)
        sprintf('\nERROR! \nNo workspaceDARTEL*.mat file found. \nUpdate ''workspaceDARTELfile''.')
        return
    end
    if length(DARTELfiles)>1
        sprintf('\nERROR! \nMore than one DARTEL*.mat file found. \nUpdate ''DARTELfile''.')
        return
    elseif isempty(DARTELfiles)
        sprintf('\nERROR! \nNo DARTEL*.mat file found. \nUpdate ''DARTELfile''.')
        return
    end
end
  
%Execute Pre-DARTEL
if execPREDARTEL==1
    for j=1:length(preDARTELfiles)
        cd(preDARbatch_dir);
        load(preDARTELfiles(j).name)
        spm_jobman('run',matlabbatch)
        clear matlabbatch
    end
end
        
%Execute DARTEL
if execDARTEL==1
    cd(DARbatch_dir);
    load(workspaceDARTELfiles.name)
    load(DARTELfiles.name)
    spm_jobman('run',matlabbatch)
    clear matlabbatch
end

