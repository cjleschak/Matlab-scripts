%Used to move downloaded behavioral data into the behav folders of subjects.

%=====================================================================%
clear all; home;
%spm('defaults','fmri');     % initializes SPM defaults for fMRI modality
%spm_jobman('initcfg');      % initializes job configurations


%=====================================================================%
% USER DEFINED
%---------------------------------------------------------------------%
%executeTAG=1;           % execute job(s) immediately? 0=no, 1=yes

studyDIR='/u/project/sanscn/data/SAFETY'; cd(studyDIR);
subID='SAS';
behavPath='/u/home/c/cjlescha/Downloads/Faces/Faces/Faces_7pt/data';

%runIDs={'BOLD_Faces1','BOLD_Faces2','BOLD_Faces3'};

taskID='Faces';
%analysisID='Faces_hpf128ar1_3runs';


% Optional inputs - leave empty if you don't want to use them
skipsub = {};               % skip these subjects, use the full subject's folder name in single quotes (i.e. {'w001' 'w002'})
subnam = {};                % do only these subjects, format is the same as skipsup (i.e. {'w001' 'w002'})


%---------------------------------------------------------------------%
d=dir([subID '*']);
if exist('subnam','var')
    if isempty(subnam)
        d=dir([subID '*']);
        for i=1:length(d)
            subnam{i}=d(i).name;
            fprintf('Adding %s\n',subnam{i})
        end
    end
else
    d=dir([subid '*']);
    for i=1:length(d)
        subnam{i}=d(i).name;
        fprintf('Adding %s\n',subnam{i})
    end
end

subnum = length(subnam);
dosubs = 1:subnum;

if size(skipsub,2) > size(skipsub,1) % keep skipsub vector vertical
    skipsub = skipsub';
end
skipped_log = {};

%%---------%%

for s = dosubs
    
    cbusub = sprintf('%s',subnam{s});
    fprintf('\nStarting subject %d of %d, %s\n',s,length(dosubs),cbusub)
    swd = fullfile(studyDIR,cbusub);  % subject's working directory
    fprintf('Subject directory is %s\n',swd)
    cd(swd)
    cd behav
    base_dir = pwd;     % subject's base behavioral data directory (behav folder)
    
    % Check subject skip list
    %-----------------------------------------------------------------%
    if sum(strcmpi(skipsub,subnam{s})) ~= 0
        fprintf('Subject is on list of subjects to skip.  Skipping...\n')
        skipped_log = [skipped_log;cbusub];
        continue;
    end;
    
    % Move behav data to sub's behav folder
    behavPath4sub=fullfile(behavPath,cbusub);
    if ~exist(behavPath4sub)
        fprintf('No 7 pt behav foldesr in downloads folder for %s.  Checking 4 pt data...\n',cbusub)
        behavPath4sub=fullfile('/u/home/c/cjlescha/Downloads/Faces/Faces/Faces_4pt/data',cbusub);
        if ~exist(behavPath4sub)
            fprintf('No behav folder in downloads folder.  Skipping %s...\n',cbusub)
            continue
        end
    end
    movefile(behavPath4sub,base_dir)
    cd(base_dir)
    movefile(cbusub,'Faces')
    
end


