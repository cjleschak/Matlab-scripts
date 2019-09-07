%=====================================================================%
%   Move Contrasts into Common Location 
%
%   This script loops through subject folders and finds contrast files
%   to move them into a common location for second level contrast.
%
%   Created 4/25/2014 by Jared Torre
%
%=====================================================================%
clear all; home;

% User Defined 
%---------------------------------------------------------------------%

% Set up directories where subject folders are located:
owd = '/u/project/sanscn/data/SAFETY'; % main study directory, no "/" at end
cd(owd);

% Subject ID: give the first starting characters:
subID = 'SAS_';

% Task directory: identify the name of the second level task directory:
task_dir = 'HAE/Habit_Acqui';

% First level directory: identify the analysis folder you want to use at the second level:
analysis_dir = 'Habit_Aqui_firstlast4_hpf128ar1';

% Include any subject folders you want skipped here (i.e. skipsub = % {'NMPRE_602' 'NMPOST_123'};):
skipsub = {};

% If you only want to run specific subjects, add them here
subnam = {};


%------------------------END USER INPUT-------------------------------%

if size(skipsub,2) > size(skipsub,1) % keep skipsub vector vertical
    skipsub = skipsub';
end


% Make directories (if needed)
%---------------------------------------------------------------------%
gwd = fullfile(owd,'_contrasts',task_dir,analysis_dir);
eval(sprintf('mkdir %s;',gwd));

% Find subject directories
%---------------------------------------------------------------------%
d=dir([subID '*']);
if exist('subnam','var')
    if isempty(subnam)
        d=dir([subID '*']);
        for i=1:length(d)
            subnam{i,1}=d(i).name;
            full_subnam{i,1}=fullfile(owd,subnam{i,1});
            fprintf('Adding %s\n',subnam{i})
        end
    end
else
    d=dir([subID '*']);
    for i=1:length(d)
        subnam{i,1}=d(i).name;
        full_subnam{i,1}=fullfile(owd,subnam{i,1});
        fprintf('Adding %s\n',subnam{i})
    end
end

subnum = length(subnam);
dosubs = 1:subnum;

% Populate Contrast List
%---------------------------------------------------------------------%
conCount = 1;
contrastInfo.name = {};
contrastInfo.file = {};
contrastInfo.sub = {};
fprintf('\nPopulating contrast lists...')
for s = dosubs
    cbusub = subnam{s};
    full_cbusub = full_subnam{s};
    sub_awd = fullfile(full_subnam{s},'analysis',task_dir,analysis_dir);
    if sum(ismember(skipsub,cbusub))>0
        fprintf('Subject %s is in exclusion list.  Skipping...\n',subnam{s})
        continue;
    end
    
    if ~exist(sub_awd,'dir')
        fprintf('Analysis directory does not exist for subject %s.  Skipping...\n',subnam{s})
        skipsub = [skipsub; cbusub];
        continue;
    end
    
    cd(sub_awd)
    d = dir('SPM.mat');
    if ~isempty(d)
        load SPM
    else
        fprintf('No SPM.mat file in subject %s analysis directory.  Skipping...\n',subnam{s})
        skipsub = [skipsub; cbusub];
        continue;
    end
        
    % add current subjects contrasts to giant list of all contrasts
    for i = 1:length(SPM.xCon)
        imageFilename = SPM.xCon(i).Vcon.fname;
        fullImageFilename = fullfile(sub_awd,imageFilename);
        contrastInfo.name{conCount,1} = SPM.xCon(i).name;
        contrastInfo.file{conCount,1} = fullImageFilename;
        contrastInfo.sub{conCount,1} = cbusub;
        conCount = conCount +1;
    end
clear SPM file_array
end
fprintf('Done.\n')

% Remove ' - All Sessions' suffix if it exists
contrastInfo.name_nosuffix = regexprep(contrastInfo.name,' - All Sessions','');

fprintf('Copying files...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPY LOOP                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for c = 1:length(contrastInfo.name_nosuffix)
    current_conname=contrastInfo.name_nosuffix{c};
    current_conname_nospaces = regexprep(current_conname,' ','_');
    new_filename = [current_conname_nospaces '_' contrastInfo.sub{c}];
    
    current_fullfile_noext = strsplit(contrastInfo.file{c},'.');
    if strcmpi(current_fullfile_noext{2},'img')
        new_fullimgname = fullfile(gwd,[new_filename '.img']);
        copyfile([current_fullfile_noext{1} '.img'],new_fullimgname);
        new_fullhdrname = fullfile(gwd,[new_filename '.hdr']);
        copyfile([current_fullfile_noext{1} '.hdr'],new_fullhdrname);
    else
        new_fullfilename = fullfile(gwd,[new_filename '.nii']);
        copyfile(contrastInfo.file{c},new_fullfilename);
    end
end
fprintf('Done\n')
skipsub