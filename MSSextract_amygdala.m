
%This script uses the .mat files created from the amygdala localizer task
%(run using the MSS script), and constructs names, onsets, and duration
%variables, and saves these for each subject in a .mat file. These subject-
%specific .mat files can then be read directly into SPM for first level 
%model specification.

% July 19, 2017 - modified CJL

%=====================================================================%
% USER DEFINED
%---------------------------------------------------------------------%
clear all;

studyDIR='/u/project/sanscn/data/SSBCS'; cd(studyDIR);
subID='SSPB';

% Optional inputs - leave empty if you don't want to use them
skipsub = {};               % skip these subjects, use the full subject's folder name in single quotes (i.e. {'w001' 'w002'})
subnam = {};
%subnam = {'SSPB_59','SSPB_79'};                % do only these subjects, format is the same as skipsup (i.e. {'w001' 'w002'})


%=====================================================================%
% Find subject directories
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
skipped_log = [];
reason = {};

%=====================================================================%
% Start Subject Loop
%---------------------------------------------------------------------%
for s = dosubs
    
    cbusub = sprintf('%s',subnam{s});
    fprintf('Starting subject %d of %d, %s\n',s,length(dosubs),cbusub)
    swd = strcat(studyDIR,filesep,cbusub,filesep,'behav');  % subject's working directory
    fprintf('Subject''s behavioral directory is %s\n',swd)
    mkdir(swd,'extracted');
    cd(swd)
    base_dir = pwd;     % subject's behav directory
    
    % Check subject skip list
    %-----------------------------------------------------------------%
    if sum(strcmpi(skipsub,subnam{s})) ~= 0
        fprintf('Subject is on list of subjects to skip.  Skipping...\n')
        reason = 'Subject is on skip list.';
        skipped_log = [skipped_log;{cbusub,reason}];
        continue;
    end;
    
    % Find .mat files
    fileList = dir('S*localizer*.txt*.mat');     %%%%%%%%%%%%%% string to locate MSS output file
    if isempty(fileList)
        fprintf('No files found for %s.  Skipping...\n', cbusub)
        reason = 'No files found.';
        skipped_log = [skipped_log;{cbusub,reason}];
        continue
    end
    
    % Start file loop
    for m = 1:length(fileList)
        load(fileList(m).name);
        fileDate = fileList(m).name(end-20:end-4);

        % Check if current file has complete data
        incomplete=0;       
        for x=1:length(run_info.onsets)
            if run_info.onsets(x)==0
                incomplete=1;
                continue
            end
        end
        
        % Skip if current file is not complete
        if incomplete==1
            fprintf('Incomplete onsets found for %s.  Skipping...\n', cbusub)
            reason = sprintf('Incomplete onsets  for file %d/%d.', m, length(fileList));
            skipped_log = [skipped_log;{cbusub,reason}];
            continue
        end
        
        % Appropriate trial numbers                     %%%%%%%%%%%%%% trial numbers of each condition
        if ~isempty(strfind(fileList(m).name,'localizer1'))
                instTrialNums = [2];
                neutralTrialNums = [3 87];
                angryTrialNums =[24 108];
                fearTrialNums = [45 129];
                happyTrialNums = [66 150]; 
                localizer='localizer1';
        elseif ~isempty(strfind(fileList(m).name,'localizer2'))
                instTrialNums = [2];
                neutralTrialNums = [45 108];
                angryTrialNums =[66 129];
                fearTrialNums = [24 150];
                happyTrialNums = [3 87];
                localizer='localizer2';
        else
            fprintf('No files found for %s.  Skipping...\n', cbusub)
            reason = 'No  files found';
            skipped_log = [skipped_log;{cbusub,reason}];
            continue
        end
        
        
        % Create structs for SPM 1st level models (names, onsetes, durations)
        
        
        %%% Names
        
        names={'inst','neutral','angry','fear','happy'};        %%%%%%%%%%%%%% condition names
        %        1        2        3      4       5

        
        %%% Onsets
        onsets={};
        onsets{1}=run_info.onsets(instTrialNums);
        onsets{2}=run_info.onsets(neutralTrialNums);
        onsets{3}=run_info.onsets(angryTrialNums);
        onsets{4}=run_info.onsets(fearTrialNums);
        onsets{5}=run_info.onsets(happyTrialNums);

        
        %%% Durations
        durations={};
        
        % Each substantive block is 20 trials (2 blocks per trial type).
        % instruction duration
        durations{1}=(run_info.onsets(3)-run_info.onsets(2)); 
        
        % neutral duration
        durations{2}=[(run_info.onsets(neutralTrialNums(1)+20)-run_info.onsets(neutralTrialNums(1))) (run_info.onsets(neutralTrialNums(2)+20)-run_info.onsets(neutralTrialNums(2)))];
        
        % angry duration
        durations{3}=[(run_info.onsets(angryTrialNums(1)+20)-run_info.onsets(angryTrialNums(1))) (run_info.onsets(angryTrialNums(2)+20)-run_info.onsets(angryTrialNums(2)))];
        
        % fear duration
        durations{4}=[(run_info.onsets(fearTrialNums(1)+20)-run_info.onsets(fearTrialNums(1))) (run_info.onsets(fearTrialNums(2)+20)-run_info.onsets(fearTrialNums(2)))];
        
        % happy duration
        durations{5}=[(run_info.onsets(happyTrialNums(1)+20)-run_info.onsets(happyTrialNums(1))) (run_info.onsets(happyTrialNums(2)+20)-run_info.onsets(happyTrialNums(2)))];

        
        % Save .mat file of names, onsets, durations for each complete data file in extracted folder within sub's behav directory.
        savefile = ['extracted/',cbusub,'_Amygdala_',localizer,'_',fileDate,'_SPM.mat'];  %%%%%%%%%%%% edit task name
        save(savefile, 'names','onsets','durations');
        clear names onsets durations key_presses run_info
    end
end    
skipped_log