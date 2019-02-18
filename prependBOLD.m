clear all; home;


%%%%%%%%%%%%%%%%%%%%%%%%%% prependBOLD.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script looks through all of your subjects' raw folders
% and prepends them with 'BOLD_' so that functional runs can be easily 
% identified later in analyses and with other scripts (like DARTEL).
% Note that it prepends 'BOLD' to both the folder AND the .nii files within.

% Created by CJL 9/12/2017

%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input study directory (location of subject folders)
owd = '/u/project/sanscn/data/SSBCS'; cd(owd); 

% input first digits/chars of subject directories (no wildcard)
subid = 'SSPB';    

% input first digits/chars of task names (WITH wildcards) that you want to prepend (e.g. {'Amyg*' 'Math*'})
taskNames = {'Amyg*' 'Math*' 'Cyberball*' 'Shapes*' 'Resting*' 'Lottery*' 'CO*'};                            

% Optional inputs - leave empty if you don't want to use them
skipsub = {};                % skip these subjects, use the full subject's folder name in single quotes (i.e. {'w001' 'w002'})
subnam = {};                 % do only these subjects, format is the same as skipsup (i.e. {'w001' 'w002'})


%%%%%%%%%%%%%%%%%%%%%%%%% END USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('subnam','var')    % if the subnam variable exists, but is empty
    if isempty(subnam) 
        d=dir([subid '*']); % then d = a list of all the directories that start with subid 
        for i=1:length(d)   % (autopopulate with all subjects in directory matching subid prefix)
            subnam{i}=d(i).name; 
            fprintf('Adding %s\n',subnam{i})
        end
    end
else                        % if the subnam variable does not exist
    d=dir([subid '*']);     % then d = a list of all the directories that start with subid
    for i=1:length(d)       % (autopopulate with all subjects in directory matching subid prefix)
        subnam{i}=d(i).name;
        fprintf('Adding %s\n',subnam{i})
    end
end                         % if the subname variable does exist and is not empty, then subnam remains as entered

subnum = length(subnam);    % how many subjects in subject list?
dosubs = 1:subnum;          % do subs 1 through however many subjects are in the list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:subnum
    sub = dosubs(s);
    cbusub = sprintf('%s',subnam{sub});

    if sum(ismember(skipsub,cbusub))>0
        fprintf('Subject %s is in exclusion list.  Skipping...\n',subnam{s})
        continue;
    end
    
    fprintf('\n\nStarting subject %d of %d, %s\n',s,length(dosubs),cbusub)
    swd = fullfile(owd,cbusub);
        fprintf('Subject directory is %s\n',swd)
    cd([swd, '/raw'])           % cd into current subject's raw folder.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for t = 1:length(taskNames)
        cd([swd, '/raw'])
        currentTask = taskNames{t};
        runList=dir(currentTask);
        if ~isempty(dir(currentTask))
            for r = 1:length(dir(currentTask))
                movefile(runList(r).name,['BOLD_' runList(r).name]);
                cd(['BOLD_' runList(r).name]);
                movefile([runList(r).name '.nii'],['BOLD_' runList(r).name '.nii']);
                cd([swd, '/raw'])
            end
        end
    end
end