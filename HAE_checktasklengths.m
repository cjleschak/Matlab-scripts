%%%%%% Determine WHICH SUBS have preprocessed HAE data %%%%%%
studyDIR='/u/project/sanscn/data/SAFETY'; cd(studyDIR);
outputDIR=fullfile(studyDIR,'_automation/Level1/HAE');

subID='SAS';
runIDs={'BOLD_Habituation', 'BOLD_Acquisition','BOLD_Extinction'};

% Optional inputs - leave empty if you don't want to use them
skipsub = {};               % skip these subjects, use the full subject's folder name in single quotes (i.e. {'w001' 'w002'})
subnam = {};                % do only these subjects, format is the same as skipsup (i.e. {'w001' 'w002'})



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


%=====================================================================%
% Start Subject Loop
%---------------------------------------------------------------------%
fprintf('\n\n-----HABITUATION------')
for s = dosubs

    cbusub = sprintf('%s',subnam{s});
    
    % Check subject skip list
    %-----------------------------------------------------------------%
    if sum(strcmpi(skipsub,subnam{s})) ~= 0
       % fprintf('Subject is on list of subjects to skip.  Skipping...\n')
        skipped_log = [skipped_log;cbusub];
        continue
    end
    
   % fprintf('\nStarting subject %d of %d, %s\n',s,length(dosubs),cbusub)
    swd = fullfile(studyDIR,cbusub);  % subject's working directory
    %fprintf('Subject directory is %s\n',swd)
    cd(swd)
    cd behav
    if any(exist('HAE'))
        cd HAE
    else
        continue;
    end

    base_dir = pwd;
    
    dirCheck=dir('*Habituation*.mat');
    if isempty(dirCheck)
        continue;
    end
    load([cbusub, '_Habituation_onsets.mat'])
    
    Onsets=[SSonsets Sonsets];
    Durs=[SSdurs Sdurs];
    sortedOnsets=sort(Onsets);
    lastOnset=sortedOnsets(end);
    index=find(Onsets==lastOnset);
    lastDur=Durs(index);
    taskLength=lastOnset+lastDur+14;
    fprintf('\n%s Habituation Task Length = %.2f s\n',cbusub,taskLength)
    %checkOns=sort(Onsets);checkOns(1:4)
    
    clear Onsets Durs sortedOnsets lastOnset index lastDur taskLength
    
end
fprintf('\n\n-----ACQUISITION------')
for s = dosubs

    cbusub = sprintf('%s',subnam{s});
    
    % Check subject skip list
    %-----------------------------------------------------------------%
    if sum(strcmpi(skipsub,subnam{s})) ~= 0
       % fprintf('Subject is on list of subjects to skip.  Skipping...\n')
        skipped_log = [skipped_log;cbusub];
        continue
    end
    
   % fprintf('\nStarting subject %d of %d, %s\n',s,length(dosubs),cbusub)
    swd = fullfile(studyDIR,cbusub);  % subject's working directory
    %fprintf('Subject directory is %s\n',swd)
    cd(swd)
    cd behav
    if any(exist('HAE'))
        cd HAE
    else
        continue;
    end

    base_dir = pwd;
    
        dirCheck=dir('*Acquisition*.mat');
    if isempty(dirCheck)
        continue;
    end
    load([cbusub, '_Acquisition_onsets.mat'])
    
    Onsets=[SSonsets Sonsets];
    Durs=[SSdurs Sdurs];
    sortedOnsets=sort(Onsets);
    lastOnset=sortedOnsets(end);
    index=find(Onsets==lastOnset);
    lastDur=Durs(index);
    taskLength=lastOnset+lastDur+14;
    fprintf('\n%s Acquisition Task Length = %.2f s\n',cbusub,taskLength)
    %checkOns=sort(Onsets);checkOns(1:4)
    
    clear Onsets Durs sortedOnsets lastOnset index lastDur taskLength
    
end
fprintf('\n\n-----EXTINCTION------')
for s = dosubs

    cbusub = sprintf('%s',subnam{s});
    
    % Check subject skip list
    %-----------------------------------------------------------------%
    if sum(strcmpi(skipsub,subnam{s})) ~= 0
       % fprintf('Subject is on list of subjects to skip.  Skipping...\n')
        skipped_log = [skipped_log;cbusub];
        continue
    end
    
   % fprintf('\nStarting subject %d of %d, %s\n',s,length(dosubs),cbusub)
    swd = fullfile(studyDIR,cbusub);  % subject's working directory
    %fprintf('Subject directory is %s\n',swd)
    cd(swd)
    cd behav
    if any(exist('HAE'))
        cd HAE
    else
        continue;
    end

    base_dir = pwd;
    
        dirCheck=dir('*Extinction*.mat');
    if isempty(dirCheck)
        continue;
    end
    load([cbusub, '_Extinction_onsets.mat'])
    
    Onsets=[SSonsets Sonsets];
    Durs=[SSdurs Sdurs];
    sortedOnsets=sort(Onsets);
    lastOnset=sortedOnsets(end);
    index=find(Onsets==lastOnset);
    lastDur=Durs(index);
    taskLength=lastOnset+lastDur+14;
    fprintf('\n%s Extinction Task Length = %.2f s\n',cbusub,taskLength)
    %checkOns=sort(Onsets);checkOns(1:4)
    
    clear Onsets Durs sortedOnsets lastOnset index lastDur taskLength
    
end