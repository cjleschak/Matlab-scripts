function Motion3_wrapper_reslice(studyID)
%% Parameters

% Folder/directory information
owd=fullfile('/u/project/sanscn/data/GIV_PRO/New_Motion',studyID);  % study directory
codeDir = '/u/project/sanscn/data/GIV_PRO/New_Motion/_automation'; % where code lives
output='/u/project/sanscn/data/GIV_PRO/New_Motion/_batches'; % dir in which to save batches
subID=[studyID '_'];      % pattern for finding subject folders (wildcard is used below, so put whatever you want to appear before the * - e.g. GIV_1 for GIV_1*)
load([codeDir,filesep,[studyID,'subjects2realign.mat']])%{'GIV_038A','GIV_038B','GIV_039A','GIV_039B','GIV_144A','GIV_146A','GIV_146B'};
doSubs=needsRealign;
%subNam = {};     % do which subjects? ('all' to do all, position vector, e.g. 1:4, to do a subset)
skipSub = {};
runID='BOLD_*';     % pattern for finding functional run folders (use wildcards)
% mbwdirID='Matched_Bandwidth_HiRes*';    % pattern for finding matched-bandwidth folder (use wildcards)
mpragedirID='MPRAGE_*';  % pattern for finding mprage folder (use wildcards)
addpath(output);
%addpath('/u/project/sanscn/data/MINERVA2/nondartel/preproc/170227');
%% Setup subjects

%CJL changed to take actual sub names as do subs, rather than position
%vector

%if isempty(subNam)
d = dir([owd '/' subID '*']);
addpath(d.folder);
for ii = 1:length(d)
    subNam{ii} = d(ii).name;
    %fprintf('Adding %s\n', subNam{ii})
end
subNam=intersect(subNam,doSubs);
%end
numSubs = length(subNam);
cd(codeDir);

% Prepare status struct
runStatus = struct('subNam',[],'status',[],'error',[]);
runStatus(numSubs).subNam = [];
runStatus(numSubs).status = [];
runStatus(numSubs).error = [];

%% Run subjects

%if you want to run just one subject, comment out pool & parfor, leave for
%uncommented -->opposite for multiple

nWorkers = maxNumCompThreads; % specify integer if desired
nWorkers = min(numSubs, nWorkers);
parpool('local', nWorkers);

%for i = 1:numSubs
parfor i = 1:numSubs
    
    runStatus(i).subNam = subNam{i};
    
%     if ismember(subNam{i}, skipSub) % Skip subjects or not
%         runStatus(i).status = 0;
%         runStatus(i).error = 'Subject in exclusion list';
%         disp(['Skipping subject ' subNam{i} ', is in exclusion list']);
%         continue
%     else % Run subject
        fprintf('\nRunning subject %s\n', subNam{i});
        try % Try running subfunction
            [runStatus(i).status, runStatus(i).error]...
                = run_motion_preprocessing_reslice(subNam{i}, owd, codeDir, output, runID, mpragedirID);
            if runStatus(i).status == 1
                fprintf('\nsubject %s successful\n', subNam{i});
            elseif runStatus(i).status == 0
                disp([runStatus(i).error ' for ' subNam{i}]);
            end
        catch % Log if fail
            runStatus(i).status = 0;
            runStatus(i).error = 'Unexpected error in run function';
            disp(['Unexpected ERROR on subject ' subNam{i}]);
        end
end
delete(gcp('nocreate'));

% Save stuff
date = datestr(now,'yyyymmdd_HHMM');
filename = [output '/runStatus_' date '.mat'];
save(filename,'runStatus');
filename = [output '/workspace_' date '.mat']; % Use this to re-do "run dartel" if it fails
save(filename);

end