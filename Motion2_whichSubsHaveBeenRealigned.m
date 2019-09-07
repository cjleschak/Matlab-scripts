function Motion2_whichSubsHaveBeenRealigned(studyID)
% Afer new subject folders have been copied over to 
% /u/project/sanscn/data/GIV_PRO/New_Motion, use this script to identify
% subject folders which need to be realigned using nondartel scripts.
% Identifies folders which need to be realigned by checking for 1) rp.txt
% file in each BOLD folder, realigned/resliced 4D BOLD image in each BOLD 
% folder, and 3) spm.ps output file in the notes folder, for each subject.

base_dir='/u/project/sanscn/data/GIV_PRO/New_Motion';
%studyID='GIV';
subID=[studyID '_*']; %include wildcards

% If you include the Run# in the runID, it will detect instances where 
% there are more than one run (e.g. false starts) - which may cause the 
% recreateMotionPlot script to fail. e.g. {'BOLD_RL_Run1*','BOLD_GIV_Run2*'}
runIDs={'BOLD_GIV_Run1*','BOLD_GIV_Run2*','BOLD_GIV_Run3*'};%,'BOLD_RL_Run1*','BOLD_RL_Run2*','BOLD_MathDecisions*Run1*','BOLD_MathDecisions*Run2*'}; 



%%%%%%%%%%%%%%%%%%
% End user input %
%%%%%%%%%%%%%%%%%%


%Determine BOLD directories.
d_sub_dirs=dir(fullfile(base_dir,studyID,subID));
needsRealign={};
RealignDone={};
hasMultRuns={};
for eachSub=1:length(d_sub_dirs)
    realigned=1;
    only1run=1;
    noBold=0;
    curSub=d_sub_dirs(eachSub).name;
    fprintf('\nChecking %s...',curSub)
    cd(fullfile(d_sub_dirs(eachSub).folder,d_sub_dirs(eachSub).name))
    cd raw
    d_BOLD_dirs=dir('BOLD*');
    if ~isempty(d_BOLD_dirs)
        for x=1:length(runIDs)
            if length(dir(char(fullfile(d_BOLD_dirs(1).folder,char(runIDs(x))))))>1
                fprintf('Warning: Multiple BOLD runs for %s\n',char(runIDs(x)))
                only1run=0;
            end
        end
        for eachDir=1:length(d_BOLD_dirs)
            cd(fullfile(d_BOLD_dirs(eachDir).folder,d_BOLD_dirs(eachDir).name))
            if ~exist(['rp_',d_BOLD_dirs(eachDir).name,'.txt'],'file')
                fprintf('Warning: No rp*.txt file for %s\n',d_BOLD_dirs(eachDir).name)
                realigned=0;
            end
            if ~exist(['r',d_BOLD_dirs(eachDir).name,'.nii'],'file')
                fprintf('Warning: No realigned/resliced .nii file for %s\n',d_BOLD_dirs(eachDir).name)
                realigned=0;
            end
        end
        cd(fullfile(d_sub_dirs(eachSub).folder,d_sub_dirs(eachSub).name))
        cd notes
        if isempty(dir('spm*.ps'))
            fprintf('Warning: No post-script file for %s\n',curSub)
            realigned=0;
        elseif length(dir('spm*.ps'))>1
            fprintf('Warning: More than 1 post-script file for %s\n',curSub)
        end
    else
        fprintf('Warning: No BOLD .nii files for %s.\n',curSub)
        noBold=1;
    end
    
    if realigned==0
       needsRealign=[needsRealign,curSub]; 
    else
        if noBold==0
           RealignDone=[RealignDone,curSub];
        end
    end
    if only1run==0
        hasMultRuns=[hasMultRuns,curSub];
    end
end
cd /u/project/sanscn/data/GIV_PRO/New_Motion/_automation

fprintf('\n\nThe following subjects appear to have been realigned already:\n')
fprintf(1, '%s\n', RealignDone{:})

if ~isempty(needsRealign)
    fprintf('\n\nThe following subjects appear to be missing realignment files:\n')
    fprintf(1, '%s\n', needsRealign{:})
    save([studyID,'subjects2realign.mat'],'needsRealign')
else
    fprintf('\n\nAll %s subjects in /u/project/sanscn/data/GIV_PRO/New_Motion have been realigned.\n', studyID)
    delete([studyID,'subjects2realign.mat'],'needsRealign')
end

if ~isempty(hasMultRuns)
    fprintf('\n\nWARNING: The following subjects appear to have multiple runs of a single task.\nPut an underscore at the front of the folder name of the run that will not be analyzed before proceeding:\n')
    fprintf(1, '%s\n', hasMultRuns{:})
end


end