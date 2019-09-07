function Motion1_copyFolders_noDicom(studyID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After subject folders are pulled over from the dicom using setup_subject
% to /u/project/sanscn/data/GIV_PRO/Raw_Data, you can run this script in
% order to copy over any folders which are not already copied into
% /u/project/sanscn/data/GIV_PRO/New_Motion. (It checks so duplicates are
% not copied.) 
% It copies over the raw and notes subdirectories of each folder, creates
% an analysis and behav directory for each subject folder.
% It does NOT copy over the dicom folder/files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Which folders to copy from Raw_Data to New_Motion? GIV or PRO?

%studyID='PRO';

%%%%%%%%%%%%%%%%%%
% End user input %
%%%%%%%%%%%%%%%%%%

%Copy subject folders from...
base_dir='/u/project/sanscn/data/GIV_PRO/Raw_Data'; cd(base_dir);

%Determine subject directories in base_dir, and new path to copy to.
d_subs_orig=dir(['/u/project/sanscn/data/GIV_PRO/Raw_Data',filesep,'*',studyID,'*']);
d_subs_new=fullfile('/u/project/sanscn/data/GIV_PRO/New_Motion',studyID);

%Check for already copied directories, add to skip list if already there.
doFolds=[];
skipFolds={};
for checkFold=1:length(d_subs_orig)
    if ~exist(fullfile(d_subs_new,d_subs_orig(checkFold).name),'dir')
        doFolds=[doFolds checkFold];
        fprintf('%s added to copy list.\n',d_subs_orig(checkFold).name)
    else
        fprintf('%s directory already exists in %s. Skipping...\n',d_subs_orig(checkFold).name,d_subs_new)
        skipFolds=[skipFolds d_subs_orig(checkFold).name];
    end
end


%Copy raw/notes folders of subject directories, and create analysis
%directory. Ignore dicom directory.
copied_list={};
for curFold=doFolds
     mkdir(fullfile(d_subs_new,d_subs_orig(curFold).name))
     fprintf('Copying %s''s raw folder...\n',d_subs_orig(curFold).name)
     copyfile(fullfile(d_subs_orig(curFold).folder,d_subs_orig(curFold).name,'raw'),fullfile(d_subs_new,d_subs_orig(curFold).name,'raw'));
     fprintf('\nDone. Creating other directories...\n')
     copyfile(fullfile(d_subs_orig(curFold).folder,d_subs_orig(curFold).name,'notes'),fullfile(d_subs_new,d_subs_orig(curFold).name,'notes'))
     mkdir(fullfile(d_subs_new,d_subs_orig(curFold).name,'analysis'))
     mkdir(fullfile(d_subs_new,d_subs_orig(curFold).name,'behav'))
     fprintf('%s copy complete!\n',d_subs_orig(curFold).name)
     copied_list=[copied_list;d_subs_orig(curFold).name];
end

%Print summary of what happened.
fprintf('\nThe following subject directories (minus dicom) were copied!\n')
fprintf(1, '%s\n', copied_list{:})
fprintf('\nThe following subject directories were NOT copied (already existed)!\n')
fprintf(1, '%s\n', skipFolds{:})
cd /u/project/sanscn/data/GIV_PRO/New_Motion/_automation
end
