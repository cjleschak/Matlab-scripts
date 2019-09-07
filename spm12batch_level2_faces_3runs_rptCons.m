%=====================================================================%
% SPM12BATCH_LEVEL2  Build SPM12 Level 2 Jobs 
%
%   This script loops through subject folders and builds level 2
%   jobs using the SPM12 matlabbatch structure. User options are 
%   defined below in the section entitled "User Defined". 
%
%   July 26, 2011 -- Updated distributed to SLIC by Bob Spunt
%
%=====================================================================%
clear all; home;
spm('defaults','fmri');   % initiatizes SPM defaults for fMRI modality
spm_jobman('initcfg');    % initializes job configurations

% User Defined 
%---------------------------------------------------------------------%
executeTAG=1;   % execute job(s) immediately? 0=no, 1=yes

studyDIR = '/u/project/sanscn/data/SAFETY/'; cd(studyDIR);       % study directory
subID = 'SAS*';                                                % first digits/chars of subject directories, with wildcard

taskID = 'Faces';
analysis_dir = 'Faces_hpf128ar1_anyRuns_rptCons_07182019';                                

groupstats_folder='_groupstats';                               % name of your groupstats folder
batch_dir = [studyDIR '_batches' filesep 'Level2' filesep];

explicit_mask = '/u/project/sanscn/data/SAFETY/_automation/ROI/FieldMap_greyMatter_20pct.nii';

% Optional inputs - leave empty if you don't want to use them
skipsub = {'SAS_358','SAS_416','SAS_433'}';               % skip these subjects, use the full subject's folder name in single quotes (i.e. {'w001' 'w002'})
subnam = {};                % do only these subjects, format is the same as skipsup (i.e. {'w001' 'w002'})


%---------------------------------------------------------------------%

% Make directories (if needed)
%---------------------------------------------------------------------%
jwd = [batch_dir taskID filesep analysis_dir];
gwd = [studyDIR groupstats_folder filesep taskID filesep analysis_dir '_2or3Runs'];
eval(sprintf('mkdir %s;',jwd));
eval(sprintf('mkdir %s;',gwd));

% Find subject directories
%---------------------------------------------------------------------%
d=dir(subID);
for i=1:length(d)
    subnam{i}=d(i).name;
    fprintf('Adding %s\n',subnam{i})
end

subnum = length(subnam);
dosubs = 1:subnum;

conNum = 0;
conCount = 1;
contrastInfo.name = {};
contrastInfo.file = {};
conNameList = {};

% Load SPM.mat for subjects
%---------------------------------------------------------------------%
for s = dosubs
    cbusub = sprintf('%s',subnam{s});
    fwd = sprintf('%s%s%sanalysis%s%s%s%s',studyDIR,cbusub,filesep,filesep,taskID,filesep,analysis_dir);
    if sum(ismember(skipsub,cbusub))>0
        fprintf('Subject %s is in exclusion list.  Skipping...\n',subnam{s})
        continue;
    end
    
    if ~exist(fwd,'dir')
        fprintf('Analysis directory does not exist for subject %s.  Skipping...\n',subnam{s})
        skipsub = [skipsub; cbusub];
        continue;
    end
    
    cd(fwd)
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
        imageFilename = sprintf('con_%04d.nii',i);
        fullImageFilename = strcat(fwd,filesep,imageFilename);
        contrastInfo.name{conCount,1} = SPM.xCon(i).name;
        contrastInfo.file{conCount,1} = fullImageFilename;
        conCount = conCount +1;
    end
    
    % check to see if current subject has most contrasts (to be used 
    % as comprehensive list of all contrasts in analysis later)
    if length(SPM.xCon) > conNum
        clear conNameList
        conNum = length(SPM.xCon);
        nameCount = 1;
        for j = 1:length(SPM.xCon)
            conNameList{nameCount} = SPM.xCon(j).name;
            nameCount = nameCount +1;
        end
    end
end

% find and remove the Omnibus test from being run at the second level
findOmni = regexp(conNameList,'Omnibus');
lineCount = 1;
omniIdx = 0;
for i = 1:length(findOmni)
    if isempty(findOmni{i})
        omniIdx(lineCount,:) = 0;
    else
        omniIdx(lineCount,:) = 1;
    end
    lineCount = lineCount +1;
end
conNameList(logical(omniIdx)) = [];
conNum = length(conNameList);

% Populate list of con image files
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTRAST LOOP                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for c = 1:conNum
    current_conname=conNameList{c};
    current_conname_nosuffix=regexprep(conNameList{c},' - All Sessions','');
    current_dwd2=[gwd filesep current_conname_nosuffix '_' analysis_dir];
    current_dwd = current_dwd2(current_dwd2 ~= ' ');
    eval(sprintf('mkdir %s;',current_dwd));
    
    % now populate list of contrast path/filenames for current contrast
    current_conidx = find(ismember(contrastInfo.name, current_conname)==1);
    conimg_filenames = contrastInfo.file(current_conidx);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE MATLABBATCH STRUCTURE                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    matlabbatch{1}.spm.stats.factorial_design.dir{1} = current_dwd;
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = strcat(conimg_filenames,',1');
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = strcat(explicit_mask,',1');
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(current_dwd,'SPM.mat');
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(current_dwd,'SPM.mat');
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Positive';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Negative';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;

    % Save job
    %-----------------------------------------------------------------%
    a=datestr(clock,31);   % returns date string of the form 'YYYY-MM-DD HH:MM:SS' e.g., 2006-12-27 15:03:37
    time_stamp = [a(6:7) a(9:10) a(3:4) '_' a(12:13) a(15:16)];   % timestamp is a function name, hence the _ in time_stamp
    filename=sprintf('%s/%s_%s_%s',jwd,current_conname_nosuffix,analysis_dir,time_stamp); 
    save(filename,'matlabbatch');   % save jobs variable
    
    % Execute job (if applicable)
    %-----------------------------------------------------------------%
    if executeTAG==1
        spm_jobman('run',matlabbatch);
    end
    
    clear matlabbatch conimg_filenames
    
end
skipsub