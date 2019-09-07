%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Faces Task - Runs 1, 2, and 3 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Should be run after preprocessing/DARTEL, and after badscan.m.
%Behavioral output .csv files which have timing info should be 
%located in each subject's 'behav' folder.

%=====================================================================%
clear all; home;
spm('defaults','fmri');     % initializes SPM defaults for fMRI modality
spm_jobman('initcfg');      % initializes job configurations


%=====================================================================%
% USER DEFINED
%---------------------------------------------------------------------%
executeTAG=1;           % execute job(s) immediately? 0=no, 1=yes

studyDIR='/u/project/sanscn/data/SAFETY'; cd(studyDIR);
subID='SAS';

runIDs={'BOLD_Faces1','BOLD_Faces2','BOLD_Faces3'};

taskID='Faces';
analysisID='Faces_hpf90ar1_3runs';

funcName='sw';              % first char(s) in functional image names (for matching)
funcFormat=2;               % format of functional images: 1=3D img/hdr, 2=4D nii

% Spreadsheet containing design matrix & conrasts
%contrastFile=[studyDIR '/_automated/Level1/Faces123Contrasts.xlsx'];
contrastFile=fullfile(studyDIR,'_automation','Level1',taskID,'all3runs','facesContrasts.xlsx');

% Other inputs
acTAG=1;                    % autocorrelation correction: 0=no, 1=yes
rpTAG=2;                    % motion regressors: 0=no, 1=yes, 2=badscan
hpf=90;                    % high-pass filter (in secs), 128 = spm default%
brainmask='/u/project/data/SAFETY/_automation/ROI/FieldMap_greyMatter_20pct.nii'; % directory where the brainmask is

% Optional inputs - leave empty if you don't want to use them
skipsub = {};               % skip these subjects, use the full subject's folder name in single quotes (i.e. {'w001' 'w002'})
subnam = {};                % do only these subjects, format is the same as skipsup (i.e. {'w001' 'w002'})



%=====================================================================%
% Make directory in which to save job files
%---------------------------------------------------------------------%
jwd=fullfile(studyDIR,'_batches','Level1',taskID,analysisID);
eval(sprintf('mkdir %s;',jwd));


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
skipped_log = {};


%=====================================================================%
% Start Subject Loop
%---------------------------------------------------------------------%
g=1;
for s = dosubs
    
    cbusub = sprintf('%s',subnam{s});
    fprintf('\nStarting subject %d of %d, %s\n',s,length(dosubs),cbusub)
    swd = fullfile(studyDIR,cbusub);  % subject's working directory
    fprintf('Subject directory is %s\n',swd)
    cd(swd)
    cd raw
    base_dir = pwd;     % subject's base data directory (raw folder)
    
    % Check subject skip list
    %-----------------------------------------------------------------%
    if sum(strcmpi(skipsub,subnam{s})) ~= 0
        fprintf('Subject is on list of subjects to skip.  Skipping...\n')
        skipped_log = [skipped_log;cbusub];
        continue;
    end;
    
    % Find run directories
    %-----------------------------------------------------------------%

    d={};
    for b=1:length(runIDs)
        d=[d dir([runIDs{b} '*'])];
    end
    
    clear run_names

    if isempty(d)
        numruns=0;
    else
        for i=1:length(d)
            run_names{i}=d{i}.name;
            numruns=length(d);
        end;
    end
    
    %If subject does not have 3 run directories, skip.
    if numruns~=3
        fprintf('%d runs found for this subject.  Skipping...\n', numruns)
        skipped_log = [skipped_log;cbusub];
        continue
    end
    
    
    % Find functional images for run(s)
    %-----------------------------------------------------------------%
    cd(swd)
    cd raw
    for i = 1:numruns
        load_dir = fullfile(base_dir,run_names{i});
        if funcFormat==1
            tmpString=sprintf('^%s.*\\.img$',funcName);
            [raw_func_filenames{i},dirs] = spm_select('List',load_dir,tmpString);
        else
            tmpString=sprintf('^%s.*\\.nii$',funcName);
            [raw_func_filenames{i},dirs] = spm_select('ExtList',load_dir,tmpString,1:10000);
        end
        eval(sprintf('filenames%d=cellstr(strcat(load_dir,filesep,raw_func_filenames{i},'',1''));',i));
        if rpTAG==2
            rp_file{i} = spm_select('List',load_dir,'^badscan.*\.txt$');
        elseif rpTAG==1
            rp_file{i} = spm_select('List',load_dir,'^rp.*\.txt$');
        end
    end
    
    % Find behav data output files for each run
    %-----------------------------------------------------------------%
    cd(swd);
    cd behav    %cd into sub's behav folder where data output files are
    cd Faces

    %Find .txt data files for this subject
    dataFilesCSV=[];
    dataFiles=[];
    for k=1:length(runIDs)
        dataFilesCSV=[dataFilesCSV dir([cbusub '_Faces*Run' num2str(k) '*.csv'])];
        copyfile(dataFilesCSV(k).name,[dataFilesCSV(k).name(1:end-3) 'txt'])
        dataFiles=[dataFiles dir([cbusub '_Faces*Run' num2str(k) '*.txt'])];
    end

    %How many data files were found for this subject?
    if isempty(dataFiles)
        numDataFiles=0;
    else
        for i=1:length(dataFiles)
            file_names{i}=dataFiles(i).name;
            numDataFiles=length(dataFiles);
        end;
    end
    
    %If subject does not have 3 .txt data files, skip.
    if numDataFiles~=3
        fprintf('%d .txt data files found for this subject.  Skipping...\n', numDataFiles)
        skipped_log = [skipped_log;cbusub];
        continue
    end
    
    
    % Make a folder for the analysis in the current subject's folder
    %-----------------------------------------------------------------%
    cwd = fullfile(swd,'analysis',taskID,analysisID);
    eval(sprintf('mkdir %s;',cwd));

      
    matlabbatch{1}.spm.stats.fmri_spec.dir{1} = cwd;            % where will SPM.mat file be written?
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';   % specified in either scans or seconds
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;           % TR, specified in seconds
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;      % default is 16, doesn't need to be changed unless slice-timing was used
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;      % default in SPM 12 is now 8 (SPM 8 default = 1)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START RUN LOOP                                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %EXTRACT DATA FROM OUTPUT FILE
    
    %Preallocate structs & clear onsets etc. from prev. participant.
    HEALTHY=repmat(struct('run',[],'trials',[],'ratings',[],'faceOnset',[], ...
        'faceOffset',[],'faceDuration',[],'rateOnset',[],'rateOffset',[], ...
        'rateDuration',[]),3,1); 
    HEALTHY(1).run=1; HEALTHY(2).run=2; HEALTHY(3).run=3;
    SICK=repmat(struct('run',[],'trials',[],'ratings',[],'faceOnset',[], ...
        'faceOffset',[],'faceDuration',[],'rateOnset',[],'rateOffset',[], ...
        'rateDuration',[]),3,1); 
    SICK(1).run=1; SICK(2).run=2; SICK(3).run=3;
    
    
    for u=1:numDataFiles                            %for each run
        curRunData=readtable(dataFiles(u).name,'Delimiter',',');    %load appropriate data file for current run
        t1=1; t2=1;                                            
        for n=1:height(curRunData)                  %for each row of the file (e.g. each trial
            if curRunData.Face_condition(n)==0      %if healthy trial
                HEALTHY(u).trials=[HEALTHY(u).trials curRunData.Trial(n)];
                HEALTHY(u).ratings=[HEALTHY(u).ratings curRunData.Rating_final(n)];
                HEALTHY(u).faceOnset=[HEALTHY(u).faceOnset curRunData.Face_onset(n)];
                HEALTHY(u).faceOffset=[HEALTHY(u).faceOffset curRunData.Face_Offset(n)];
                HEALTHY(u).faceDuration=[HEALTHY(u).faceDuration curRunData.Face_Duration(n)];
                HEALTHY(u).rateOnset=[HEALTHY(u).rateOnset curRunData.Rating_Onset(n)];
                HEALTHY(u).rateOffset=[HEALTHY(u).rateOffset curRunData.Rating_Offset(n)];
                HEALTHY(u).rateDuration=[HEALTHY(u).rateDuration curRunData.Rating_Duration(n)];
                t1=t1+1;
            elseif curRunData.Face_condition(n)==1   %if sick trial
                SICK(u).trials=[SICK(u).trials curRunData.Trial(n)];
                SICK(u).ratings=[SICK(u).ratings curRunData.Rating_final(n)];
                SICK(u).faceOnset=[SICK(u).faceOnset curRunData.Face_onset(n)];
                SICK(u).faceOffset=[SICK(u).faceOffset curRunData.Face_Offset(n)];
                SICK(u).faceDuration=[SICK(u).faceDuration curRunData.Face_Duration(n)];
                SICK(u).rateOnset=[SICK(u).rateOnset curRunData.Rating_Onset(n)];
                SICK(u).rateOffset=[SICK(u).rateOffset curRunData.Rating_Offset(n)];
                SICK(u).rateDuration=[SICK(u).rateDuration curRunData.Rating_Duration(n)];
                t2=t2+1;
            end
        end
    end

    %Compute and save aggregate ratings
    %preallocate if first subject
    if g==1
        RATINGS=repmat(struct('SubjectID',[],'healthyRun1',[],'sickRun1',[], ...
            'healthyRun2',[],'sickRun2',[],'healthyRun3',[],'sickRun3',[], ...
            'healthyOverall',[],'sickOverall',[]),subnum,1);
    end
    
    for kk=1:3
        if ~isfloat(HEALTHY(kk).ratings)
            HEALTHY(kk).ratings=str2double(HEALTHY(kk).ratings);
        end
        if ~isfloat(SICK(kk).ratings)
            SICK(kk).ratings=str2double(SICK(kk).ratings);
        end
    end
    
    RATINGS(g).SubjectID=cbusub;
    RATINGS(g).healthyRun1=nanmean(HEALTHY(1).ratings);
    RATINGS(g).healthyRun2=nanmean(HEALTHY(2).ratings);
    RATINGS(g).healthyRun3=nanmean(HEALTHY(3).ratings);
    RATINGS(g).sickRun1=nanmean(SICK(1).ratings);
    RATINGS(g).sickRun2=nanmean(SICK(2).ratings);
    RATINGS(g).sickRun3=nanmean(SICK(3).ratings);
    RATINGS(g).healthyOverall=nanmean([HEALTHY.ratings]);
    RATINGS(g).sickOverall=nanmean([SICK.ratings]);
    g=g+1;
   
    for r = 1:numruns
        
        % Specify scans and directory for current run
        %-------------------------------------------------------------%
        eval(sprintf('current_filenames=filenames%d;',r));
        current_load_dir = fullfile(base_dir,run_names{r});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DESIGN SPECIFICATION MODULE                                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = current_filenames;   % load all scans of each run of the current task for that subject

        % Predictors
        %-------------------------------------------------------------%
        
        condNum=0;
        
        % 1. healthy faces 1, 2, 3 (based on iteration of r loop)
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = ['healthyFaces' num2str(r)];
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = HEALTHY(r).faceOnset;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = HEALTHY(r).faceDuration;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;

        % 2. healthy ratings 1, 2, 3 (based on iteration of r loop)
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = ['healthyRatings' num2str(r)];
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = HEALTHY(r).rateOnset;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = HEALTHY(r).rateDuration;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;
        
        % 3. sick faces 1, 2, 3 (based on iteration of r loop)
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = ['sickFaces' num2str(r)];
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SICK(r).faceOnset;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SICK(r).faceDuration;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;

        % 4. sick ratings 1, 2, 3 (based on iteration of r loop)
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = ['sickRatings' num2str(r)];
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SICK(r).rateOnset;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SICK(r).rateDuration;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;
        
       
        
        % Additional session-specific parameters
        %-------------------------------------------------------------%
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {''};

        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        if rpTAG==1 || rpTAG==2
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg{1} = strcat(current_load_dir,filesep,rp_file{r});
        elseif rpTAG==0
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg{1} = '';
        end;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = hpf;
        
    end;

    % Final inputs to design specification (session non-specific)
    %-----------------------------------------------------------------%
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs(1) = 0; % time derivative (0=no, 1=yes)
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs(2) = 0; % dispersion derivative (0=no, 1=yes)
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask{1} = brainmask;
    if acTAG==1
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    else
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DESIGN ESTIMATION MODULE                                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(cwd,'SPM.mat');
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONTRAST MANAGER MODULE                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Figure out # of motion regressors added from badscan to go in between
    % each run's contrast weights
    for r=1:numruns
        load(fullfile(base_dir,run_names{r},rp_file{r}))
        badscanCols(r)=size(badscan,2);
    end
    rpMatrix1=zeros(1,badscanCols(1));
    rpMatrix2=zeros(1,badscanCols(2));
    rpMatrix3=zeros(1,badscanCols(3));
    
        
    % Import file with contrast info
    contrastData = importdata(contrastFile);
    conNum=length(contrastData.textdata)-1;
    contrastNames = contrastData.textdata(2:conNum+1,1);

  
    % Contrast Loop
    matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(cwd,'SPM.mat');
     
    for conIter = 1:conNum
        conName = char(contrastNames(conIter));
        run1Matrix = contrastData.data(conIter,3:6);
        run2Matrix = contrastData.data(conIter,7:10);
        run3Matrix = contrastData.data(conIter,11:14);
        fullContrast = [run1Matrix rpMatrix1 run2Matrix rpMatrix2 run3Matrix rpMatrix3];
        
        matlabbatch{3}.spm.stats.con.consess{conIter}.tcon.name = conName;
        matlabbatch{3}.spm.stats.con.consess{conIter}.tcon.weights = fullContrast;
        matlabbatch{3}.spm.stats.con.consess{conIter}.tcon.sessrep = 'none';  % replicate contrasts over mutliple sessions?
      
        clear conName run1Matrix run2Matrix run3Matrix fullContrast
    end
    
    matlabbatch{3}.spm.stats.con.delete = 0;                        % option to delete existing contrasts (0=no, 1=yes)

    % Save job
    %-----------------------------------------------------------------%
    a=datestr(clock,31);   % returns date string of the form 'YYYY-MM-DD HH:MM:SS' e.g., 2006-12-27 15:03:37
    time_stamp = [a(6:7) a(9:10) a(3:4) '_' a(12:13) a(15:16)];   % timestamp is a function name, hence the _ in time_stamp
    filename=fullfile(jwd,[cbusub '_' time_stamp]);
    save(filename,'matlabbatch');
    
    
    % Execute job (if applicable)
    %-----------------------------------------------------------------%
    if executeTAG==1
        spm_jobman('run',matlabbatch);
    end
    
    % Clear some variables before the next iteration of the loop
    %-----------------------------------------------------------------%
    clear matlabbatch *Ons rt* *Dur
    
end;
a=datestr(clock,31);   % returns date string of the form 'YYYY-MM-DD HH:MM:SS' e.g., 2006-12-27 15:03:37
time_stamp = [a(6:7) a(9:10) a(3:4) '_' a(12:13) a(15:16)];   % timestamp is a function name, hence the _ in time_stamp
filename2=fullfile(studyDIR,'_automation/Level1',taskID,['facesRatings_' time_stamp]);
    save(filename2,'RATINGS')
skipped_log