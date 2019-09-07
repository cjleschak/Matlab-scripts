%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  HAE --- Acquisition  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Should be run after preprocessing/DARTEL, and after badscan.m.


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

runIDs={'BOLD_Acquisition'};

taskID='HAE/Acquisition';
analysisID='Acquisition_gsrSubs_minConds_6sTrials_last3_1sShock';

funcName='sw';              % first char(s) in functional image names (for matching)
funcFormat=2;               % format of functional images: 1=3D img/hdr, 2=4D nii

% Spreadsheet containing design matrix & conrasts
contrastFile=fullfile(studyDIR,'_automation','Level1',taskID,'AcquiContrasts_minConds.xlsx');

% Other inputs
acTAG=1;                    % autocorrelation correction: 0=no, 1=yes
rpTAG=2;                    % motion regressors: 0=no, 1=yes, 2=badscan
hpf=128;                    % high-pass filter (in secs), 128 = spm default%
brainmask=fullfile(spm('Dir'),'toolbox','FieldMap','brainmask.nii'); % directory where the brainmask is

% Optional inputs - leave empty if you don't want to use them
load(fullfile(studyDIR,'_automation/Level1/HAE/whichSubs_acquiGSR.mat'));
skipsub = {};               % skip these subjects, use the full subject's folder name in single quotes (i.e. {'w001' 'w002'})
subnam = {};                % do only these subjects, format is the same as skipsup (i.e. {'w001' 'w002'})

% Shock specific inputs
early=1:7;   % which trials to model as early?
late=8:10;   % which trials to use for late trials?
shockTime=1; % Length of time (s) to model as shock?


%=====================================================================%
% Make directory in which to save job files
%---------------------------------------------------------------------%
jwd=fullfile(studyDIR,'_batches','Level1',taskID,analysisID);
eval(sprintf('mkdir %s;',jwd));


%Check specified trials.
trials=[early late];
if ~all(trials==[1:10])
    fprintf('\nERROR--Early & late trials must refer to 10 trials.\n')
    return
end
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
%g=1;
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
    cd behav
    if exist('HAE','dir')
        cd HAE
    else
        skipped_log = [skipped_log;cbusub];
        continue;
    end

    %Check for, and load, .mat data file.
    dirCheck=dir('*Acquisition*.mat');
    if isempty(dirCheck)
        continue;
    end
    load([cbusub, '_Acquisition_onsets.mat'])
    
    %Comment out if including those with messed up timing.
%     allDurs=[SSplusDurs SSminusDurs SplusDurs SminusDurs];
%     if any(allDurs~=6)
%         skipped_log=[skipped_log;cbusub];
%         fprintf('\nSome durations not 6 seconds for %s. Skipping...\n',cbusub)
%         continue
%     end
    
    %Determine shock onsets/durs
    shockOns=[SSplusOnsets+SSplusDurs-shockTime SplusOnsets+SplusDurs-shockTime];
    noShockOns=[SSminusOnsets+SSminusDurs-shockTime SminusOnsets+SminusDurs-shockTime];

    ShockDurs=shockTime;
    
   
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
            
        % Import file with contrast info
        contrastData = importdata(contrastFile);
        contrastData2=readtable(contrastFile);
        conNum=length(contrastData.textdata)-1;
        contrastNames = contrastData.textdata(2:conNum+1,1);

        condNames = contrastData2.Properties.VariableNames(3:end);
        
         
        % 1. SS plus early trials
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SSplusOnsets(early);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SSplusDurs(early)-shockTime;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;

        % 2. SS plus late trials
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SSplusOnsets(late);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SSplusDurs(late)-shockTime;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;
        
        % 3. SS minus early trials
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SSminusOnsets(early);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SSminusDurs(early)-shockTime;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;

        % 4. SS minus late trials
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SSminusOnsets(late);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SSminusDurs(late)-shockTime;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;
        
        % 5. S plus early trials
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SplusOnsets(early);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SplusDurs(early)-shockTime; %Use only first 5.5 seconds for consistentcy since modeling extra timing errors into junk.
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;

        % 6. S plus late trials
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SplusOnsets(late);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SplusDurs(late)-shockTime;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;
        
        % 7. S minus early trials
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SminusOnsets(early);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SminusDurs(early)-shockTime;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;

        % 8. S minus late trials
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = SminusOnsets(late);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = SminusDurs(late)-shockTime;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;
        
        % 9. SHOCK
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = shockOns;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = ShockDurs;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;

        % 10. No SHOCK
        condNum=condNum+1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = char(condNames(condNum));
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = noShockOns;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = ShockDurs;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;

%          % 5. junk
%          tempjunkOnsets=[SSonsets+5.5 Sonsets+5.5];
%          tempjunkDurs=[SSdurs-5.5 Sdurs-5.5];
%          index=tempjunkDurs~=0;
%          g=1;
%          
%          for ii=1:length(index)
%              if index(ii)==1
%                 junkOnsets(g)=tempjunkOnsets(ii);
%                 junkDurs(g)=tempjunkDurs(ii);
%              end
%          end
%         
%          if exist('junkOnsets','var') && exist('junkDurs','var')
%             condNum=condNum+1;
%             matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).name = 'junk';
%             matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).onset = junkOnsets;
%             matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).duration = junkDurs;
%             matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).tmod = 0;
%             matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).pmod = struct('name', {}, 'param', {}, 'poly', {});
%             matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(condNum).orth = 1;
%          else
%              fprint('\nNote that %s did not have junk onsets.\n',cbusub')
%          end
%         
%          clear junkOnsets junkDurs tempjunkOnsets tempjunkDurs
%          
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

    % Contrast Loop
    matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(cwd,'SPM.mat');
     
    for conIter = 1:conNum
        conName = char(contrastNames(conIter));
        fullContrast = contrastData.data(conIter,3:2+condNum);
        
        matlabbatch{3}.spm.stats.con.consess{conIter}.tcon.name = conName;
        matlabbatch{3}.spm.stats.con.consess{conIter}.tcon.weights = fullContrast;
        matlabbatch{3}.spm.stats.con.consess{conIter}.tcon.sessrep = 'none';  % replicate contrasts over mutliple sessions?
      
        clear conName fullContrast
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