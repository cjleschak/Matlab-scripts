%% Fear Conditioning Task %%
% This script will first index through the subject's raw folders
% in the Hoffman2 study folder (studyDIR), and determine which subjects have
% data for the runIDs specified (runIDs).
clear all;

runWhichSubs=0;
runDataExtract=1;

%%%%%% Determine WHICH SUBS have preprocessed HAE data %%%%%%
studyDIR='/u/project/sanscn/data/SAFETY'; cd(studyDIR);
outputDIR=fullfile(studyDIR,'_automation/Level1/HAE');

subID='SAS';
runIDs={'BOLD_Habituation', 'BOLD_Acquisition','BOLD_Extinction'};

% Optional inputs - leave empty if you don't want to use them
skipsub = {};               % skip these subjects, use the full subject's folder name in single quotes (i.e. {'w001' 'w002'})
subnam = {};                % do only these subjects, format is the same as skipsup (i.e. {'w001' 'w002'})

%=====================================================================%
% Find subject directories
%---------------------------------------------------------------------%
if runWhichSubs==1
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
    % Make sure Sub variables are cleared.
    HabitSubs=[];
    AcquiSubs=[];
    ExtinSubs=[];
    
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
       

        %Check directories
        dHabit=dir([runIDs{1} '*']);
        dAcqui=dir([runIDs{2} '*']);
        dExtin=dir([runIDs{3} '*']);

        %If subjects have usable BOLD Data, add to list of subs for each.
        if ~isempty(dHabit)
            HabitSubs=[HabitSubs;cbusub];
        end
        if ~isempty(dAcqui)
            AcquiSubs=[AcquiSubs;cbusub];
        end
        if ~isempty(dExtin)
            ExtinSubs=[ExtinSubs;cbusub];
        end

    end
    whichSubsFile=fullfile(outputDIR,'whichSubs_HAE.mat');
    save(whichSubsFile,'HabitSubs','AcquiSubs','ExtinSubs')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Extract onsets from behav files %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if runDataExtract==1
    %trigDiffs=[];
    warning('off','MATLAB:table:ModifiedVarnames')
    orderVersionInfo=readtable(fullfile(studyDIR,'_automation/Level1/HAE/HAE order info.xlsx'));
    
    % Check for whichSubs.mat file, created by first half of script.
    if ~exist(fullfile(outputDIR,'whichSubs_HAE.mat'))
        fprintf('Cannot find whichSubs.mat file. Check that whichSubs was run, or recheck filepath.')
        return
    end
    
    % Specify which rows of output files correspond to which run.
    load(fullfile(outputDIR,'whichSubs_HAE.mat'))
    for rr=1:length(runIDs)
        if strcmp(runIDs(rr),'BOLD_Habituation')
            whichSubsVar=cellstr(HabitSubs);
            curRunRows=1:16;
        elseif strcmp(runIDs(rr),'BOLD_Acquisition')
            whichSubsVar=cellstr(AcquiSubs);
            curRunRows=17:56;
        elseif strcmp(runIDs(rr),'BOLD_Extinction')
            whichSubsVar=cellstr(ExtinSubs);
            curRunRows=58:81;
        else
            fprintf('Error - current run ID not asociated with whichSubs file.')
            return
        end
        curRunName=char(runIDs(rr)); curRunName=curRunName(6:end);
        
        % Start subject loop.
        for ss=1:length(whichSubsVar)
            cbusub=char(whichSubsVar(ss));
                
           % Check run-these-subjects-only list.
                if exist('subnam') && ~isempty(subnam)
                    if ~any(strcmpi(cbusub,subnam))
                        fprintf('\nSkipping %s. Sub not on subnam var.\n',cbusub)
                        continue
                    end
                end
                % Check skip-these-subjects list.
                if any(strcmpi(cbusub,skipsub))
                    fprintf('\nSkipping %s. On skip list.\n',cbusub)
                    continue
                end
                
             cd(fullfile(studyDIR,cbusub,'behav/HAE'))
             subsData=readtable(['SST_fMRI_merged-' cbusub(5:end) '-1.xlsx']);
             
             % Determine subject's order.
             curOrderInfoRow=find(strcmp(orderVersionInfo.SubjectID,cbusub));
             curOrder=orderVersionInfo.FullOrder(curOrderInfoRow);
             
             % Determine onset of task
             if strcmp(runIDs(rr),'BOLD_Habituation')
                taskOnset=str2num(char(subsData.B1_OnsetTime(1)));
                A1rows=[3 10 12 14];
                A2rows=[];
                A3rows=[2 7 8 16];
                B1rows=[1 4 6 13];
                B2rows=[];
                B3rows=[5 9 11 15];
                A1='Aminus';
                B1='Bminus';
                A3='Aplus';
                B3='Bplus';
             elseif strcmp(runIDs(rr),'BOLD_Acquisition')
                taskOnset=str2num(char(subsData.A1_OnsetTime(17)));
                A1rows=[17 20 26 29 32 38 39 42 47 51];
                A2rows=[18 23 24 28 33 36 43 44 46 56];
                A3rows=[];
                B1rows=[19 21 25 27 35 41 45 49 52 53];
                B2rows=[22 30 31 34 37 40 48 50 54 55];
                B3rows=[];
                A1='Aminus';
                B1='Bminus';
                A2='Aplus';
                B2='Bplus';
             elseif strcmp(runIDs(rr),'BOLD_Extinction')
                taskOnset=str2num(char(subsData.A1_OnsetTime(58)));
                A1rows=[58 65 70 75 77 79];
                A2rows=[];
                A3rows=[60 64 67 69 71 74];
                B1rows=[61 66 68 73 78 80];
                B2rows=[];
                B3rows=[59 62 63 72 76 81];    
                A1='Aminus';
                B1='Bminus';
                A3='Aplus';
                B3='Bplus';
             else
                fprintf('\nError for %s - can''t identify task onset.\n',cbusub)
                return
             end
             
        
             %A1 and B1 = not reinforced
             %A2 anad B2 = reinforced during acqui.
             %A3 and B3 = not reinforced
             
             %NRa1 = A1 = Aminus.jpg = 6s
             %NRb1 = B1 = Bminus.jpg = 6s
             %NRa2 = A3 = Aplus.jpg = 6s
             %NRb2 = B3 = Bplus.jpg = 6s
             
             %Ra2 = A2 = Aplus.jpg = 5.5s
             %Rb2 = B2 = Bplus.jpg = 5.5s
             

             if strcmp(curOrder,'1') || strcmp(curOrder,'3')
                 SS='A'; S='B';
             elseif strcmp(curOrder,'2') || strcmp(curOrder,'4')
                 SS='B'; S='A'; 
             else
                 fprintf('\nError for %s - current order unknown.\n',cbusub)
                 return
             end
             
 
                 % Reset onset vars.
                 A1_onsets=[]; A2_onsets=[]; A3_onsets=[];
                 B1_onsets=[]; B2_onsets=[]; B3_onsets=[];
                 A1_durs=[];   A2_durs=[];   A3_durs=[];
                 B1_durs=[];   B2_durs=[];   B3_durs=[];
                 
                 %A1 and B1 in all 3 runs of task.
                 for llA1=1:length(A1rows)
                     A1_onsets(llA1)=[(str2num(char(subsData.A1_OnsetTime(A1rows(llA1))))-taskOnset)/1000];
                     if ~isempty(strmatch('A1_Duration',subsData.Properties.VariableNames))
                         A1_durs(llA1)=[(str2num(char(subsData.A1_Duration(A1rows(llA1)))))/1000];
                     else
                         A1_durs(llA1)=6;
                         if llA1==1
                             fprintf('\nNo A1 duration for %s - %s.\n',cbusub,curRunName)
                         end
                     end
                 end
                 for llB1=1:length(B1rows)
                     B1_onsets(llB1)=[(str2num(char(subsData.B1_OnsetTime(B1rows(llB1))))-taskOnset)/1000];
                      if ~isempty(strmatch('B1_Duration',subsData.Properties.VariableNames))
                         B1_durs(llB1)=[(str2num(char(subsData.B1_Duration(B1rows(llB1)))))/1000];
                      else
                         B1_durs(llB1)=6;
                         if llB1==1
                             fprintf('\nNo B1 duration for %s - %s.\n',cbusub,curRunName)
                         end
                     end
                 end
                 
                 %A2 and B2 in acqui only.
                 if strcmp(runIDs(rr),'BOLD_Acquisition')
                     for llA2=1:length(A2rows)
                         A2_onsets(llA2)=[(str2num(char(subsData.A2_OnsetTime(A2rows(llA2))))-taskOnset)/1000];
                         if ~isempty(strmatch('A2_Duration',subsData.Properties.VariableNames))
                             A2_durs(llA2)=[(str2num(char(subsData.A2_Duration(A2rows(llA2)))))/1000]+.5;
                         else
                             A2_durs(llA2)=6;
                             if llA2==1
                                 fprintf('\nNo A2 duration for %s - %s.\n',cbusub,curRunName)
                             end
                         end
                     end
                     
                     for llB2=1:length(B2rows)
                         B2_onsets(llB2)=[(str2num(char(subsData.B2_OnsetTime(B2rows(llB2))))-taskOnset)/1000];
                         if ~isempty(strmatch('B2_Duration',subsData.Properties.VariableNames))
                             B2_durs(llB2)=[(str2num(char(subsData.B2_Duration(B2rows(llB2)))))/1000]+.5;
                         else
                             B2_durs(llB2)=6;
                             if llB2==1
                                 fprintf('\nNo B2 duration for %s - %s.\n',cbusub,curRunName)
                             end
                         end
                     end
                     
                 end
                 
                 %A3 and B3 in habit and ext only.
                 if strcmp(runIDs(rr),'BOLD_Habituation') || strcmp(runIDs(rr),'BOLD_Extinction')
                     for llA3=1:length(A3rows)
                         A3_onsets(llA3)=[(str2num(char(subsData.A3_OnsetTime(A3rows(llA3))))-taskOnset)/1000];
                         if ~isempty(strmatch('A3_Duration',subsData.Properties.VariableNames))
                             A3_durs(llA3)=[(str2num(char(subsData.A3_Duration(A3rows(llA3)))))/1000];
                         else
                             A3_durs(llA3)=6;
                             if llA3==1
                                 fprintf('\nNo A3 duration for %s - %s.\n',cbusub,curRunName)
                             end
                         end
                     end
                     for llB3=1:length(B3rows)
                         B3_onsets(llB3)=[(str2num(char(subsData.B3_OnsetTime(B3rows(llB3))))-taskOnset)/1000];
                         if ~isempty(strmatch('B3_Duration',subsData.Properties.VariableNames))
                             B3_durs(llB3)=[(str2num(char(subsData.B3_Duration(B3rows(llB3)))))/1000];
                         else
                             B3_durs(llB3)=6;
                             if llB3==1
                                 fprintf('\nNo B3 duration for %s - %s.\n',cbusub,curRunName)
                             end
                         end
                     end
                     
                     SSplusOnsets=eval([SS '3_onsets']);
                     SSplusDurs=eval([SS '3_durs']);
                     SplusOnsets=eval([S '3_onsets']);
                     SplusDurs=eval([S '3_durs']);
                 else
                     SSplusOnsets=eval([SS '2_onsets']);
                     SSplusDurs=eval([SS '2_durs']);
                     SplusOnsets=eval([S '2_onsets']);
                     SplusDurs=eval([S '2_durs']);
                 end

                 SSminusOnsets=eval([SS '1_onsets']);
                 SSminusDurs=eval([SS '1_durs']);
                 SminusOnsets=eval([S '1_onsets']);
                 SminusDurs=eval([S '1_durs']);
                 
                 SSonsets=[SSplusOnsets SSminusOnsets];
                 Sonsets=[SplusOnsets SminusOnsets];
                 SSdurs=[SSplusDurs SSminusDurs];
                 Sdurs=[SplusDurs SminusDurs];
                 
                 savefile=[cbusub '_' curRunName '_onsets.mat'];
                 save(savefile,'SSplusOnsets','SSminusOnsets','SplusOnsets', ...
                 'SminusOnsets','SSplusDurs','SSminusDurs','SplusDurs', ...
                 'SminusDurs','SSonsets','Sonsets','SSdurs','Sdurs')
                 
         end
           
%% Determine differences in trigger time and first stim onset.
%  (Just for troubleshooting. Can comment out.)
%              if rr==1
%                  if ~isempty(strmatch('trigger1_RTTime',subsData.Properties.VariableNames))
%                     trigDiff=(str2num(char(subsData.B1_OnsetTime(1)))-str2num(char(subsData.trigger1_RTTime(1))))/1000;
%                     
%                     trigDiffs=[trigDiffs trigDiff];
%                  end
%              end
    end
end

% trigDiffs
warning('on','MATLAB:table:ModifiedVarnames')