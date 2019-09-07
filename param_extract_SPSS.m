function param_extract(varargin)
%=====================================================================%
% BETA EXTRACTION
%
%   This script allows the user to select a set of ROIs and a set of
%   contrast images, extracting parameter estimates.
%
%   Created 1/23/2015 by Jared Torre
%   Amended 7/30/2018 by CJL - added SPSS compatible output.
%
%=====================================================================%
%%%% BEGIN USER INPUT %%%%%

%Folder to save ROI output file in:
studyDIR='/u/project/sanscn/data/SAFETY/';
subIDpre='SAS';

taskID='HAE/Habit_Acqui';
analysisID='Habit_Aqui_firstlast4_hpf128ar1';

outFolder=fullfile(studyDIR,'_param_ests',taskID,analysisID);

%Optional inputs to make your life easier:
roiFolder='/u/project/sanscn/cjlescha/ROIs_2019';       %if leaving blank, leave as ''
conFolder=fullfile(studyDIR,'_contrasts',taskID,analysisID); %if leaving blank, leave as ''


%%%% END USER INPUT %%%%
%======================================================================


javaaddpath('/u/project/sanscn/cjlescha/toolboxes/xlwrite/jxl.jar')
javaaddpath('/u/project/sanscn/cjlescha/toolboxes/xlwrite/MXL.jar')
addpath(genpath('/u/project/sanscn/cjlescha/toolboxes/xlwrite'))

startDIR = pwd;
eval(sprintf('mkdir %s;',outFolder));

if ~isempty(roiFolder)
    cd(roiFolder);
end
roifiles = uigetvol('Select the ROI .img or .nii files...', 1);

if ~isempty(conFolder)
    cd(conFolder);
end
confiles = uigetvol('Select the contrast/beta .img or .nii files...', 1);

roifiles = cellstr(roifiles)';
confiles = cellstr(confiles)';

[roipath, roiname, roi_e] = cellfun(@(x) fileparts(x), roifiles, 'UniformOutput', false);
[conpath, conname, con_e] = cellfun(@(x) fileparts(x), confiles, 'UniformOutput', false);

nroi = length(roifiles);
ncon = length(confiles);

fprintf('\nThe following %s ROI(s) will be used:',num2str(nroi))
for i = 1:nroi
    fprintf('\n%s',roiname{i})
end

fprintf('\n\nThe following %s con/beta image(s) will be used:',num2str(ncon))
for i = 1:ncon
    fprintf('\n%s',fullfile(conpath{i},conname{i}))
end

% Convert the .hdr/.img ROIs to .nii

roi_hdrs_in_idx = strcmpi(roi_e,'.img');

if sum(roi_hdrs_in_idx) ~= 0
    fprintf('\n\nConverting %s ROI images to .nii...',num2str(sum(roi_hdrs_in_idx)))
    converted2nii_roifiles = cellfun(@(x,y) fullfile(x,strcat(y,'.nii')), roipath(roi_hdrs_in_idx), roiname(roi_hdrs_in_idx), 'UniformOutput', false);
    cellfun(@(x,y) img2nii(x,y), roifiles(roi_hdrs_in_idx),converted2nii_roifiles);
    roifiles(roi_hdrs_in_idx) = converted2nii_roifiles;
    fprintf('Done!');
end

% Create list of new resliced output filenames
rtmp_roifiles = cellfun(@(x,y) fullfile(x,strcat('rtmp_',y,'.nii')),roipath,roiname,'UniformOutput',false);

% Reslice ROI .nii images into sample contrast space
%==================================================%
fprintf('\n\nReslicing ROI image(s) into correct space for data image(s)...')
% flags = struct('interp', 1, ... % b-spline
%     'mask', 0, ...              % do not mask
%     'mean', 0, ...              % do not write mean image
%     'hold', -1, ...             % i don't think this is used anymore
%     'which', 1, ...             % reslice 2nd-nth only
%     'wrap', [0 0 0]', ...        % the default; don't know what this is
%     'prefix', 'rtmp_' ...
%     );
% P = [confiles(1);rtmp_roifiles];
% spm_reslice(P,flags)

% read in reference image
RefHead = spm_vol(confiles{1});
RefData = spm_read_vols(RefHead);
mat=RefHead.mat;
dim=RefHead.dim;

for r = 1:length(roifiles)
    % read in image to reslice
    SourceHead = spm_vol(roifiles{r});
    SourceData = spm_read_vols(SourceHead);
    
    % do the reslicing
    [x1,x2,x3] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
    d           = [1*[1 1 1]' [1 1 0]'];
    C = spm_bsplinc(SourceHead, d);
    v = zeros(dim);
    M = inv(SourceHead.mat)*mat; % M = inv(mat\SourceHead.mat) in spm_reslice.m
    y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
    y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
    y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
    out    = spm_bsplins(C, y1,y2,y3, d);
    
    %Revised by YAN Chao-Gan 121214. Apply a mask from the source image: don't extend values to outside brain.
    tiny = 5e-2; % From spm_vol_utils.c
    Mask = true(size(y1));
    Mask = Mask & (y1 >= (1-tiny) & y1 <= (SourceHead.dim(1)+tiny));
    Mask = Mask & (y2 >= (1-tiny) & y2 <= (SourceHead.dim(2)+tiny));
    Mask = Mask & (y3 >= (1-tiny) & y3 <= (SourceHead.dim(3)+tiny));
    
    out(~Mask) = 0;
    outmat = mat;
    
    OutHead=SourceHead;
    OutHead.mat      = mat;
    OutHead.dim(1:3) = dim;
%    [p n e] = fileparts(SourceHead.fname);
%    newname = sprintf('%s_%dx%dx%d%s',n,dim,e);
%    OutHead.fname = [p filesep newname];
    OutHead.fname = rtmp_roifiles{r};
    spm_write_vol(OutHead,out);
end
fprintf('Done!')

%Splits conname into conname & sub ID
splitCons=split(conname,['_' subIDpre]);      %[TOKEN,REMAIN] = strtok(conname,subIDpre);
for xx=1:length(splitCons)
    splitCons(xx,2)=strcat(subIDpre,char(splitCons(xx,2)));
end

connameOnly=splitCons(1:end,1);
subIDonly=splitCons(1:end,2);

% Set up data matrix
data = cell(length(confiles)+1,length(roifiles)+4);
data(:,1) = ['Full File'; confiles];
data(:,2) = ['Con Name + Sub ID'; conname];
data(:,3) = ['Con Name'; cellstr(connameOnly)];
data(:,4) = ['Subject ID'; cellstr(subIDonly)];
data(1,5:end) = roiname;

% Extract the parameters
for r = 1:nroi
    fprintf('\nExtracting from %s...',roiname{r})
    roi = rtmp_roifiles{r};
    hdr = spm_vol(roi); img = spm_read_vols(hdr);
    roiIDX{r} = find(img);
    
    for c = 1:ncon
        fprintf('%s',num2str(ceil(c/ncon*100)))
        if c ~=ncon            
            for i = 1:length(num2str(ceil(c/ncon*100)))
                fprintf('\b')
            end
        end
        conimg = confiles{c};
        hdr = spm_vol(conimg); img = spm_read_vols(hdr);
        data{c+1,r+4} = nanmean(img(roiIDX{r}));
    end
end
fprintf('\nAll done.')

a=datestr(clock,31);   % returns date string of the form 'YYYY-MM-DD HH:MM:SS' e.g., 2006-12-27 15:03:37
time_stamp = [a(6:7) a(9:10) a(3:4) '_' a(12:13) a(15:16)];   % timestamp is a function name, hence the _ in time_stamp
xlwrite(strcat(outFolder,filesep,analysisID,'_roidata_',time_stamp,'.xls'),data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Reorient for SPSS compat.   %%%%%%
%%%%%   Added by CJL, 7/30/2018   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subsInData=unique(data(2:end,4));
consInData=unique(data(2:end,3));
ncon=length(consInData);
subsRows={};
for jj=1:length(subsInData)
    subsRows{jj,1}=subsInData(jj);
    subsRows{jj,2}=find(ismember(data(:,4),subsInData(jj)));
end

%Set up SPSS variable names
conRow = ['.';sort(repmat(consInData,nroi,1))]';
roiRow = ['SubID';repmat(roiname,ncon,1)]';

for gg=1:length(conRow)
    SPSSvars(gg)=strcat(conRow(gg),'_',roiRow(gg));
    SPSSvars(gg)=regexprep(SPSSvars(gg),'-','_');    %replace dash with underscore (valid char for SPSS)
    SPSSvars(gg)=regexprep(SPSSvars(gg),' ','_');    %replace space with underscore (valid char for SPSS)
    SPSSvars(gg)=regexprep(SPSSvars(gg),'(','');     %delete parentheses
    SPSSvars(gg)=regexprep(SPSSvars(gg),')','');
    SPSSlabels(gg)=strcat(conRow(gg),' ',roiRow(gg));
    curVar=cell2char(SPSSvars(gg));
    if length(cell2char(SPSSvars(gg)))>64
        curVarShort=curVar(1:64);    %take only first 64 chars
    else
        curVarShort=curVar(1:end);
    end
    SPSSvarsShort(gg)=cellstr(curVarShort);
end



% Set up data matrix
data_SPSS = cell(length(subsInData)+4,((nroi*ncon)+1));
data_SPSS(:,1) = ['.';'SubID'; subsInData; '.';'.'];
data_SPSS(1,:) = ['.';sort(repmat(consInData,nroi,1))]';
data_SPSS(2,:) = ['SubID';repmat(roiname,ncon,1)]';
data_SPSS((length(subsInData)+3),:) = ['SPSS Label Names' SPSSlabels(2:end)];

SPSSvars=['SPSS Full Var Names' SPSSvars(2:end)]';
SPSSvarsShort=['SPSS Short Var Names' SPSSvarsShort(2:end)]';
SPSSlabels=['SPSS Label Names' SPSSlabels(2:end)]';
VarInfo=[SPSSlabels SPSSvars SPSSvarsShort];


%Reorient data from orig. PE output variable 'data'.
for kk=1:length(subsInData) %for each sub
    curRows=subsRows{kk,2};
    
    for ll=1:length(curRows) %for each row (contrast)
        curData=data(curRows(ll),:);
        curPEs=curData(1,(length(curData)-(nroi-1)):end);    
        curCon=curData(1,3);
        curConIndex=find(strcmp(curCon,data_SPSS(1,:)));
        
        for mm=1:nroi %for each ROI of that contrast
            if ~cellfun(@isempty,curPEs(1)) %if PE is not empty
                data_SPSS(kk+2,curConIndex(1))=curPEs(mm);
            else
                data_SPSS(kk+2,curConIndex(1))={'.'};     %%%% How to register non-existant PEs
            end
            curConIndex=curConIndex+1;
        end
    end
end

%determine # of columns exceeds max for excel (256), accounting for extra
%sub Column in each XL file.
numXLfiles=ceil(length(data_SPSS)/256);
if numXLfiles*256<length(data_SPSS)+numXLfiles-1
    numXLfiles=numXLfiles+1;
end

data_SPSSxl={};
curColStart=2;
curColEnd=255;

save('roi_workspace')
for yy=1:numXLfiles
    
    if curColEnd>size(data_SPSS,2)
        curColEnd=size(data_SPSS,2);
    end
  
    yyData={data_SPSS(:,1) data_SPSS(:,curColStart:curColEnd)};
    yyData=[yyData{:}];   
    data_SPSSxl{yy}=yyData;
    SPSS_xlFilename=strcat(outFolder,filesep,analysisID,'_roidata_SPSS_',time_stamp,'.xls');
    sheet_name=strcat('Sheet',string(yy),'of',string(numXLfiles));
    xlwrite(char(SPSS_xlFilename),yyData,char(sheet_name));
    clear yyData
    curColStart=curColEnd+1;
    curColEnd=curColEnd+254;
end

xlwrite(char(SPSS_xlFilename),VarInfo,'VariableInfo');

filename=strcat(outFolder,filesep,analysisID,'_roidata_',time_stamp,'.mat');
save(filename,'data','data_SPSS','data_SPSSxl','roiname','conname', 'connameOnly', 'nroi','ncon','VarInfo');

fprintf('\nROI data files output in %s\n',strcat(outFolder))
cd(outFolder);


%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% End CJL amendment %%%
%%%%%%%%%%%%%%%%%%%%%%%%%% 


% Cleanup the resliced and converted ROI files
if exist('converted2nii_roifiles')
    cellfun(@(x) delete(x), converted2nii_roifiles)
end
cellfun(@(x) delete(x), rtmp_roifiles)

end

% - SUBFUNCTIONS
function vol = uigetvol(message, multitag)
% UIGETVOL Dialogue for selecting image volume file
%
%   USAGE: vol = uigetvol(message, multitag)
%       
%       message = to display to user
%       multitag = (default = 0) tag to allow selecting multiple images 
%
if nargin < 2, multitag = 0; end
if nargin < 1, message = 'Select Image File'; end
if ~multitag
    [imname, pname] = uigetfile({'*.img; *.nii', 'Image File'; '*.*', 'All Files (*.*)'}, message);
else
    [imname, pname] = uigetfile({'*.img; *.nii', 'Image File'; '*.*', 'All Files (*.*)'}, message, 'MultiSelect', 'on');
end
if isequal(imname,0) || isequal(pname,0)
    vol = [];
    return
else
    vol = fullfile(pname, strcat(imname, ',1'));
end
%if size(vol,1)==1, vol = vol'; end
end

function img2nii(imgfile,outputname)

V=spm_vol(imgfile);
ima=spm_read_vols(V);
V.fname=fullfile(outputname);
spm_write_vol(V,ima);

end

function y = nanmean(x,dim)
% FORMAT: Y = NANMEAN(X,DIM)
% 
%    Average or mean value ignoring NaNs
%
%    This function enhances the functionality of NANMEAN as distributed in
%    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
%    identical name).  
%
%    NANMEAN(X,DIM) calculates the mean along any dimension of the N-D
%    array X ignoring NaNs.  If DIM is omitted NANMEAN averages along the
%    first non-singleton dimension of X.
%
%    Similar replacements exist for NANSTD, NANMEDIAN, NANMIN, NANMAX, and
%    NANSUM which are all part of the NaN-suite.
%
%    See also MEAN

% -------------------------------------------------------------------------
%    author:      Jan Glï¿½scher
%    affiliation: Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
%    
%    $Revision: 1.1 $ $Date: 2004/07/15 22:42:13 $

if isempty(x)
	y = NaN;
	return
end

if nargin < 2
	dim = min(find(size(x)~=1));
	if isempty(dim)
		dim = 1;
	end
end

% Replace NaNs with zeros.
nans = isnan(x);
x(isnan(x)) = 0; 

% denominator
count = size(x,dim) - sum(nans,dim);

% Protect against a  all NaNs in one dimension
i = find(count==0);
count(i) = ones(size(i));

y = sum(x,dim)./count;
y(i) = i + NaN;



% $Id: nanmean.m,v 1.1 2004/07/15 22:42:13 glaescher Exp glaescher $
end


