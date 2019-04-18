function processStudy(varargin);
%   processStudy is a MATLAB function designed to provide "hands-off"
%   preprocesxing to diffusion MRI images.
%
%   processStudy(studyDir) procprocesses the input directory 'studyDir' and
%   returns diffusion/kurtosis maps
%
%   Author: Siddhartha Dhiman
%   Date Created: 02/14/2019
%   Matlab 2018b
warning off;

%% Input Parsing
%   Set default values
defaultTBSS = 'off';

%   Initialize input parser
p = inputParser;
addOptional(p,'studyPath');
addParameter(p,'tbss',defaultTBSS);

%   Parse input
parse(p,width,varargin{:});
%% Conda Variables
envName = 'py36';

%% Define Paths
addpath(fullfile(pwd,'dependencies'));

%% Locate Study Directory
%   Opens GUI to locate directory
studyPath = uigetdir(pwd,'Select subject directory');

%% Check and Decompress Files if Present
studyDirFolders = dir(studyPath);
compressedFiles = vertcat(dir(fullfile(studyPath,'**/*.zip')),...
    dir(fullfile(studyPath,'**/*.tar')),...
    dir(fullfile(studyPath,'**/*.gz')));

if length(compressedFiles) >= 1
    parfor i = 1:length(compressedFiles)
        %   Check file extension
        [~,name,ext] = fileparts(fullfile(compressedFiles(i).folder,...
            compressedFiles(i).name));
        mkdir(fullfile(compressedFiles(i).folder,name));
        try
            if ext == '.zip';
                unzip(fullfile(compressedFiles(i).folder,...
                    compressedFiles(i).name),...
                    fullfile(compressedFiles(i).folder));
            elseif ext == '.tar'
                untar(fullfile(compressedFiles(i).folder,...
                    compressedFiles(i).name),...
                    fullfile(compressedFiles(i).folder));
            elseif ext == '.gz'
                gunzip(fullfile(compressedFiles(i).folder,...
                    compressedFiles(i).name),...
                    fullfile(compressedFiles(i).folder));
            else
            end
        catch
            fprintf('Corruption detected: skipping uncompression of %s',name);
            continue
        end
    end
end

%   Snap picture of original directory
studyDirFolders = dir(studyPath);

% Clean it up
for i = 1:length(studyDirFolders)
    rmPattern = [".","..",".DS_Store"];
    %   Check for '.'
    if any(contains(studyDirFolders(i).name,rmPattern));
        rmIdx(i) = 1;
    else
        %   If nothing found, don't mark for deletion
        rmIdx(i) = 0;
    end
end
studyDirFolders(rmIdx ~= 0) = [];   %   Apply deletion filter

%   Assign Output
[fp,~,~] = fileparts(studyPath);
outPath = fullfile(fp,'Processing');
if exist(outPath) == 7;
    rmdir(outPath,'s');
    mkdir(outPath);
else
    mkdir(outPath);
end

%% Sort Dicoms
studyDir = vertcat(dir(fullfile(studyPath,['**' filesep '*'])),...
    dir(fullfile(studyPath,['**' filesep '*.dcm']))); %   Recursive directory listing
rmPattern = {'.','.DS_Store'};   %   Remove files beginning with

% Clean-up Dicom Directory Listing
rmIdx = zeros(1,length(studyDir));
for i = 1:length(studyDir)
    %  Check for directories
    if studyDir(i).isdir == 1
        rmIdx(i) = 1;
        
        %   Check for files starting with'.'
    elseif any(startsWith(studyDir(i).name,'.'));
        rmIdx(i) = 1;
        
    else
        %   If nothing found, don't mark for deletion
        rmIdx(i) = 0;
    end
end
studyDir(rmIdx ~= 0) = [];   %   Apply deletion filter
nFiles = length(studyDir);
fprintf('Found %d possible images in %s\n',nFiles,studyPath);
%   Initialize Parallel Data Queue
parQ = parallel.pool.DataQueue;
%   Initialize progress waitbar
parWaitBar = waitbar(0,'Initializing sorting algorithm...',...
    'Name','Sorting dicoms');
%   After receiving new data, update_progress() will be called
fprintf('Sorting dicoms...');
afterEach(parQ,@parProgress);
n_completed = 0;

parfor i = 1:nFiles
    try
        tmp = dicominfo(fullfile(studyDir(i).folder,studyDir(i).name));
    catch
        continue
    end
    %   Replace all '.' in protocol names with '_'
    tmp.ProtocolName = strrep(tmp.ProtocolName,'.','_');
    %   Replace spaces
    tmp.ProtocolName = strrep(tmp.ProtocolName,' ', '_');
    
    %   Define copy/move folder
    defineFolder = fullfile(outPath,tmp.PatientID,...
        sprintf('%s_%d',tmp.ProtocolName,tmp.SeriesNumber));
    
    if ~exist(defineFolder,'dir')
        mkdir(defineFolder);
    end
    
    if ~contains(studyDir(i).name,'.dcm')
        newName = [studyDir(i).name '.dcm'];
    else
        newName = studyDir(i).name;
    end
    copyfile(tmp.Filename,fullfile(defineFolder,newName));
    send(parQ,i);
end
delete(parWaitBar);
fprintf('done\n')

%% Create Subject List
clear rmIdx
subList = dir(outPath);
rmPattern = [".","..",".DS_Store"];
for i = 1:length(subList)
    %   Check for non-directories
    if subList(i).isdir ~= 1
        rmIdx(i) = 1;
        
        %   Check for '.'
    elseif any(contains(subList(i).name,rmPattern));
        rmIdx(i) = 1;
        
    else
        %   If nothing found, don't mark for deletion
        rmIdx(i) = 0;
    end
end
subList(rmIdx ~= 0) = [];   %   Apply deletion filter
subList = {subList(:).name};

%% Unique Protocol List
seriesNames = dir(fullfile(outPath,subList{1}));
if length(seriesNames) > 1
    for i = 2:length(subList)
        tmp = dir(fullfile(outPath,subList{i}));
        seriesNames = vertcat(seriesNames,tmp);
    end
    for i = 1:length(seriesNames)
        seriesList{i} = seriesNames(i).name;
    end
else
    seriesList = seriesNames.name;
end

%   Clean up to remove '.' and '..'
seriesList(strcmp(seriesList(:),'.')) = [];
seriesList(strcmp(seriesList(:),'..')) = [];

%   Remove sequence number
%   Find positions of '_' in strings
for i = 1:length(seriesList)
    tmpIdx = strfind(seriesList{i},'_');
    seriesList(i) = extractBetween(seriesList{i},1,tmpIdx(end)-1);
end

%   Create unique list and count occurence
[seriesUniList,~,seriesUniIdx] = unique(seriesList);
seriesUniCnt = histc(seriesUniIdx,1:numel(seriesUniList));

%   Locate unique protcols which have the same numer of occurences as number
%   of subjects
cmnSeriesIdx = find(seriesUniCnt == numel(subList));

%   Refine unique list down to protcols which only exists in all subjects
seriesUniList = seriesUniList(cmnSeriesIdx);
fprintf('Found the following commonly occuring imaging protocols:\n');

[seriesIdx seriesSel] = listdlg('PromptString','Select relevant sequences:',...
    'SelectionMode','multiple',...
    'OKString','Run Designer',...
    'CancelString','Cancel',...
    'ListSize',[300 500],...
    'ListString',seriesUniList);

%% Create Designer Output Directory
desPath = fullfile(fp,'Output');
mkdir(desPath);
for i = 1:length(subList)
    mkdir(fullfile(desPath,subList{i}));
end

repInput = repmat('%s,',1,numel(seriesIdx));
% Remove comma at the end
repInput = repInput(1:end-1);

%% Run Designer
for i = 1:length(subList)
    subSeriesDir = dir(fullfile(outPath,subList{i}));
    for j = 1:length(subSeriesDir)
        subSeries{j} = subSeriesDir(j).name;
    end
    %   Clean up to remove '.' and '..'
    subSeries(strcmp(subSeries(:),'.')) = [];
    subSeries(strcmp(subSeries(:),'..')) = [];
    
    for j = 1:length(seriesIdx)
        %   Match modified strings for readibility with actual directory
        matchIdx = startsWith(subSeries,seriesUniList(seriesIdx(j)));
        %   If there are more than one strings, pick the one that most
        %   closely matches based on euclidean distance
        if nnz(matchIdx) > 1
            matchLoc = find(matchIdx);
            for k = 1:length(matchLoc)
                strDist(k) = sqrt(strlength(seriesUniList(seriesIdx(j)))...
                    ^2 + strlength(subSeries(matchLoc(k)))^2);
            end
            [~,matchIdx] = min(strDist);
            seriesUniPath(j) = fullfile(outPath,subList{i},...
                subSeries(matchLoc(matchIdx)));
        else
            seriesUniPath(j) = fullfile(outPath,subList{i},...
                subSeries(matchIdx));
        end
        
    end
    desInput = sprintf(repInput,seriesUniPath{:});
    desOutput = sprintf(fullfile(desPath,subList{i}));
    
    %     runDesigner = sprintf('./designer.sh %s %s',...
    %         desInput,desOutput);
    %     runDesigner = sprintf('source activate %s && python designer.py -denoise -extent 5,5,5 -degibbs -rician -mask -prealign -eddy -rpe_header -smooth 1.25 -DKIparams -DTIparams %s %s',...
    %         envName,desInput,desOutput);
    runDesigner = sprintf('source activate %s && python designer.py -denoise -extent 5,5,5 -degibbs -mask -rician -prealign -rpe_header -eddy -median -smooth 1.25 -DKIparams -DTIparams %s %s',...
        envName,desInput,desOutput);
    [s,t] = system(runDesigner,'-echo');
    
    %% Create SNR Plots
    mkdir(fullfile(desOutput,'QC'));
    
    %   Load all files
    bvalPath = dir(fullfile(desOutput,'*bval'));
    bvalPath = fullfile(bvalPath.folder,bvalPath.name);
    bval = round(load(bvalPath));
    
    %     bvecPath = dir(fullfile(desOutput,'*bvec'));
    %     bvecPath = fullfile(bvecPath.folder,bvecPath.name);
    
    procPath = dir(fullfile(desOutput,'dwi_designer.nii'));
    procPath = fullfile(procPath.folder,procPath.name);
    Iproc = niftiread(procPath);
    
    rawPath = dir(fullfile(desOutput,'dwi_raw*'));
    rawPath = fullfile(rawPath.folder,rawPath.name);
    Iraw = niftiread(rawPath);
    
    noisePath = dir(fullfile(desOutput,'fullnoisemap*'));
    noisePath = fullfile(noisePath.folder,noisePath.name);
    Inoise = niftiread(noisePath);
    Inoise(find(isnan(Inoise))) = 0;

    brainPath = dir(fullfile(desOutput,'brain_mask*'));
    brainPath = fullfile(brainPath.folder,brainPath.name);
    brainMask = logical(niftiread(brainPath));
    
    %     csfPath = dir(fullfile(desOutput,'CSFmask*'));
    %     csfPath = fullfile(csfPath.folder,csfPath.name);
    %     csfMask = ~logical(niftiread(csfPath));
    %
    %     histoMask = logical(brainMask .* csfMask);
    
    %   Create logical indexes of bvalues
    listBval = unique(bval);
    bIdx = zeros(length(listBval),length(bval));
    for j = 1:length(listBval)
        bIdx(j,:) = bval == listBval(j);
    end
    bIdx = logical(bIdx);
    
    %   Form index vector for mean 4D array
    nB0 = find(bIdx(1,:));
    meanIdx = horzcat(repmat(listBval(1),[1 numel(nB0)]),listBval(2:end));
    
    for j = 1:length(listBval)
        for k = 1:length(meanIdx)
            if k <= length(nB0);
                tmpProc = Iproc(:,:,:,nB0(k)) ./ Inoise;
                snrProc(:,k) = tmpProc(brainMask);
                
                tmpRaw = Iraw(:,:,:,nB0(k)) ./ Inoise;
                snrRaw(:,k) = tmpRaw(brainMask);
                
            elseif meanIdx(k) == listBval(j)
                tmpProc = mean(Iproc(:,:,:,bIdx(j,:)),4)./Inoise;
                snrProc(:,k) = tmpProc(brainMask);
                
                tmpRaw = mean(Iraw(:,:,:,bIdx(j,:)),4)./Inoise;
                snrRaw(:,k) = tmpRaw(brainMask);
            else
                continue;
            end
        end
    end
    
    %   Histogram counts --------------------------------------------------
    %   Generate subplot layout such that plots are in nPlots x 2 layout.
    %   The firsy column contains plots and the second contains a single
    %   legend. Plot layout sets the size of subplot and plotLocations
    %   determines where to plot them on a nPlots x 2 grid.
    plotLayout = [numel(listBval) 3];
    plotLocation = linspace(1,(numel(listBval)+4),numel(listBval));
    
    figure; fig = gcf;
    
    set(fig,'PaperUnits','inches','PaperPosition',[.25 .25 10 12],...
        'InvertHardcopy','off','Color','white','Visible','off');
    
    for j = 1:length(listBval)
        [Nproc Eproc] = histcounts(snrProc(:,meanIdx == listBval(j)),...
            'Normalization','Probability');
        [Nraw Eraw] = histcounts(snrRaw(:,meanIdx == listBval(j)),...
            'Normalization','Probability');
        
        Nproc = smooth(Nproc);
        Nraw = smooth(Nraw);
        
        %   Compute median value of each bin
        for k = 1:length(Nproc)
            Mproc(k) = median([Eproc(k) Eproc(k+1)]);
        end
        for k = 1:length(Nraw)
            Mraw(k) = median([Eraw(k) Eraw(k+1)]);
        end
        
        %   Plot graphs ---------------------------------------------------
        
        %   Plotting colors
        c1 = [86,187,131]/255;   % processed color
        c2 = [78,173,241]/255;   % raw color
        c3 = [235,235,235]/255;  % background color
        
        subplot(numel(listBval),3,[plotLocation(j) plotLocation(j)+1])
        hold on;
        pArea = area(Mproc,Nproc,...
            'EdgeColor',c1,'LineWidth',3,...
            'FaceColor',c1,'FaceAlpha',0.55);
        rArea = area(Mraw,Nraw,...
            'EdgeColor',c2,'LineWidth',3,...
            'FaceColor',c2,'FaceAlpha',0.55);
        hold off;
        xlabel('SNR');
        ylabel('%Voxels');
        
        if listBval(j) == 0 
            xlim([0 500]);
        elseif listBval(j) == 1000
            xlim([0 80]);
        elseif listBval(j) == 2000
            xlim([1 40]);
        else
            xlim([1 20]);
        end
        
        title(sprintf('B%d',listBval(j)));
        grid on; box on; legend off;
        %   Set axial colors
        ax = gca;
        set(ax,'Color',c3,...
            'GridColor','white','GridAlpha',1,'MinorGridAlpha',0.15,...
            'fontname','helvetica','FontWeight','bold','fontsize',14);
        clear Nproc Eproc Nraw Eraw Mproc Mraw
    end
    legPlot = subplot(numel(listBval),3,numel(listBval)*3)
    posLeg = get(legPlot,'Position');
    leg1 = legend([pArea,rArea],...
        {'After Preprocessing','Before Preprocessing'});
    set(leg1,'Position',posLeg);
    title(leg1,'Legend');
    axis(legPlot,'off');
    print(fullfile(desOutput,'QC','SNR_Plots'),'-dpng','-r800');
    clear tmpRaw tmpProc snrProc snrRaw
end

%% PARFOR Progress Calculation
    function parProgress(~)
        if ~exist('n_completed','var')
            n_completed = 0;
        else
            n_completed = n_completed + 1;
        end
        %   Calculate percentage
        parPercentage = n_completed/nFiles*100;
        %   Update waitbar
        waitbar(n_completed/nFiles,parWaitBar,...
            sprintf('%0.1f%% completed\nHang in there...',parPercentage));
    end
end

