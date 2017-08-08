% IMPORT DATA
pathname = 'H:\Academia\BumbleBees2016\Behav_Ovaries\Behav\A\test\';
cd(pathname)

k = dir('*.avi');
vidNames = {k.name}';
nVid = size(vidNames,1);

revFile = 'reverseFrames.xlsx';
[num, text, raw] = xlsread(revFile);

scalesFile = fopen('vidScales.txt');
C_text = textscan(scalesFile,'%s',2,'Delimiter',',');
C_data = textscan(scalesFile,'%s %n', 'Delimiter',',');
fclose(scalesFile);
fileS = C_data{1};
scales = C_data{2};

parpool(4)

for video = 1:nVid
    vidFile = vidNames{video};
    colony = vidFile(7);
%     ovaFile = strcat(colony,'ovariescore.csv');
    trackFile  = strcat(vidFile,'_tracked.mat');
    sizFile = strcat('Colony',colony,'sizes.csv');

    S = load(strcat(pathname,trackFile));
    vid = VideoReader(strcat(pathname,vidFile));
%     scores = csvread(strcat(pathname,ovaFile));
    sizes = csvread(strcat(pathname,sizFile));
    
    coordinates = S.trackingData;
    allTags = S.taglist;
%     scoredTags = sizes(1,:);
    % coordPlain = coordinates(:,toKeep,:);
    
    nFrames = vid.NumberOfFrames;
    frameWidth = vid.Width;
    frameHeight = vid.Height;
    % nPixels = frameWidth*frameHeight;
    popSize = size(allTags,1);
%     nScored = size(scores,1);
    
    %% CHANGE ZEROs TO NaNs & CORRECT TRACKING DATA
    
    coordinates(coordinates==0) = NaN;
    
    for indiv = 1:popSize
        
        untracked = numel(find(isnan(coordinates(:,indiv,1))));
        
        if untracked > 0.1 * nFrames
            coordinates(:,indiv,:) = NaN;
        else 
            position = [coordinates(:,indiv,1), coordinates(:,indiv,2)];
            firstPos = position(1,:);
            movement = abs(position - firstPos);
            meanMovement = nanmean(movement);
            
            if meanMovement < 1
                coordinates(:,indiv,:) = NaN;
            end
        end
            
    end
    
    allVideos = raw(2:size(raw,1),1);
    % nChar = size(vidFile,2);
    myVid = find(all(char(allVideos) == vidFile(1:end-4),2));
    shiftFrame = cell2mat(raw(myVid+1,2));
    shiftLine = cell2mat(raw(myVid+1,3));
    delete = cell2mat(raw(myVid+1,4));
    
    coordCorr = coordinates;
    centerXs = coordCorr(:,:,1);
    frontXs = coordCorr(:,:,3);
    
    centerXs((shiftFrame+1):end,:) = mod(centerXs((shiftFrame+1):end,:)+ frameWidth - shiftLine,frameWidth);
    frontXs((shiftFrame+1):end,:) = mod(frontXs((shiftFrame+1):end,:)+ frameWidth - shiftLine,frameWidth);
    
    coordCorr(:,:,1) = centerXs;
    coordCorr(:,:,3) = frontXs;
    
    %% FILL-IN SMALL TRACKING GAPS
    
    fRate = nFrames/(60*60); % number of frames per second (1H recording)
    maxSecondGap = 2; % max time gap to interpolate (in seconds)
    maxFrameGap = fRate * maxSecondGap; % max time gap to interpolate (in nOfFrames)
    % tic
    coordFixed = fixShortNanGaps(coordCorr, maxFrameGap);
    % toc
    
    %% ADAPT SIZES WITH SCALE
    
    vidScale = scales(strcmp(fileS, vidFile));
    tagEdge = 0.11;
    len = vidScale * sizes(2,:) / tagEdge;
    width = vidScale * sizes(3,:) / tagEdge;
    beeShape = [len;width];
    
    %%
    % INTERACTIONS QUANTIFICATION: PARPOOL
    
    interactions = nan(popSize, popSize, nFrames);
    sizeFactor = 1.2; % bee size + 20% of the body
    visualize = 0;
    
    cXs = coordFixed(:,:,1);
    cYs = coordFixed(:,:,2);
    fXs = coordFixed(:,:,3);
    fYs = coordFixed(:,:,4);
    
    % create a local cluster object
%     pc = parcluster('local');
    
    % explicitly set the JobStorageLocation to the temp directory that was
    % created in your sbatch script
%     pc.JobStorageLocation = strcat('/scratch/cguerin/', getenv('SLURM_JOB_ID'));
    
    % start the parallel pool with 12 workers
    
    
    % tic
    parfor frm = 1:nFrames
        %%
        %frm = 1;
        
        fImage = read(vid,frm);
        
        tagPos = nan(4,popSize);
        
        tagPos(1,:) = cXs(frm,:);
        tagPos(2,:) = cYs(frm,:);
        tagPos(3,:) = fXs(frm,:);
        tagPos(4,:) = fYs(frm,:);
        
        interactions(:,:,frm) = calculateInteractions(tagPos,beeShape,sizeFactor,fImage,visualize);
        
        if visualize == 1
            hold off
        end
        
        %writeVideo(v,getframe)
        
        %clearvars fImage
    end
    %toc
    
    save(strcat(vidFile, '_interactions.mat'), 'interactions', 'tagNumbers');
end