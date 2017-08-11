% IMPORT DATA
pathname = 'H:\Academia\BumbleBees2016\Behav_Ovaries\Behav\Odyssey\allFiles\track\'; 
cd(pathname)

k = dir('*tracked.mat'); % list tracking data in path
fileNames = {k.name}'; 
nFiles = size(fileNames,1); % number of files

revFile = 'reverseFrames.xlsx'; % file for coordinates correction of corrupted vids
[num, text, raw] = xlsread(revFile);

scalesFile = fopen(strcat(pathname,'vidScales.txt')); % file of video scaling
C_text = textscan(scalesFile,'%s',2,'Delimiter',',');
C_data = textscan(scalesFile,'%s %n', 'Delimiter',',');
fclose(scalesFile);
fileS = C_data{1};
scales = C_data{2};

% FOR PROBA

probBody = csvread('probBodyDistribCenteredAll.csv');
meanSize = csvread('meanBeeSizeCentAll.csv');

bwProb = probBody > 0.2;

region = regionprops(bwProb, 'Centroid','BoundingBox');
midLine = ceil(region.Centroid(2));
frontMap = zeros(size(probBody));
frontMap(1:midLine,:) = 1;
frontMap = logical(frontMap);
%

% start the parallel pool with n workers
% pc = parcluster('local');
% pc.JobStorageLocation = strcat('/scratch/cguerin/', getenv('SLURM_JOB_ID'));
pc = parpool(4);
pc.IdleTimeout = 120;
%%
for file = 2
    
    trackFile = fileNames{file};    
	[colony, date, time, chbr, trim] = strread(trackFile, '%s %s %s %s %s', 'delimiter','_');
    chbr = cellstr(chbr{1}(1:end-4));
    S = load(strcat(pathname,trackFile));   
	sizesFile = strcat(colony,'sizes.csv');
	sizes = csvread(strcat(pathname,char(sizesFile)));
   
    coordinates = S.trackingData;
    taglist = S.taglist;
    popSize = size(taglist,1);
    
    nFrames = size(coordinates,1);
	if strcmp(chbr{1},'FC')
		frameWidth = 2448;
		frameHeight = 2048;
	else
		frameWidth = 1288;
		frameHeight = 964;		
	end
    
    % CHANGE ZEROs TO NaNs & CORRECT TRACKING DATA
    
    coordinates(coordinates==0) = NaN;
	
	for indiv = 1:popSize
    
		untracked = numel(find(isnan(coordinates(:,indiv,1))));
    
		if untracked > 0.99 * nFrames
			coordinates(:,indiv,:) = NaN; % bee wasn't tracked on this frame
		else
			position = [coordinates(:,indiv,1), coordinates(:,indiv,2)];
			firstPos = repmat(position(1,:), size(position,1), 1);
			movement = abs(position - firstPos);
			meanMovement = nanmean(movement);
        
			if meanMovement < 6
				coordinates(:,indiv,:) = firstPos; % bee didn't actually move
			end
		end
	end

    allVideos = raw(2:size(raw,1),1);
    %nChar = size(vidFile,2);
    myVid = find(all(char(allVideos) == trackFile(1:end-16),2));
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
    
    % FILL-IN SMALL TRACKING GAPS
    
    fRate = nFrames/(60*60); % number of frames per second (1H recording)
    maxSecondGap = 2; % max time gap to interpolate (in seconds)
    maxFrameGap = fRate * maxSecondGap; % max time gap to interpolate (in nFrames)
    coordFixed = fixShortNanGaps(coordCorr, maxFrameGap);

    
    % ADAPT SIZES WITH SCALE
    
    vidScale = scales(strcmp(fileS, trackFile(1:end-12)));
    tagEdge = 0.11;
    len = vidScale * sizes(2,:) / tagEdge;
    width = vidScale * sizes(3,:) / tagEdge;
    beeShape = [len;width];
	
	% FOR PROBABILITY CALCULATION
	imScale = (meanSize(1,1) * vidScale) / (tagEdge * meanSize(2,1)); % prob distri: 2*2 inches
	meanSize(3,:) = meanSize(1,:) * vidScale / tagEdge;

	stepSize = floor(size(probBody)/45);
	probBody(isnan(probBody)) = 0;
	probBodySmooth = smooth2a(probBody, stepSize(1), stepSize(2));

	probBodyResized = imresize(probBodySmooth, imScale, 'method', 'nearest');
	frontMapResized = imresize(probBodySmooth, imScale, 'method', 'nearest');

	len = vidScale * sizes(2,:) / tagEdge;
	width = vidScale * sizes(3,:) / tagEdge;
	beeShape = [len;width];
	%

    % INTERACTIONS QUANTIFICATION: PARPOOL
    
%     nFrames = ceil(nFrames/6);
    
    interEllipses = zeros(popSize, popSize, nFrames);
    sizeFactor = 1.2; % bee size + 20% of the body
	
	interProbabilities = nan(popSize, popSize, nFrames);
	orientations = nan(popSize, popSize, nFrames);
    
    cXs = coordFixed(:,:,1);
    cYs = coordFixed(:,:,2);
    fXs = coordFixed(:,:,3);
    fYs = coordFixed(:,:,4);
    %%
    parfor frm = 1:nFrames

        fImage = [frameWidth,frameHeight]; % image size
        tagPos = nan(4,popSize);
        tagPos(1,:) = cXs(frm,:);
        tagPos(2,:) = cYs(frm,:);
        tagPos(3,:) = fXs(frm,:);
        tagPos(4,:) = fYs(frm,:);
        
        interEllipses(:,:,frm) = interactionEllipses(tagPos,beeShape,sizeFactor,fImage);
		[interProbabilities(:,:,frm), orientations(:,:,frm)] = interactionProbas(tagPos,beeShape,probBodyResized,frontMapResized,meanSize,fImage);
        
    end
    
    
    save(strcat(trackFile, '_interactions.mat'), 'interEllipses', 'interProbabilities', 'orientations','taglist');
end