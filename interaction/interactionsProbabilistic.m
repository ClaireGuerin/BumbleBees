% IMPORT DATA

[vidFile, pathname] = uigetfile('F:\Academia\BumbleBees2016\*.avi','Select video file');
cd(pathname)

probBody = csvread('probBodyDistribCenteredAll.csv');
meanSize = csvread('meanBeeSizeCentAll.csv');

bwProb = probBody > 0.2;
% imshow(bwProb)

region = regionprops(bwProb, 'Centroid','BoundingBox');
midLine = ceil(region.Centroid(2));
frontMap = zeros(size(probBody));
frontMap(1:midLine,:) = 1;
frontMap = logical(frontMap);

[colony, date, time, chbr] = strread(vidFile, '%s %s %s %s', 'delimiter','_');
chbr = cellstr(chbr{1}(1:end-4));
trackFile = strcat(vidFile,'_tracked.mat');
revFile = 'reverseFrames.xlsx';
[num, text, raw] = xlsread(revFile);
scalesFile = fopen(strcat(pathname,'vidScales.txt'));
C_text = textscan(scalesFile,'%s',2,'Delimiter',',');
C_data = textscan(scalesFile,'%s %n', 'Delimiter',',');
fclose(scalesFile);
fileS = C_data{1};
scales = C_data{2};
sizFile = strcat(colony{1},'sizes.csv');
% ovaFile = strcat(colony,'ovariescore.csv');

S = load(trackFile);
vid = VideoReader(vidFile);
%scores = csvread(strcat(pathname,ovaFile));
sizes = csvread(sizFile);

coordinates = S.trackingData;
allTags = S.taglist;
% scoredTags = sizes(1,:);
% coordPlain = coordinates(:,toKeep,:);

nFrames = vid.NumberOfFrames;
frameWidth = vid.Width;
frameHeight = vid.Height;
popSize = size(allTags,1);
%nScored = size(scores,1);

%% CHANGE ZEROs TO NaNs & CORRECT TRACKING DATA

coordinates(coordinates==0) = NaN;

for indiv = 1:popSize
    
    untracked = numel(find(isnan(coordinates(:,indiv,1))));
    
    if untracked > 0.9 * nFrames
        coordinates(:,indiv,:) = NaN;
    else
        position = [coordinates(:,indiv,1), coordinates(:,indiv,2)];
        firstPos = repmat(position(1,:), size(position,1), 1);
        movement = abs(position - firstPos);
        meanMovement = nanmean(movement);
        
        if meanMovement < 10
            coordinates(:,indiv,:) = NaN;
        end
    end
end

allVideos = raw(2:size(raw,1),1);
nChar = size(vidFile,2);
myVid = find(strcmp(allVideos,vidFile(1:end-4)));
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
maxFrameGap = fRate * maxSecondGap; % max time gap to interpolate (in nFrames)
tic
coordFixed = fixShortNanGaps(coordinates, maxFrameGap);
disp(['It took ', num2str(toc), ' seconds to fix gaps'])

%% SMOOTH AND ADAPT SIZE DISTRIBUTION WITH SCALE

vidScale = scales(strcmp(fileS, vidFile));
tagEdge = 0.11;
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

%%
% INTERACTIONS

interactions = nan(popSize, popSize, nFrames);
orientations = nan(popSize, popSize, nFrames);
visualize = 1;
cXs = coordFixed(:,:,1);
cYs = coordFixed(:,:,2);
fXs = coordFixed(:,:,3);
fYs = coordFixed(:,:,4);

% create a local cluster object
%pc = parcluster('local');

% explicitly set the JobStorageLocation to the temp directory that was
% created in your sbatch script
%pc.JobStorageLocation = strcat('/scratch/cguerin/', getenv('SLURM_JOB_ID'));
v = VideoWriter('interProbsVis.avi');
open(v)

samp = 1:(nFrames/100);
% i = 1;

tic
for frm = samp %1:(nFrames/100)
    %%
    % frm = 1;
    
    fImage = rgb2gray(read(vid,frm));
    
    tagPos = nan(4,popSize);
    
    tagPos(1,:) = cXs(frm,:);
    tagPos(2,:) = cYs(frm,:);
    tagPos(3,:) = fXs(frm,:);
    tagPos(4,:) = fYs(frm,:);
    
    [interactions(:,:,frm), orientations(:,:,frm)] = calculateInteractionsProbs(tagPos,beeShape,probBodyResized,frontMapResized,meanSize,fImage,visualize);
    
%     subplot(2,2,i);
%     plotInteraction(fImage, tagPos,bodyStack);
%     
    if visualize == 1
        hold off
    end
    
    writeVideo(v,getframe(gcf))
%     i = i+1;
    
    %clearvars fImage
end

close(v)

save(strcat(vidFile, '_interactions.mat'), 'interactions', 'allTags');
disp(['Interactions calculated and saved in ', num2str(toc), ' sec.'])