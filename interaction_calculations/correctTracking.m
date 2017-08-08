% written by James Crall to work on one tracking file
% modified by Claire Guerin to work on all tracking files in a folder

% IMPORT FILES
path = = strcat('COPY AND PASTE YOUR PATHNAME HERE\');
cd(path)

trackFiles = dir('*.avi');
trackNames = {trackFiles.name}'; 
nFiles = size(trackNames,1); % number of videos

for fff = 1:nFiles
	
	matFile = trackNames{fff};
	S = load(strcat(path,matFile));
	vidFile = strrep(matFile,'_tracked.mat','');
	vid = VideoReader(strcat(path,vidFile));
	colony = vidFile(1:7);
	sizes = csvread([colony,'sizes.csv']);
	[num, text, raw] = xlsread('reverseFrames.xlsx');

	coordinates = S.trackingData;
	taglist = S.taglist;
	toKeep = ismember(allTags, sizes(1,:));
	coordPlain = coordinates(:,toKeep,:); % only keep tracking data of bees for which we have ovary scores (and sizes)

	videos = raw(2:size(raw,1),1);
	nChar = size(vidFile,2);
	myVid = find(all(char(videos) == vidFile(1:nChar-4),2));

	coordinates(coordinates == 0) = NaN; %Remove zeros
	movWidth = vid.Width;

	frame = cell2mat(raw(myVid+1,2));
	line = cell2mat(raw(myVid+1,3));
	delete = cell2mat(raw(myVid+1,4));

	coordCorr = coordinates;
	xCenter = coordCorr(:,:,1); %isolate x-coordinates
	xFront = coordCorr(:,:,3);

	%load back into coordCoor:
	coordCorr(:,:,1) = xCenter;
	coordCorr(:,:,3) = xFront;
	
	save(strcat(myVid, '_trackedcorr.mat'), 'coordCorr', 'taglist');

end