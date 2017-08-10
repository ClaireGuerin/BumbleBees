% IMPORT DATA
pathname = '/n/regal/debivort_lab/claire/bees/track/track'; 
cd(pathname)

k = dir('*tracked.mat'); % list tracking data in path
fileNames = {k.name}'; 
nFiles = size(fileNames,1); % number of files

[revFile,fileS, scales] = getFiles(pathname);

pc = parcluster('local');
pc.JobStorageLocation = strcat('/scratch/cguerin/', getenv('SLURM_JOB_ID'));
parpool(pc, 32)

for file = 1:nFiles

	trackFile = fileNames{file};    
	[colony, date, time, chbr, trim] = strread(trackFile, '%s %s %s %s %s', 'delimiter','_');
    chbr = cellstr(chbr{1}(1:end-4));
    S = load(strcat(pathname,trackFile)); 

	[sizes, shiftFrame, shiftLine, delete] = getInfo(pathname,trackFile,colony)	
   
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
	
	coordCorr = trackingDataCorrection(coordinates, shiftFrame, shiftLine, frameWidth);
	
    % FILL-IN SMALL TRACKING GAPS
    
    fRate = nFrames/(60*60); % number of frames per second (1H recording)
    maxSecondGap = 2; % max time gap to interpolate (in seconds)
    maxFrameGap = fRate * maxSecondGap; % max time gap to interpolate (in nFrames)
    coordFixed = fixShortNanGaps(coordCorr, maxFrameGap);

	% ADAPT SIZES WITH SCALE
    
    vidScale = scales(strcmp(fileS, trackFile(1:end-12)));
    tagEdge = 0.11;