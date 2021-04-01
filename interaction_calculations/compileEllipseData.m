%% Load data

pathname = strcat('H:\Academia\BumbleBees2016\Behav_Ovaries\Behav\');
cd(pathname)

videoFiles = dir('*interactions.mat');
videoList = sortrows({videoFiles.name}');
nVideos = size(videoList,1);

maxNIndivPerColony = 100;
maxNComb = size(nchoosek(1:maxNIndivPerColony,2),1);
nDays = 4;
nRows = nVideos * maxNComb;
dataToStore = {'id1';'length1';'width1';'ov1';'id2';'length2';'width2';'ov2';'colony';'chamber';'day';'time';'int.ellps'};
nColumns = size(dataToStore,1);

compileDat = nan(nRows,nColumns);
index = 1;
durVid = 1; % videos duration in hours -> for data normalization

%% Compile data

for vNum = 1:nVideos
    
    vName = videoList{vNum};
    [colonyName, recDate, recTime, chbr] = strread(vName, '%s %s %s %s', 'delimiter','_');
    colonyN = double(colonyName{1}(end)); % ASCII A = 65, B = 66 or C = 67
    chbr = double(chbr{1}(1));% ASCII N(est) = 78 or F(oraging) = 70
    recDate = str2double(recDate{1}(1:2)); % extract only the day (5th,6th,7th or 10th of Jan)
    timeSS = str2double(recTime{1}(5:6)) / (60*60);
    timeMM = str2double(recTime{1}(3:4)) / 60;
    timeHH = str2double(recTime{1}(1:2));
    timeDec = timeSS+timeMM+timeHH;
    
    S = load(vName);
	
    ellpsInter = S.interactions;
    
    scores = csvread(strcat(colonyName{1}(end), 'ovariescore.csv'));
    sizes = csvread(strcat(colonyName{1},'sizes.csv'));
    popSize = size(sizes,2);
    
    combos = nchoosek(1:popSize,2);
    nComb = size(combos,1);
    
    for i = 1:nComb
        pair = combos(i,:);
        fstB = pair(1);
        secB = pair(2);
        
        ind1 = sizes(1,fstB);
        ind2 = sizes(1,secB);
        
        scored1 = scores(:,1) == ind1;
        scored2 = scores(:,1) == ind2;
        
        inter1 = S.tagNumbers == ind1;
        inter2 = S.tagNumbers == ind2;
        
        interactionE = nanmean(nanmean(nanmean(ellpsInter(inter1,inter2,:),3))); % mean number of interactions
        
        compileDat(index,1) = sizes(1,fstB); % Individual tag number 1
        compileDat(index,2) = sizes(2,fstB); % Length 1
        compileDat(index,3) = sizes(3,fstB); % Width 1
        
        if any(scored1)
            compileDat(index,4) = scores(scored1,2); % Ovary score 1
        else
            compileDat(index,4) = NaN; % Ovary score 1
        end
        
        compileDat(index,5) = sizes(1,secB); % Individual tag number 2
        compileDat(index,6) = sizes(2,secB); % Length 2
        compileDat(index,7) = sizes(3,secB); % Width 2
        
        if any(scored2)
            compileDat(index,8) = scores(scored2,2); % Ovary score 2
        else
            compileDat(index,8) = NaN; % Ovary score 2
        end
        
        compileDat(index,9) = colonyN; % Colony
        compileDat(index,10) = chbr; % Chamber
        compileDat(index,11) = recDate; % Day
        compileDat(index,12) = timeDec; % Time of the day
        compileDat(index,13) = interactionE; % Social interactions rate
        
        index = index + 1;
        
    end
end

out = compileDat(any(compileDat,2),:);
csvwrite([pathname,'compiledEllipseDataForR.csv'],out)

fileID = fopen([pathname,'compiledEllipseDataForR_headers.txt'],'w');
formatSpec = '%s\r\n';
[nrows,ncols] = size(dataToStore);
for row = 1:nrows
    fprintf(fileID,formatSpec,dataToStore{row,:});
end