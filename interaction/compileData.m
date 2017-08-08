%% Load data

pathname = strcat('H:\Academia\BumbleBees2016\Behav_Ovaries\Behav\Odyssey\allFiles\');
cd(pathname)
% ovaFiles = dir('*ovariescore.csv');
% ovaList = {ovaFiles.name}';
% sizeFiles = dir('*sizes.csv');
% sizeList = {sizeFiles.name}';

methodName = {'ellipses';'probabilities'};

videoFiles = dir('*.avi');
videoList = sortrows({videoFiles.name}');
nVideos = size(videoList,1);

maxNIndivPerColony = 100;
maxNComb = size(nchoosek(1:maxNIndivPerColony,2),1);
nDays = 4;
nRows = nVideos * maxNComb;
dataToStore = {'id1';'length1';'width1';'ov1';'id2';'length2';'width2';'ov2';'colony';'chamber';'day';'time';'int.ellps';'int.proba'};
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
    
    ellpsS = load([vName,'_inter',methodName{1},'.mat']);
    probaS = load([vName,'_inter',methodName{2},'.mat']);
    
%     % ensure the tag lists in both files are in the same order
    assert(all(ellpsS.allTags == probaS.allTags),'tags not in the same order!')
    
%     if any(ellpsS.allTags ~= probaS.allTags)
%         
%         ellpsTagsInd = 1:numel(ellpsS.allTags);
%         toSort = [probaS.allTags';ellpsTagsInd].';
%         sortedByProba = sortrows(toSort);
%         tagInd = sortedByProba(:,2);
%     else
%         tagInd = 1:numel(ellpsS.allTags);
%     end

    ellpsInter = ellpsS.interactions;
    probaInter = probaS.interactions;
    
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
        
        interE1 = ellpsS.allTags == ind1;
        interE2 = ellpsS.allTags == ind2;
        
        interP1 = probaS.allTags == ind1;
        interP2 = probaS.allTags == ind2;
        
        interactionE = nanmean(nanmean(nanmean(ellpsInter(interE1,interE2,:),3))); % mean number of interactions
        interactionP = nanmean(nanmean(nanmean(probaInter(interP1,interP2,:),3))); % mean probability of interaction
        
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
        compileDat(index,14) = interactionP; % Proba to interact
        
        index = index + 1;
        
    end
end

out = compileDat(any(compileDat,2),:);
csvwrite([pathname,'compiledDataForR.csv'],out)

fileID = fopen([pathname,'compiledDataForR_headers.txt'],'w');
formatSpec = '%s\r\n';
[nrows,ncols] = size(dataToStore);
for row = 1:nrows
    fprintf(fileID,formatSpec,dataToStore{row,:});
end
%%
% for i = 1:popSize
%
%     n = 0;
%
%     for j = 1:popSize
%         n = n + nansum(S.interactions(i,j,:));
%     end
%
%     if isnan(n)
%         n = 0;
%     end
%
%     allData(index,1) = scores(i,1); % Individual tag numbers
%     allData(index,2) = colony; % Colony
%     allData(index,3) = chamber; % Chamber
%     allData(index,4) = day; % Day
%     allData(index,5) = timeDec; % Time of the day
%     allData(index,6) = sizes(2,sizes(1,:) == scores(i,1)); % Length
%     allData(index,7) = sizes(3,sizes(1,:) == scores(i,1)); % Width
%     allData(index,8) = n / nFrames; % Social interactions rate
%     allData(index,9) = n; % Social interactions counts
%
%     allData(index,10) = scores(i,2); % Ovary scores
%     index = index + 1;
% end
% end
%
% end
%
% csvwrite('ovExp.csv',allData)

%%
% Example of Multinomial Regression for Ordinal Responses
% (https://www.mathworks.com/help/stats/mnrfit.html#btpyj65)

load carbig
X = [Acceleration Displacement Horsepower Weight];

[B,dev,stats] = mnrfit(X,miles,'model','ordinal');
% first 3 elements of B = intercept terms
% last 4 elements of B = coefficents of the covariates, assumed common
% across all categories
% i.e.  parallel regression (proportional odds model) = different intercept
% but common slopes among categories

[B(1:3)'; repmat(B(4:end),1,3)];
%For example, the coefficient estimate of 0.1048 indicates that a unit
% change in acceleration would impact the odds of the mpg of a car being
% less than or equal to 19 versus more than 19, or being less than or equal
% to 29 versus greater than 29, or being less than or equal to 39 versus
% greater than 39, by a factor of exp(0.01048) given all else is equal.

stats.p % Assess the significance of the coefficients.
% The  $p$-values of 0.035, 0.0000, and 0.0118 for engine displacement,
% horsepower, and weight of a car, respectively, indicate that these
% factors are significant on the odds of mpg of a car being less than or
% equal to a certain value versus being greater than that value.