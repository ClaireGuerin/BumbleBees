%% FIND AVERAGE BEE

pathname = 'E:\Academia\BumbleBees2016\Vid_Behav_Ovaries\Odyssey';
cd(pathname)

k = dir('*sizes.csv');
sizFiles = {k.name}';
nCol = size(sizFiles,1);
nIndivPerCol = 20;

sizes = nan(nIndivPerCol*nCol, 4);
indices = [1:20;21:40;41:60];

for i = 1:nCol
    f = sizFiles{i};
    ind = indices(i,:);
    colSizes = csvread(f);
    sizes(ind,1) = i;
    sizes(ind,2:4) = colSizes';
end

meanLength = mean(sizes(:,3));
meanWidth = mean(sizes(:,4));

sizeDiff = mean([abs(sizes(:,3) - meanLength), abs(sizes(:,4) - meanWidth)], 2);

nBeeSamp = 5;
allMinDiff = nan(1,nBeeSamp);
allMinIndex = nan(1,nBeeSamp);
tmpSizeDiff = sizeDiff;

for bSamp = 1:nBeeSamp
    [minDiff, minIndex] = min(tmpSizeDiff);
    tmpSizeDiff(minIndex) = 1000;
    allMinDiff(1,bSamp) = minDiff;
    allMinIndex(1,bSamp) = minIndex;
end

averageB = sizes(allMinIndex,:);
beeTagNumber = averageB(:,2);
% disp(['Bee ',num2str(beeTagNumber),' in colony ',Bcolony])

disp(['On average, a bee is ', num2str(meanLength), ' inches long and ', num2str(meanWidth), ' inches wide.'])

% maxB = [max(sizes(:,3));max(sizes(:,4))];

%% CALCULATE DISTRIBUTION PROBABILITY

rng default % for reproducibility

tagScale = 0.11;
nSample = 20;
nBees = numel(beeTagNumber);

allBWFora = nan(500,500,8*nSample*nBees);
allBWNest = nan(500,500,8*nSample*nBees);

i = 0;

for Bs = 1:nBees
    
    bTag = beeTagNumber(Bs);
    
    sizName = sizFiles{averageB(Bs,1)};
    sizes = csvread(sizName);
    sizeB = sizes(2:3,sizes(1,:) == bTag);
    
    Bcolony = sizName(1);
    j = dir(strcat('Colony',Bcolony,'*.avi'));
    vidFiles = {j.name}';
    nVid = size(vidFiles, 1);
    
    scalesFile = fopen('vidScales.txt');
    C_text = textscan(scalesFile,'%s',2,'Delimiter',',');
    C_data = textscan(scalesFile,'%s %n', 'Delimiter',',');
    fclose(scalesFile);
    fileS = C_data{1};
    scales = C_data{2};

    for video = 1:nVid
        
        vidName = vidFiles{video};
        [colony, date, time, chbr] = strread(vidName, '%s %s %s %s', 'delimiter','_');
        chbr = string(chbr{1}(1:end-4));
        mov = VideoReader(vidName);
        nFrames = mov.NumberOfFrames;
        
        vidScale = scales(strcmp(fileS, vidName));
        sizeBPX = sizeB * vidScale / tagScale;
        
        trackName = strcat(vidName,'_tracked.mat');
        S = load(trackName);
        coordinates = S.trackingData;
        tags = S.taglist;
        bee = find(tags == bTag);
        [ZeroFrame, ZeroCoord] = find(coordinates(:,bee,:) == 0);
        [NanFrame, NanCoord] = find(isnan(coordinates(:,bee,:)));
        
        allFrames = 1:nFrames;
        detectFrames = [];
        
        for frame = allFrames
            
            if all(frame ~= [ZeroFrame;NanFrame])
                detectFrames = [detectFrames,frame];
            end
        end
        
        nSampleAdjusted = min(nSample,numel(detectFrames));
        sampleFrames = datasample(detectFrames, nSampleAdjusted, 'replace', false);
        
        for catchFrame = sampleFrames
            %
            %         catchFrame = sampleFrames(1);
            i = i + 1;
            
            im = read(mov,catchFrame);
            coord = coordinates(catchFrame,bee,:);
            
            xfront = coord(:,:,3);
            yfront = coord(:,:,4);
            xcenter = coord(:,:,1);
            ycenter = coord(:,:,2);
            
            %         slope = (yfront - ycenter) / (xfront - xcenter);
            %         angle = atan(slope);
            angle = atan2(yfront - ycenter,xfront - xcenter);
            angle = angle .* (angle >= 0) + (angle + 2 * pi) .* (angle < 0);
            
            pxLength = sizeB(1)*vidScale/tagScale;
            xbarycenter = xcenter + 1/4*pxLength*cos(angle+pi);
            ybarycenter = ycenter + 1/4*pxLength*sin(angle+pi);
            
            %         imshow(im)
            %         hold on
            %         plot(xbarycenter, ybarycenter, 'ro' )
            %         plot(xcenter,ycenter, 'bo')
            %
            %         xmin = coord(:,:,1) - 0.5*vidScale/tagScale;
            %         ymin = coord(:,:,2) - 0.5*vidScale/tagScale;
            imCrop = imcrop(im, [xbarycenter - 0.5*vidScale/tagScale, ybarycenter - 0.5*vidScale/tagScale, vidScale/tagScale, vidScale/tagScale]);
            
            %         angle = atan2(coord(:,:,4) - coord(:,:,2), coord(:,:,3) - coord(:,:,1));
            rotationAngleDegrees = rad2deg(pi/2 + angle);
            rotatedIm = imrotate(imCrop, rotationAngleDegrees);
            figure
            %         imshow(rotatedIm)
            
            BW = roipoly(rotatedIm);
            close
            
            if chbr == 'FC'
                allBWFora(1:size(rotatedIm,1),1:size(rotatedIm,2),i) = BW;
            elseif chbr == 'NC'
                allBWNest(1:size(rotatedIm,1),1:size(rotatedIm,2),i) = BW;
            end
            
            % figure(1)
            % subplot(2,2,1)
            % imshow(im)
            % hold on
            % plot(coord(:,:,1), coord(:,:,2), 'ro')
            % line([coord(:,:,1), coord(:,:,3)],[coord(:,:,2), coord(:,:,4)], 'LineWidth', 2)
            % hold off
            % subplot(2,2,2)
            % imshow(imCrop)
            % subplot(2,2,3)
            % imshow(rotatedIm)
            
        end
    end
end
%%

% save('distriProbaMats.mat', 'allBWFora', 'allBWNest', '-v7.3');

nEffFora = size(allBWFora,3) - size(find(isnan(allBWFora(1,1,:))),1); % 204
nEffNest = size(allBWNest,3) - size(find(isnan(allBWNest(1,1,:))),1); % 60

probFora = sum(allBWFora,3, 'omitnan')/nEffFora;
probNest = sum(allBWFora,3, 'omitnan')/nEffNest;
probTot = (sum(allBWFora,3, 'omitnan') + sum(allBWFora(1:size(allBWFora,1),1:size(allBWFora,2),:),3, 'omitnan')) / (nEffFora + nEffNest);

% probTot = csvread('probBodyDistrib.csv');

bwProb = probTot > 0.4;
imshow(bwProb)

region = regionprops(bwProb, 'centroid');
limits = ceil(region.Centroid * 2);
cropProbTot = probTot(1:limits(2), 1:limits(1));
marginTop = ceil((limits(2) - limits(1))/2);
marginBot = floor((limits(2) - limits(1))/2);
imshow(cropProbTot)
hold on
plot(region.Centroid(1),region.Centroid(2), 'ro')
% [x, y] = ginput(4);
% close ALL
hold off

% csvwrite('probBodyDistrib.csv', cropProbTot);
% 
% probLength = sqrt((x(1) - x(2))^2 + (y(1) - y(2))^2);
% probWidth = sqrt((x(3) - x(4))^2 + (y(3) - y(4))^2);
% meanSize = [meanLength, meanWidth; probLength, probWidth; region.Centroid(1), region.Centroid(2)];
% csvwrite('meanBeeSize.csv', meanSize)

probDensRows = mean(cropProbTot,2);
probDensCols = mean(cropProbTot,1);

figure

subplot(1,2,1)
hold on
plot(1:size(cropProbTot,1),probDensRows, '-b', 'userdata', 'Image Height')
plot(1:size(cropProbTot,2),probDensCols, '--r', 'userdata', 'Image Length')
hold off
legend(get(gca, 'children'), get(get(gca, 'children'), 'userdata')) 
title('Observed probabilities')

subplot(1,2,2)
hold on
[fRows,xiRows] = ksdensity(probDensRows,1:size(cropProbTot,1));
plot(xiRows, fRows, '-b', 'userdata', 'Image Height')
[fCols,xiCols] = ksdensity(probDensCols,1:size(cropProbTot,2));
plot(xiCols, fCols, '--r', 'userdata', 'Image Length')
hold off
legend(get(gca, 'children'), get(get(gca, 'children'), 'userdata')) 
title('Kernel density function')

countsTot = sum(allBWFora,3, 'omitnan') + sum(allBWFora(1:size(allBWFora,1),1:size(allBWFora,2),:),3, 'omitnan');
countsTotCrop = countsTot(1:limits(2), 1:limits(1));
countsCol = sum(countsTotCrop,1);
columns = 1:numel(countsCol);
columnsStd = (columns - mean(columns)) / std(columns);
countsRow = sum(countsTotCrop,2);
rows = 1:numel(countsRow);
rowsStd = (rows - mean(rows)) / std(rows);

figure

subplot(1,2,1)
hold on
plot(rowsStd,countsRow, '-b', 'userdata', 'Image Height')
plot(columnsStd,countsCol, '--r', 'userdata', 'Image Length')
hold off
legend(get(gca, 'children'), get(get(gca, 'children'), 'userdata')) 
title('Observed counts - Standardized')

subplot(1,2,2)
hold on
[fRows,xiRows] = ksdensity(countsRow,rowsStd);
plot(xiRows, fRows, '-b', 'userdata', 'Image Height')
[fCols,xiCols] = ksdensity(countsCol,columnsStd);
plot(xiCols, fCols, '--r', 'userdata', 'Image Length')
hold off
legend(get(gca, 'children'), get(get(gca, 'children'), 'userdata')) 
title('Kernel density function')
