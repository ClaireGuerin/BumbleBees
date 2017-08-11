function [interactionsMat, orientationsMat] = interactionProbas(tag,bsize,distri,orientation,mnSiz,im)

map = hot;
colorSteps = size(map,1);
rangeValues = linspace(0, 1, colorSteps);

probWidth = mnSiz(2,2);
probLength = mnSiz(2,1);

xcenter = tag(1,:);
ycenter = tag(2,:);
xfront = tag(3,:);
yfront = tag(4,:);

wid = bsize(1,:);
len = bsize(2,:);
widthFactor = wid/ probWidth;
lengthFactor = len / probLength;
sizeFactor = nanmean([widthFactor; lengthFactor],1);

fWidth = im(1);
fHeight = im(2);
backGround = zeros(fHeight,fWidth);

% ORIENTATION

angle = atan2(yfront - ycenter, xfront - xcenter);
angle = angle .* (angle >= 0) + (angle + 2 * pi) .* (angle < 0); % to have the angle between 0 and 2pi.
rotationAngleDegrees = rad2deg(3*pi/2 - angle);

% INTERACTION DETECTION

bHere = find(~isnan(angle));
nPresentB = numel(bHere);

nPop = numel(xcenter);
interactionsMat = zeros(nPop,nPop);
orientationsMat = nan(nPop,nPop);

if nPresentB > 1
    
    allBeesDistri = repmat(backGround, 1, 1, nPresentB);
    allBeesOrient = nan(fHeight,fWidth, nPresentB);
    
    for bee = 1:nPresentB
        
        beeDistri = imresize(distri,sizeFactor(bHere(bee)), 'method', 'nearest');
        beeOrient = imresize(orientation,sizeFactor(bHere(bee)), 'method', 'nearest');
        imRot = imrotate(beeDistri, rotationAngleDegrees(bHere(bee)));
        imOrient = imrotate(beeOrient, rotationAngleDegrees(bHere(bee)));
        imCenter = ceil(size(imRot)/2);
        
        xOffset = ceil(xcenter(bHere(bee)) - imCenter(1));
        yOffset = ceil(ycenter(bHere(bee)) - imCenter(2));
        imXPos = 1 + xOffset: xOffset + size(imRot,2);
        imYPos = 1 + yOffset: yOffset + size(imRot,1);
        inRangeX = imXPos > 0 & imXPos <= size(im,2);
        inRangeY = imYPos > 0 & imYPos <= size(im,1);
        
        allBeesDistri(imYPos(inRangeY), imXPos(inRangeX), bee) = imRot(inRangeY,inRangeX);
        allBeesOrient(imYPos(inRangeY), imXPos(inRangeX), bee) = imOrient(inRangeY,inRangeX);
        
    end
    
    combos = nchoosek(1:nPresentB,2);
    
    for pair = 1:size(combos,1)
        bee1 = combos(pair,1);
        bee2 = combos(pair,2);
        interbee = allBeesDistri(:,:,bee1).*allBeesDistri(:,:,bee2);
        orientbee = allBeesOrient(:,:,bee1) + allBeesOrient(:,:,bee2);
        
        intB1 = bHere(bee1);
        intB2 = bHere(bee2);
        P = reshape(interbee, numel(interbee),1);
        Prob = max(P);
        interactionsMat(intB1,intB2) = Prob;

        [M, F] = mode(orientbee);
        [~, MaxInd] = max(F);
        orientationsMat(intB1,intB2) = M(MaxInd);
        
    end

end