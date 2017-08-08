% function [allBeesDistri] = calculateInteractionsProbs(tag,bsize,distri,orientation,mnSiz,im,vis)
function [interactionsMat, orientationsMat] = calculateInteractionsProbs(tag,bsize,distri,orientation,mnSiz,im,vis)
% tag = tagPos;
% bsize = beeShape;
% distri = probBodyResized;
% orientation = frontMapResized;
% mnSiz = meanSize;
% im = fImage;
% vis = visualize;

map = hot;
% mapPlot = map((map(:,2) ~= 0 | map(:,3) ~= 0),:);
colorSteps = size(map,1);
rangeValues = linspace(0, 1, colorSteps);

% meanMeas = mnSiz(1,:);
probWidth = mnSiz(2,2);
probLength = mnSiz(2,1);
% centProb = mnSiz(3,:); % centroid pos
% sizePx = mnSiz(3,:);

xcenter = tag(1,:);
ycenter = tag(2,:);
xfront = tag(3,:);
yfront = tag(4,:);

wid = bsize(1,:);
len = bsize(2,:);
widthFactor = wid/ probWidth;
lengthFactor = len / probLength;
sizeFactor = nanmean([widthFactor; lengthFactor],1);

fHeight = size(im,1);
fWidth = size(im,2);
backGround = zeros(fHeight,fWidth);

% ORIENTATION

angle = atan2(yfront - ycenter, xfront - xcenter);
angle = angle .* (angle >= 0) + (angle + 2 * pi) .* (angle < 0); % to have the angle between 0 and 2pi.
% xbarycenter = xcenter + 1/4*sizePx(1)*cos(angle+pi);
% ybarycenter = ycenter + 1/4*sizePx(1)*sin(angle+pi);
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
    
    if vis

        subplot(1,2,1)
        imshow(im);
        hold on
        plot(xcenter, ycenter, '.', 'Color', map(1,:), 'MarkerSize', 15)
        subplot(1,2,2)
        imagesc(nansum(allBeesDistri,3));
        caxis([0, 2])
        axis off
%         axis equal
        axis image
        colormap(map)
        
    end
    
%     if vis == 1
%         
%         mask = nanmax(allBeesDistri,[],3);
%         colImage = cat(3, zeros(size(im)), zeros(size(im)), ones(size(im))); %Single color image, this one is blue
%         maxOpacity = 0.7; %Set the maximum opacity values across all pixels (must be 0-1)
%         imshow(im);
%         hold on
%         h = imshow(colImage);
%         set(h, 'AlphaData', mask.*maxOpacity);
% 
% %         heatMapData = 
% %         imagesc(heatMapData);
% %         hold on
% % %         positiveProb = heatMapData ~= 0;
% % %         alphaMap = double(positiveProb);
% % %         alphaMap(positiveProb) = 0.2;
% % %         alpha(h, alphaMap)
% %         caxis([0, 2])
% %         axis off
% % %         axis equal
% %         axis image
% %         colormap(map)
% %         h = imshow(im);
% %         set(h, 'AlphaData', 0.2)
% %         plot(xcenter, ycenter, '.', 'Color', map(1,:), 'MarkerSize', 15)
%     
%     end
    
    combos = nchoosek(1:nPresentB,2);
    
    for pair = 1:size(combos,1)
        bee1 = combos(pair,1);
        bee2 = combos(pair,2);
        interbee = allBeesDistri(:,:,bee1).*allBeesDistri(:,:,bee2);
        orientbee = allBeesOrient(:,:,bee1) + allBeesOrient(:,:,bee2);
        % 0 = posterior-posterior;
        % 1 = posterior-anterior;
        % 2 = anterior-anterior;
        % NaN = no px overlap
        
%         subplot(1,2,1)
%         imagesc(allBeesDistri(:,:,bee1))
%         axis equal
%         subplot(1,2,2)
%         imagesc(allBeesDistri(:,:,bee2))
%         axis equal
        
        intB1 = bHere(bee1);
        intB2 = bHere(bee2);
        P = reshape(interbee, numel(interbee),1);
%         P = 1-prod(1-P);
        Prob = max(P);
        interactionsMat(intB1,intB2) = Prob;

        [M, F] = mode(orientbee);
        [~, MaxInd] = max(F);
        orientationsMat(intB1,intB2) = M(MaxInd);
        
        if vis && P > 0
            
            valDiff = rangeValues - P;
            [~, indx] = min(abs(valDiff));
            subplot(1,2,1)
            plot(xcenter(:,[intB1,intB2]), ycenter(:,[intB1, intB2]), '-', 'Color', map(indx,:), 'LineWidth', 1)
            txt = num2str(P);
            text(mean(xcenter(:,[intB1,intB2])), mean(ycenter(:,[intB1, intB2])), txt, 'Color', map(indx,:))
            
        end
    end

end

%%

%[x, y] = ind2sub([frameWidth,frameHeight],px);
%presence = x^2./a.^2 + y^2./b.^2 <= 1 + interaction_range;
%temp = full(pxVal,px);
%pxVal(px,:) = presence;
%pxVal = sparse(pxVal);
%clearvars im xcenter ycenter xfront yfront a b slope angle x1 x2 y1 y2 t x y w x y bHere nPresentB BW combos interbee