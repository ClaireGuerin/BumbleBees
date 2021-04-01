function meanDist = pairwiseDistancesForNetwork(coordinates, scaleFactor)

nFrames = size(coordinates,1);
popSize = size(coordinates,2);
nComb = popSize*(popSize-1)/2;

distances = nan(nFrames, nComb);

for frm = 1:nFrames
    x = coordinates(frm,:,1) ;
    y = coordinates(frm,:,2) ;
    
    distances(frm,:) = pdist([x',y']) * scaleFactor;
end

meanDist = nanmean(distances,1);
end
