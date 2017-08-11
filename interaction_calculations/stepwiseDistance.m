function stepwiseDist = stepwiseDistance(trackingDataSmooth,scaleFactor,fRate)

nFrames = size(trackingDataSmooth,1);
nPop = size(trackingDataSmooth,2);
nCoord = size(trackingDataSmooth,3);
stepwiseDist = nan([nFrames,nPop,nCoord/2]);

for timeStep = 1:(nFrames - 1)
    x0 = trackingDataSmooth(timeStep,:,[1,3]);
    x1 = trackingDataSmooth(timeStep+1,:,[1,3]);
    y0 = trackingDataSmooth(timeStep,:,[2,4]);
    y1 = trackingDataSmooth(timeStep+1,:,[2,4]);
    
    stepwiseDist(timeStep,:,:) = sqrt((x1 - x0).^2 + (y1 - y0).^2) * scaleFactor * fRate;
end
end