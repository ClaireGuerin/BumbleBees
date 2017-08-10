function [coordCorr] = trackingDataCorrection(coordinates, shiftFrame, shiftLine, frameWidth)

coordinates(coordinates==0) = NaN;
popSize = size(coordinates,2);

for indiv = 1:popSize

	untracked = numel(find(isnan(coordinates(:,indiv,1))));

	if untracked > 0.99 * nFrames
		coordinates(:,indiv,:) = NaN; % bee wasn't tracked on this frame
	else
		position = [coordinates(:,indiv,1), coordinates(:,indiv,2)];
		firstPos = repmat(position(1,:), size(position,1), 1);
		movement = abs(position - firstPos);
		meanMovement = nanmean(movement);
	
		if meanMovement < 6
			coordinates(:,indiv,:) = firstPos; % bee didn't actually move
		end
	end
end

coordCorr = coordinates;
centerXs = coordCorr(:,:,1);
frontXs = coordCorr(:,:,3);

centerXs((shiftFrame+1):end,:) = mod(centerXs((shiftFrame+1):end,:)+ frameWidth - shiftLine,frameWidth);
frontXs((shiftFrame+1):end,:) = mod(frontXs((shiftFrame+1):end,:)+ frameWidth - shiftLine,frameWidth);

coordCorr(:,:,1) = centerXs;
coordCorr(:,:,3) = frontXs;

end