function [beePositions, beeRadius] = spotBees(currentImage, backImage, intensityThreshold)

imDiff = abs(currentImage - uint8(backImage));

% Threshold into binary image
imBW = imDiff > intensityThreshold;

imshow(currentImage);
title('Please indicate approximate bee size');
[x, y] = ginput(2);
close ALL
radius = sqrt(sum([(x(1)-x(2)).^2 (y(1)-y(2)).^2]));

se = strel('disk',double(uint8(radius/10)));
erodedI = imerode(imBW,se);
dilatedI = imdilate(erodedI,se);

s = regionprops(dilatedI,'centroid');   % Calculate centroids for connected 
                                        % components in the image using 
                                        % regionprops.
centroids = cat(1, s.Centroid);         % Concatenate structure array 
                                        % containing centroids into a 
                                        % single matrix.

beePositions = centroids;
beeRadius = radius;