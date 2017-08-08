function [beePositions, beeRadius] = spotBees(currentImage, backImage, intensityThreshold, radius)

imDiff = abs(currentImage - uint8(backImage));

% Threshold into binary image
imBW = imDiff > intensityThreshold;


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