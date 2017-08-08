[filename pathname]  = uigetfile('*avi', 'Select video file for tracking');
mov = VideoReader(strcat(pathname, '/',filename));
cd(pathname);
[tagfile tagpath] = uigetfile('*.csv', 'select taglist file');
taglist = csvread(strcat(tagpath, '/', tagfile), 1,0);
nframes = mov.NumberOfFrames;

% [camparfile camparpath] = uigetfile('*.mat', 'Load camera parameters');
% load(strcat(camparpath, camparfile));
%%
tags = taglist(:,1);
xcent = zeros(nframes, numel(tags));
ycent = nan(nframes, numel(tags));
frontx = nan(nframes, numel(tags));
fronty = nan(nframes, numel(tags));

ntags = numel(tags);
tic;
parpool(2);
%
parfor i = 1:nframes
    %%
    try
        im = read(mov,i);
        R = locateCodes(im, 'threshMode', 1, 'bradleyThreshold', 5, 'bradleyFilterSize', [12 12], 'vis', 0);
        %%
        rtags = [R.number];
        for j = 1:ntags
            
            rt = R(rtags == tags(j));
            
            if numel(rt) == 1
                xcent(i,j) = rt.Centroid(1);
                ycent(i,j) = rt.Centroid(2);
                frontx(i,j) = rt.frontX;
                fronty(i,j) = rt.frontY;
            end
        end
    catch
        continue
    end
    
    i
end
toc
%
trackingData = nan(nframes, ntags, 4);
trackingData(:,:,1) = xcent;
trackingData(:,:,2) = ycent;
trackingData(:,:,3) = frontx;
trackingData(:,:,4) = fronty;
%
save(strcat(filename, '_tracked.mat'), 'trackingData', 'taglist');

% %% Undistort data
% trackingDataUndistorted = trackingData;
% trackingDataUndistorted(trackingDataUndistorted == 0) = NaN;
% 
% trackingDataUndistortedcentU = trackingDataUndistorted(:,:,1:2);
% 
% for i = 1:size(trackingDataUndistortedcentU,1)
%     %%
%     i
%     for j = 1:size(trackingDataUndistortedcentU,2)
%         %%
%         %j
%         cent = [trackingDataUndistorted(i,j,1) trackingDataUndistorted(i,j,2)];
%         if ~isnan(cent(1))
%             centu = undistortPoints(cent, cameraParams);
%             trackingDataUndistortedcentU(i,j,1:2) = centu;
%         end
%         
%     end
% end
% trackingDataUndistorted = trackingDataUndistortedcentU*.2023/3032;
% 
% save(strcat(filename, '_undistorted.mat'), 'trackingData', 'trackingDataUndistorted', 'taglist');

%% Check tracking Data
for i = 101:200
    im = read(mov,i);
    xc = trackingData(i,:,1);
    yc = trackingData(i,:,2);
    imshow(im);
    hold on;
    plot(xc, yc, 'g.');
    text(xc+20, yc-10, num2str(tags), 'Color', 'g');
    hold off;
    pause(0.05);
end
