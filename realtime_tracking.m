%% Set up video parameters

imaqreset;

%vid = videoinput('pointgrey', 1, 'F7_Mono8_1280x960_Mode0');
vid = videoinput('pointgrey', 1, 'F7_Mono8_2448x2048_Mode0'); %high-res camera (5.0 mp point grey blackfly)

%
%preview(vid);

src = getselectedsource(vid);

vid.FramesPerTrigger = 1;

src.ShutterMode = 'Manual';

src.Shutter = 80;

src.SharpnessMode = 'Manual';

src.Sharpness = 500;

src.GammaMode = 'Manual';

src.Gamma = 1.5;

src.GainMode = 'Manual';

src.Gain = 15;

src.ExposureMode = 'Manual';

src.Exposure = 1;

src.Brightness = 6.4453;

vid.FramesPerTrigger = inf;
triggerconfig(vid,'manual');

% %
% src = getselectedsource(vid);
% 
% vid.FramesPerTrigger = inf;
% 
% src.Exposure = 2.413635;
% src.Brightness = 5;
% src.Gain = 10;
% src.Shutter = 16.621;
% src.Gamma = 1.4;
% 
% triggerconfig(vid,'manual');
% 
% vidRes = vid.VideoResolution;
% nBands = vid.NumberOfBands;
% frameRate = src.FrameRate; %not a propety for high-res cam?
% 
 start(vid)

%% Background Calculation

BAD = 10;   % background acquisition duration (in seconds)
depth = 5;    % number of images to use for creating background
backIm = calculateBackground(vid,BAD,depth);
imshow(uint8(backIm));

%% Record experiment & track bees simultaneously

expDuration = 60;  % duration of experiment in seconds
intThresh = 5;      % intensity threshold used to turn grey-scale,   
                    % background subtracted image into BW image 
                    % (B = bee, W = bg)
                    
%Gather user input on approximate bee size (maybe switch this to a fixed
%value when we're running this a lot?
%%
im = peekdata(vid,1);
imshow(im);
title('Please indicate approximate bee size');
[x, y] = ginput(2);
close ALL
radius = sqrt(sum([(x(1)-x(2)).^2 (y(1)-y(2)).^2]));

%%
tic
while toc < expDuration

    hold on
    im = peekdata(vid,1);
   
    [beePos, beeRad] = spotBees(im, backIm, intThresh, radius);
    for i = 1:size(beePos,1)
        beeX = beePos(i,1);
        beeY = beePos(i,2);
        subIm = imcrop(im, [beeX-beeRad*2 beeY-beeRad*2 beeRad*4 beeRad*4]); %[xmin ymin width height]
        codes = locateCodes(subIm, 'vis', 0); % identify bee codes without showing cropped image
        
        if isempty(codes) == 0
            plot(beeX, beeY, 'r*');
            txt1 = num2str(codes.number);
            text(beeX,beeY,txt1, 'color', 'r')
        end
    end 
    
    hold off
    flushdata(vid)
end

stop(vid)

%% Test timing
backImU = uint8(backIm); %Convert background image to integer, faster caculation
cropSize = 100;
while 1
    %%
    tic
    im = peekdata(vid,1);
    imd = abs(im - backImU);
    imbw = imd > intThresh;
    
            codes = locateCodes(im, 'vis', 0); % identify bee codes without showing cropped image

   toc
   
  
end
