%% Set up video parameters

imaqreset;
vid = videoinput('pointgrey', 1, 'F7_Mono8_1280x960_Mode0');
src = getselectedsource(vid);

vid.FramesPerTrigger = inf;

src.Exposure = 2.413635;
src.Brightness = 20.507812;
src.Gain = 10;
src.Shutter = 16.621;
src.Gamma = 1.4;

triggerconfig(vid,'manual');

vidRes = vid.VideoResolution;
nBands = vid.NumberOfBands;
frameRate = src.FrameRate;

start(vid)

%% Background Calculation

BAD = 20;   % background acquisition duration (in seconds)
depth = 30;    % number of images to use for creating background
backIm = calculateBackground(vid,BAD,depth);
imshow(backIm);

%% Record experiment & track bees simultaneously

expDuration = 120;  % duration of experiment in seconds
intThresh = 5;      % intensity threshold used to turn grey-scale,   
                    % background subtracted image into BW image 
                    % (B = bee, W = bg)

tic
while toc < expDuration
    im = peekdata(vid,1);
    [beePos, beeRad] = spotBees(im, backIm, intThresh);
    imshow(im)
    
    hold on
   
    for i = 1:size(beePos,1)
        beeX = beePos(i,1);
        beeY = beePos(i,2);
        subIm = imcrop(im, [beeX-beeRad*2 beeY-beeRad*2 beeRad*4 beeRad*4]); %[xmin ymin width height]
        codes = locateCodes(subImage, vis, 0); % identify bee codes without showing cropped image
        plot(beeX, beeY, 'r*');
        txt1 = codes.number;
        text(beeX,beeY,txt1, Color, [250 0 0])
    end 
    
    hold off
end
