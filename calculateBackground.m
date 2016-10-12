function backImage = calculateBackground(video, duration, bgDepth)
% create a background image calculated over a specific duration of video 
% pre-acquisition

herz = duration/bgDepth;
resolution = video.VideoResolution;
backCalc = nan(resolution(2),resolution(1),0);
i = 1;

tic
while toc < duration

    im = peekdata(video,1);
    backCalc(:,:,i) = im;
    pause(herz)
    flushdata(video);
    i = i + 1;
end

backImage = median(backCalc,3);