function backImage = backgroundCalculation(video, duration, herz)
% calculate background over frames within
% background acquisition time.

resolution = video.VideoResolution;
backCalc = nan(resolution(2),resolution(1),0);
i = 1;

tic
while toc < duration

    im = peekdata(video,1);
    backCalc(:,:,i) = im;
    pause(1000/herz)
    flushdata(video);
    i = i + 1;
end

backImage = median(backCalc,3);