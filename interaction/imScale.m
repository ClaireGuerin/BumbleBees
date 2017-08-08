%%
% IMPORT DATA

pathname = strcat('E:\Academia\BumbleBees2016\Vid_Behav_Ovaries\Odyssey\');
cd(pathname)

k = dir(strcat(pathname,'*.avi'));
filenames = {k.name}';
nFiles = size(filenames,1);
scales = nan(nFiles,1);

chambers = {};

%%
% SCALING IMAGES

for vidFile = 1:2
% for vidFile = 1:nFiles
    vid = VideoReader(strcat(pathname,filenames{vidFile}));
    chambers{vidFile} = filenames{vidFile}(end-5:end-4);
    scaleIm = read(vid,1);
    imshow(scaleIm)
    title('Outline a tag edge');
    [xScale, yScale] = ginput(2);
    close ALL
    scale = sqrt((xScale(1)-xScale(2))^2 + (yScale(1)-yScale(2))^2);
    scales(vidFile,1) = scale;
end
%%
% save('vidScales.mat', 'filenames', 'scales');
% vidScalesChambers = table(chambers', scales(1:2,1));
% writetable(vidScalesChambers)