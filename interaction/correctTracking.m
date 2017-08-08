%%
% Import the files needed
[matFile, path]  = uigetfile('*mat', 'Select tracking data');
cd(path)
S = load(strcat(path,matFile));
vidFile = strrep(matFile,'_tracked.mat','');
vid = VideoReader(strcat(path,vidFile));
colony = vidFile(1:7);
sizes = csvread([colony,'sizes.csv']);
[num, text, raw] = xlsread('reverseFrames.xlsx');

coordinates = S.trackingData;
allTags = S.taglist;
toKeep = ismember(allTags, sizes(1,:));
coordPlain = coordinates(:,toKeep,:); %only keep tracking data of bees for which we have ovary scores (and sizes)

%%

videos = raw(2:size(raw,1),1);
nChar = size(vidFile,2);
myVid = find(all(char(videos) == vidFile(1:nChar-4),2));

%%

zz = myVid; 

coordinates(coordinates == 0) = NaN; %Remove zeros
movWidth = vid.Width;

frame = cell2mat(raw(zz+1,2));
line = cell2mat(raw(zz+1,3));
delete = cell2mat(raw(zz+1,4));

coordCorr = coordinates;
xs = coordCorr(:,:,1); %isolate x-coordinates


subplot(2,1,1);
plot(xs);
hold on
plot([frame frame], [0 movWidth], 'r-.','LineWidth', 2);
title('uncorrected');
hold off



subplot(2,1,2);
xs((frame+1):end,:) = mod(xs((frame+1):end,:)+ movWidth - line,movWidth);
plot(xs);

hold on
plot([frame frame], [0 movWidth], 'r-.','LineWidth', 2);
title('corrected');
hold off

%load back into coordCoor:
coordCorr(:,:,1) = xs;


%% Verify correction on video

nFrames = vid.NumberOfFrames;
nIndividuals = numel(allTags);
colorsOrder = jet(nIndividuals);
colorsShuffle = colorsOrder(randperm(nIndividuals),:);

fig = figure(1);clf;
set(fig,'defaulttextinterpreter','latex','Color','w');
%
%for i = (frame - 30):(frame+40)
for i = 10149
    %%
    subplot(1,2,1);
    im = read(vid,i);
    im = rgb2gray(im);
    imshow(im);
    
    hold on;
    scatter(nanmean(coordinates((i - 10):i,:,1)), nanmean(coordinates((i - 10):i,:,2)), [], colorsOrder, 'filled');
    hold off
    %title(strcat('uncorrected, frame:', num2str(i)))
    title('corrupted frame')
    
    subplot(1,2,2);
    
    %Correct image
    if i >= frame
        c1 = im(:,1:floor(line));
        c2 = im(:,ceil(line):end);
        imc = [c2 c1];
    else
        imc = im;
    end
    imshow(imc);
    
    hold on;
    scatter(nanmean(coordCorr((i - 10):i,:,1)), nanmean(coordCorr((i - 10):i,:,2)), [], colorsOrder, 'filled');
    hold off
    %title(strcat('corrected, frame:', num2str(i)))
    title('corrected coordinates')
    drawnow
    pause(0.1);
end
% %%
% [I1, I2, I3] = find(coordPlain(frame+1:nFrames,:,1)>line);
% [J1, J2, J3] = find(coordPlain(frame+1:nFrames,:,1)<line);
% %%
% for pos = size(I1,1)
%     coordCorr(frame+I1(pos),I2(pos),1) = coordPlain(frame+I1(pos),I2(pos),1) - line;
% end
%
% for pos = size(J1,1)
%     coordCorr(frame+J1(pos),J2(pos),1) = coordPlain(frame+J1(pos),J2(pos),1) + (frameWidth - line);
% end
%

export_fig 'H:\Academia\MEME\S3 - HARVARD\Report\Figs\revFrames.bmp' -m2


