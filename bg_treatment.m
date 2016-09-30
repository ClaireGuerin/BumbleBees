[filename, pathname] = uigetfile('*.avi'); % user input: movie file name and path

mov = VideoReader([pathname filename]); % import file

%%
nframes = mov.NumberOfFrames; % number of frames in the video

for i = 1:nframes % for each frame
%%
    im = read(mov,i); % read frame
    imshow(im); % show frame

    % insert tracking code
    figure(1);
    codes = locateCodes(im);
    figure(2);
    tic
    codes = locateCodes(im, 'threshMode', 1, 'bradleyFilterSize', [10 10], 'bradleyThreshold', 3);
    toc % processing delay
    drawnow
end

im = rgb2gray(im); % convert image from color to grey scale

%frame1 = read(mov,1);
%col = [];
%allStd = []

%%
%for i = 1:size(frame1,1);
%    for j = 1:size(frame1,2);
%        for k = 1:nframes;
%            col(k) = im(i,j);
%        end
%        allStd(i,j) = std(col);
%    end
%end
%%

nBackFrames = 20;
backFrameInds = round(linspace(1,nframes,nBackFrames)); % regularly spaced-out frame indices

%Create empty matrix for background calculation
backCalc = nan(mov.Height, mov.Width,nBackFrames);
for i = 1:nBackFrames
    im = read(mov,backFrameInds(i));
    imshow(im);
    im = rgb2gray(im);
    backCalc(:,:,i) = im;
end

backImage = median(backCalc,3); % try also mode or mean --> which one renders a better background mask?

%% Example background subtraction for single image

frame = 200;
imCur  = read(mov,frame);
disp(['Frame #',num2str(frame)])
imCur = rgb2gray(imCur);

imDiff = abs(imCur - uint8(backImage));

intThresh = 5;

%Threshold into binary image
imBW = imDiff > intThresh;

%%Test how bg is altered by erosion and dilation
test = uint8(backImage);
se = offsetstrel('ball',5,5);
erodedI = imerode(test,se);
dilatedI = imdilate(test,se);
erodilatedI = imdilate(dilatedI,se);
figure(1);
subplot(2,2,1);
imshow(test); title('original background')
subplot(2,2,2);
imshow(erodedI); title('eroded background')
subplot(2,2,3);
imshow(dilatedI); title('dilated background')
subplot(2,2,4);
imshow(erodilatedI); title('eroded & dilated background');
%%

%%Compare efficiency in tag recognition under different image treatments
imDiff2 = abs(imCur - uint8(erodilatedI));
figure(2);
subplot(3,2,1);
imshow(imCur); title('Original Image');
subplot(3,2,2);
tic
codes = locateCodes(imCur); title('Identified tags')
toc
subplot(3,2,3)
imshow(imDiff); title('Filtered Img')
subplot(3,2,4);
tic
codes = locateCodes(imDiff); title('Identified tags')
toc
subplot(3,2,5)
imshow(imDiff2); title('Filtered Img, Eroded+Dilated Bg')
subplot(3,2,6);
tic
codes = locateCodes(imDiff2); title('Identified tags')
toc
%%

%%Extract chunks of images containing single bees for individual tag
%%identification
%%1: Create Background

nBackFrames = 20;
backFrameInds = round(linspace(1,nframes,nBackFrames)); % regularly spaced-out frame indices

%Create empty matrix for background calculation
backCalc = nan(mov.Height, mov.Width,nBackFrames);
for i = 1:nBackFrames
    im = read(mov,backFrameInds(i));
    imshow(im);
    im = rgb2gray(im);
    backCalc(:,:,i) = im;
end

backImage = median(backCalc,3); % try also mode or mean --> which one renders a better background mask?

%%2: Spot bees

frame = 200;
imCur  = read(mov,frame);
disp(['Frame #',num2str(frame)])
imCur = rgb2gray(imCur);

imDiff = abs(imCur - uint8(backImage));
%%
intThresh = 10;

%Threshold into binary image
imBW = imDiff > intThresh;
figure(1)
imshow(imCur)
figure(2)
imshow(imDiff)
figure(3)
imshow(imBW)

se = strel('line',11,90);
erodedI = imerode(imBW,se);
dilatedI = imdilate(erodedI,se);

s = regionprops(dilatedI,'centroid'); %Calculate centroids for connected components in the image using regionprops.
centroids = cat(1, s.Centroid); %Concatenate structure array containing centroids into a single matrix.

%Display binary image with centroid locations superimposed.
imshow(dilatedI) 
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
hold off