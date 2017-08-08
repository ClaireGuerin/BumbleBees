[FileName, PathName, FilterIndex] = uigetfile('*.avi','Select the video file','MultiSelect','on');
nFiles = size(FileName,2);

%%

parpool(4)

for video = 1:nFiles
    vidName = FileName{video};
    [date, time, position, chbr] = strread(vidName, '%s %s %s %s', 'delimiter','_');
    chbr = cellstr(chbr{1}(1:end-4));
    vid = VideoReader([PathName,vidName]);
    im = readFrame(vid);
    
    sizes = linspace(1,20,20);
    thresh = linspace(0.5,10,20);
    comb = combvec(sizes, thresh);
    n = size(comb, 2);
    test = nan(1,n);
    
    parfor i = 1:n
        tempFilterSize = comb(1,i);
        tempBradleyThreshold = comb(2,i);
        tempCodes = locateCodes(im, 'vis', 0, 'threshMode', 1, 'bradleyFilterSize', [tempFilterSize,tempFilterSize], 'bradleyThreshold', tempBradleyThreshold);
        
        test(1,i) = size(tempCodes, 1);
    end
    
    [Y,I] = max(test);
    bestComb = comb(:,test == Y);
    bestFsize = mode(bestComb(1,:));
    bestThresh = mode(bestComb(2,:));
    disp(strcat('Position= ', position{1}, '; Chamber= ', chbr{1}))
    disp([bestThresh,bestFsize])
end

% codes = locateCodes(im, 'threshMode', 1, 'bradleyFilterSize', [bestComb(1,1),bestComb(1,1)], 'bradleyThreshold', bestComb(2,1));

