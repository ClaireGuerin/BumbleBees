function [sizes, shiftFrame, shiftLine, delete] = getInfo(pathname,trackFile,colony)

sizesFile = strcat(colony,'sizes.csv');
sizes = csvread(strcat(pathname,char(sizesFile)));
allVideos = raw(2:size(raw,1),1);
%nChar = size(vidFile,2);
myVid = find(all(char(allVideos) == trackFile(1:end-16),2));
shiftFrame = cell2mat(raw(myVid+1,2));
shiftLine = cell2mat(raw(myVid+1,3));
delete = cell2mat(raw(myVid+1,4));

end