[filename, pathname] = uigetfile('*.xlsx', 'Select reverse frame excel spreadsheet info');
cd(pathname)

[num, text, raw] = xlsread(filename);
videos = cell2mat(raw(2:size(raw,1),1));
frames = cell2mat(raw(2:size(raw,1),2));
lines = cell2mat(raw(2:size(raw,1),3));
delete = cell2mat(raw(2:size(raw,1),4));

for i = 1:size(raw,1)-1
    vidFile = videos(i,:);
    colony = vidFile(7);
    vid = VideoReader(strcat(colony,'\',vidFile,'.avi'));
    nFrames = vid.NumberOfFrames;
    im = read(vid,frames(i)+10);
    imshow(im)
    [x,y] = ginput(1);
    limit = x;
    close ALL
    xlswrite(filename,limit,'Feuil1',strcat('C',num2str(i+1)))
end