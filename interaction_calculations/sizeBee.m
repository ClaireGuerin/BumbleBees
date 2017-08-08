% IMPORT DATA
colony = 'ColonyA';
path = strcat('F:\Academia\BumbleBees2016\Vid_Behav_Ovaries\Odyssey\');
cd(path)
ovaFile = 'ovariescore.csv';
scores = csvread(strcat(path,ovaFile),1,0);
popSize = size(scores,1);

%%
% GET BEE SIZES - USER INPUT

chamber = {'FC', 'NC'};
len = zeros(1,popSize);
width = zeros(1,popSize);
tagEdge = 0.11; % length in inches of a tag edge (for size scaling)

for i = 1:2
    chbr = char(chamber(i));
    [matFile, pathname]  = uigetfile('*mat', strcat('Select tracking data -',chbr));
    vidFile = strrep(matFile,'_tracked.mat','');
    S = load(strcat(path,matFile));
    vid = VideoReader(strcat(path,vidFile));
    coordinates = S.trackingData;
    allTags = S.taglist;
    toKeep = ismember(allTags, scores(:,1));
    coord = coordinates(:,toKeep,:);
    
    scaleIm = read(vid,1);
    imshow(scaleIm)
    title('Outline a tag edge', 'color', 'r');
    [xScale, yScale] = ginput(2);
    close ALL
    scale = sqrt((xScale(1)-xScale(2))^2 + (yScale(1)-yScale(2))^2);
    
    for bee = 1:popSize
        if len(bee) == 0 || width(bee) == 0
            detect = find(coord(:,bee,:) ~= 0);
            [frame, pos] = ind2sub([size(coord,1),size(coord,3)],detect(1));
            im = read(vid,frame);
            xCentroid = coord(frame,bee,1);
            yCentroid = coord(frame,bee,2);
            r = 20;
            th = 0:pi/50:2*pi;
            xunit = r * cos(th) + xCentroid;
            yunit = r * sin(th) + yCentroid;
            imshow(im);
            hold on
            h = plot(xunit, yunit, 'y', 'LineWidth', 2);
            hold off
            title(strcat('Draw body length and width of bee',num2str(bee)));
            [x, y] = ginput(4);
            close ALL
            len(1,bee) = tagEdge * sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2) / scale;
            width(1,bee) = tagEdge * sqrt((x(3)-x(4))^2 + (y(3)-y(4))^2) / scale;
        end
    end  
end

sizes = [transpose(scores(:,1));len;width];
output = strcat(path,'sizes.csv');
csvwrite(output, sizes);