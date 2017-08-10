%% Claire Guerin - 08/08/2017
% Import data

[trackFile,pathname] = uigetfile('*.mat','Select the tracking file');
cd(pathname)
S = load(trackFile);

[colony, date, time, chbr, trim] = strread(trackFile, '%s %s %s %s %s', 'delimiter','_');
chbr = cellstr(chbr{1}(1:end-4));

[num, text, raw,fileS, scales] = getFiles(pathname);
[sizes, shiftFrame, shiftLine, delete] = getInfo(pathname,trackFile,colony,raw);

coordinates = S.trackingData;
taglist = S.taglist;
popSize = size(taglist,1);

nFrames = size(coordinates,1);
if strcmp(chbr{1},'FC')
    frameWidth = 2448;
    frameHeight = 2048;
else
    frameWidth = 1288;
    frameHeight = 964;
end

%% CHANGE ZEROs TO NaNs & CORRECT TRACKING DATA

coordCorr = trackingDataCorrection(coordinates, shiftFrame, shiftLine, frameWidth);

%% FILL-IN SMALL TRACKING GAPS

fRate = nFrames/(60*60); % number of frames per second (1H recording)
maxSecondGap = 2; % max time gap to interpolate (in seconds)
maxFrameGap = fRate * maxSecondGap; % max time gap to interpolate (in nFrames)
coordFixed = fixShortNanGaps(coordCorr, maxFrameGap);

%% ADAPT SIZES WITH SCALE

vidScale = scales(strcmp(fileS, trackFile(1:end-12)));
tagEdge = 0.11;
scaleFactor = tagEdge * vidScale;

stepwiseDistMat = stepwiseDistance(coordFixed, scaleFactor);

%% CALCULATE PAIRWISE DISTANCES

nodeNames = taglist;
indivPairs = nchoosek(1:popSize,2);

meanDist = pairwiseDistancesForNetwork(coordFixed, scaleFactor);

G = graph(indivPairs(:, 1),indivPairs(:, 2));
G.Nodes.Names = string(nodNames);
G.Nodes.Colors = pairColors;

G.Edges.Weight = meanDist';
G2 = rmedge(G,find(isnan(meanDist)));
G2.Edges.NormWeight = G2.Edges.Weight/sum(G2.Edges.Weight);
LWidths = 5*G2.Edges.Weight/max(G2.Edges.Weight);

subplot(3,4,network)
p = plot(G2,'Layout','circle','EdgeLabel',ceil(G2.Edges.Weight),'LineWidth',LWidths,'EdgeColor', [105/250,105/250,105/250],'NodeLabel',cellstr(G2.Nodes.Names),'NodeColor',G.Nodes.Colors);
title(graphTitle)
axis equal
set(gca,'xtick',[],'ytick',[])
p.MarkerSize = 10;