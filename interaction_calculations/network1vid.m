%% Claire Guerin - 08/08/2017
% Import data

[trackFile,pathname] = uigetfile('*.mat','Select the tracking file');
cd(pathname)
S = load(trackFile);

[colony, date, time, chbr, trim] = strread(trackFile, '%s %s %s %s %s', 'delimiter','_');
chbr = cellstr(chbr{1}(1:end-4));

[num, txt, raw,fileS, scales] = getFiles(pathname);
[sizes, shiftFrame, shiftLine, delete] = getInfo(pathname,trackFile,colony,raw);

ovscores = csvread([colony{1}(end),'ovariescore.csv']);

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

stepwiseDistMat = stepwiseDistance(coordFixed, scaleFactor, fRate);

%% CALCULATE PAIRWISE DISTANCES

nodNames = string(taglist);
indivPairs = nchoosek(1:popSize,2);
% queenID = 2123;
% queenIndex = find(taglist == queenID);
% nodNames(queenIndex) = string('Q');
[~,scoredInd] = ismember(ovscores(:,1),taglist);
ovaryScore = nan(size(taglist));
ovaryScore(scoredInd) = ovscores(:,2);
uniqueScores = unique(ovscores(:,2))';
scoreColor = summer(numel(uniqueScores));
nodeColor = nan([popSize,3]);

for node = 1:popSize
    indivScore = ovaryScore(node);
    if isnan(indivScore)
        nodeColor(node,:) = [0,0,0];
    else
        nodeColor(node,:) = scoreColor(uniqueScores == indivScore,:);
    end
end

meanDist = pairwiseDistancesForNetwork(coordFixed, scaleFactor);

%% CREATE GRAPH

G = graph(indivPairs(:, 1),indivPairs(:, 2));
G.Nodes.Names = nodNames;
G.Nodes.Score = ovaryScore;
G.Nodes.Colors = nodeColor;

infAddOn = 0.0000000001;
G.Edges.Weight = abs(max(meanDist)+infAddOn - meanDist');
G2 = rmedge(G,find(isnan(meanDist)));
% G2.Edges.NormWeight = G2.Edges.Weight/sum(G2.Edges.Weight);
LWidths = 5*G2.Edges.Weight/max(G2.Edges.Weight);

fig1 = figure(1);clf;
set(fig1,'defaulttextinterpreter','latex','Color','w');
% p = plot(G2,'Layout','circle','EdgeLabel',ceil(G2.Edges.Weight),'LineWidth',LWidths,'EdgeColor', [1/3,1/3,1/3],'NodeLabel',cellstr(G2.Nodes.Names),'NodeColor',G.Nodes.Colors);
p = plot(G2,'Layout','circle','LineWidth',LWidths,'EdgeColor', [1/3,1/3,1/3],'NodeLabel',cellstr(G2.Nodes.Names),'NodeColor',G.Nodes.Colors);
% highlight(p,queenIndex)
axis equal
title('A. Degree centrality scores - weighted','FontSize',30)
set(gca,'xtick',[],'ytick',[])
p.MarkerSize = 10;
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
legendLabels = sprintf('%s \n %d \n %d \n %d \n %d','NaN',uniqueScores);
p.DisplayName = legendLabels;
% legendNames = cellfun(@num2str, num2cell(uniqueScores), 'UniformOutput', false);
% legend('Location','southwest')

% txt(-1,-1,'test')

% layout(p,'layered')
% 
% headers = {'Event' 'Node' 'Edge'};
% bfsearchT = cell2table(cell(0,3));
% bfsearchT.Properties.VariableNames = headers;
% for node = 1:popSize
%     T = bfsearch(G2,node,'allevents');
%     bfsearchT = [bfsearchT;T];
% end
% 
% for eachscore = 1:numel(uniqueScores)
%     scoretotest = uniqueScores(eachscore);
%     zeroScores = find(ovaryScore == scoretotest);
%     findZeros = ismember(table2array(bfsearchT(bfsearchT.Event == 'startnode',2)),zeroScores);
%     disp(bfsearchT(findZeros,:))
% end

deg_ranks = centrality(G2, 'degree', 'Importance', G2.Edges.Weight);
edges = linspace(min(deg_ranks),max(deg_ranks),7);
bins = discretize(deg_ranks,edges);
p.MarkerSize = 2*bins;

fig2 = figure(2);clf;
set(fig2,'defaulttextinterpreter','latex','Color','w');
p2 = plot(G2,'Layout','circle','LineWidth',LWidths,'EdgeColor', [1/3,1/3,1/3],'NodeLabel',cellstr(G2.Nodes.Names),'NodeColor',G.Nodes.Colors);
% highlight(p2,queenIndex)
axis equal
set(gca,'xtick',[],'ytick',[])
p2.MarkerSize = 10;
p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legendLabels = sprintf('%s \n %d \n %d \n %d \n %d','NaN',uniqueScores);
p2.DisplayName = legendLabels;
wcc = centrality(G2,'closeness','Cost',G2.Edges.Weight);
p2.NodeCData = wcc;
colormap jet
colorbar
caxis([0,2*10^(-6)])
title('B. Closeness centrality scores - weighted','FontSize',30)

% fig3 = figure(3);clf;
% set(fig3,'defaulttextinterpreter','latex','Color','w');
% scatter(ones(1,5),1:5,200,[0,0,0;scoreColor],'filled')
% ylim([0,7])
% axis equal
% addlegend = cellfun(@num2str, [NaN,num2cell(uniqueScores)], 'UniformOutput', false);
% text(ones(1,5)+1,1:5, addlegend,'FontSize',30,'Color','black', 'BackgroundColor', 'w')
% text(0.5,6,'Ovarian Development Score','FontSize',30,'Color','black', 'BackgroundColor', 'w')
% export_fig 'networklegend.bmp' -m2


%% COMPARE WITH PROBA WEIGHT
pathname = 'H:\Academia\BumbleBees2016\Behav_Ovaries\Behav\Odyssey\allFiles\track\';
cd(pathname)
S2 = load([pathname, trackFile, '_interactions.mat']);
ellps = S2.interEllipses;
meanEllps = nanmean(ellps,3);
prob = S2.interProbabilities;
meanProb = nanmean(prob,3);

% taglist = S2.taglist;
% nodNames = string(S2.taglist);
% popSize = numel(S2.taglist);
% indivPairs = nchoosek(1:popSize,2);
% queenIndex = find(taglist == queenID);
% % nodNames(queenIndex) = string('Q');
% [~,scoredInd] = ismember(ovscores(:,1),taglist);
% ovaryScore = nan(size(taglist));
% ovaryScore(scoredInd) = ovscores(:,2);
% uniqueScores = unique(ovscores(:,2))';
% scoreColor = summer(numel(uniqueScores));
% nodeColor = nan([popSize,3]);
% 
% for node = 1:popSize
%     indivScore = ovaryScore(node);
%     if isnan(indivScore)
%         nodeColor(node,:) = [0,0,0];
%     else
%         nodeColor(node,:) = scoreColor(uniqueScores == indivScore,:);
%     end
% end


G = graph(indivPairs(:, 1),indivPairs(:, 2));
G.Nodes.Names = nodNames;
G.Nodes.Score = ovaryScore;
G.Nodes.Colors = nodeColor;

meanProb(1:10,1:10)
[I,J] = find(meanEllps > 0);
ellpsWeight = [indivPairs,nan([size(indivPairs,1),1])];
assocPairs1 = [I,J];
assocPairs2 = [J,I];

for n = 1:size(assocPairs1,1)
    test1 = find(sum(ismember(ellpsWeight(:,1:2),assocPairs1(n,:)),2) == 2);
    test2 = find(sum(ismember(ellpsWeight(:,1:2),assocPairs2(n,:)),2) == 2);
    
    if numel(test1) == 1
        ellpsWeight(test1,3) = meanEllps(assocPairs1(n,1),assocPairs1(n,2));
    end
    
end

G3 = G;
G3.Edges.Weight = ellpsWeight(:,3);
G4 = rmedge(G3,find(isnan(ellpsWeight(:,3))));
% G2.Edges.NormWeight = G2.Edges.Weight/sum(G2.Edges.Weight);
LWidths4 = 5*G4.Edges.Weight/max(G4.Edges.Weight);

fig1 = figure(1);clf;
set(fig1,'defaulttextinterpreter','latex','Color','w');
p = plot(G4,'Layout','circle','LineWidth',LWidths4,'EdgeColor', [1/3,1/3,1/3],'NodeLabel',cellstr(G4.Nodes.Names),'NodeColor',G4.Nodes.Colors);
% highlight(p,queenIndex)
axis equal
title('C. Degree centrality scores - weighted','FontSize',30)
set(gca,'xtick',[],'ytick',[])
p.MarkerSize = 10;
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
legendLabels = sprintf('%s \n %d \n %d \n %d \n %d','NaN',uniqueScores);
p.DisplayName = legendLabels;

deg_ranks = centrality(G4, 'degree', 'Importance', G4.Edges.Weight);
edges = linspace(min(deg_ranks),max(deg_ranks),7);
bins = discretize(deg_ranks,edges);
p.MarkerSize = 2*bins;

fig2 = figure(2);clf;
set(fig2,'defaulttextinterpreter','latex','Color','w');
p2 = plot(G4,'Layout','circle','LineWidth',LWidths4,'EdgeColor', [1/3,1/3,1/3],'NodeLabel',cellstr(G4.Nodes.Names),'NodeColor',G4.Nodes.Colors);
% highlight(p2,queenIndex)
axis equal
set(gca,'xtick',[],'ytick',[])
p2.MarkerSize = 10;
p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legendLabels = sprintf('%s \n %d \n %d \n %d \n %d','NaN',uniqueScores);
p2.DisplayName = legendLabels;
wcc = centrality(G4,'closeness','Cost',G4.Edges.Weight);
p2.NodeCData = wcc;
colormap jet
colorbar
caxis([0,14])
title('D. Closeness centrality scores - weighted','FontSize',30)