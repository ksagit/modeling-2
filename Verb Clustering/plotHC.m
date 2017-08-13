function plotHC(tree,Q,cWin,cc)

%%% Plots the hierarchical tree obtained from runHC.m
%%% as well as modularity values for each partition licensed by the tree,
%%% and shows the partition that optimizes modularity

%%% INPUT:
%%% tree = a Nx3 matrix, containing information about the structure of 
%%%        the hierarchical tree (type "help linkage" for details)
%%% Q = Nx1 vector with modularity (or modularity-density) values for each parition
%%%     licensed by tree (first value is the modularity for singleton clustering;
%%%     last value is modularity for a single cluster solution)
%%% cWin = Nx1 vector with cluster assignments for the highest modularity
%%% cc = cophenetic correlation (how faithfully tree captures the original data)

%%% Idan Blank, Aug 8 2017; EvLab rulz!

N = length(Q);  % number of nodes

f = figure;
clf reset
set(f, 'units', 'normalized', 'position',[0.1 0.1 0.8*0.6 0.8]); 

%% Modularity plot %%
subplot(10,1,1:2)
hold on
xVals = zeros(N,1);
xVals(1) = 0.9*tree(1,3);
for i = 1:(N-2)
    xVals(i+1) = tree(i,3) + 0.5*(tree(i+1,3)-tree(i,3));
end
xVals(end) = 1.025*tree(end,3);
[xVals,i] = sort(xVals,'ascend');
Q = Q(i);
plot(xVals, Q, '-ko');
set(gca, 'xlim', [0.8*min(tree(:,3)), 1.05*max(tree(:,3))]);
set(gca, 'ylim', [min(Q)-sign(min(Q))*0.5*min(Q), 1.05*max(Q)]);
ylabel('Modularity');
set(gca, 'xtick', []);

%% Find best clustering and set up an appropriate color scheme %%
nClusters = length(unique(cWin));                           % number of clusters (ideally would be equal to ind)
colors = colormap('prism');                                 % "prism" coloring scheme
colors = colors(1:6,:);                                     % there are only 6 unique colors in this color scheme (repeating pattern)
nColors = size(colors,1);
nRep = ceil(nClusters/nColors);                             % how many times the color scheme needs to be replicated (with slight changes)
if nRep == 1
    inds = [1 4 5 3 6 2];
    colors = colors(inds(1:nClusters),:);
else
    for i = 1:(nRep-1)
        for j = 1:nColors
            j1 = j;
            j2 = mod(j,nColors)+1;
            color1 = colors(j1,:);
            color2 = colors(j2,:);
            colorDiff = color2-color1;
            colors(i*nColors+j,:) = color1 + (i/nRep)*colorDiff;
        end       
    end
end
for i = 2:nRep                                              % shuffle the order of each nColors colors
    inds = (i-1)*nColors+(1:nColors);
    currColors = colors(inds,:);
    colors(inds(1),:) = currColors(mod(i-1,nColors)+1,:);
    colors(inds(2:end),:) = currColors(setdiff(randperm(nColors),mod(i-1,nColors)+1),:);
end
newColors= colors;

%% Tree plot %%
subplot(10,1,3:10)
set(gca,'fontname','calibri','fontsize',5);
colorThres = xVals(Q == max(Q));                                                % for coloring different clusters in different colors
d = dendrogram(tree, 0, 'colorthreshold', colorThres, 'orientation', 'right');  % d is a vector of line handles
hold on

%% Search through brunches of the tree until you find the color of each leaf %%
oldColors = zeros(N,3);                               % current cluster colors to be changed later
branchInd = 1;
row = 1;
while (branchInd <= length(d)) && (row <= N)
    x = get(d(branchInd), 'xData');                   % plotting information for current branch    
    if x(1) == 0
        oldColors(row,:) = get(d(branchInd), 'Color');
        row = row+1;
    end
    if x(4) == 0
        oldColors(row,:) = get(d(branchInd), 'Color');
        row = row+1;
    end
    branchInd = branchInd + 1;
end
goodCols = sum(oldColors,2)>0;           % find where the black color is in newColors
oldColors = oldColors(goodCols,:);

row = 1;
while row <= size(oldColors,1)
    if sum(sum(repmat(oldColors(row,:),row-1,1) == oldColors(1:(row-1),:),2) == 3) > 0
        oldColors = oldColors([1:(row-1), (row+1):end],:);
    else
        row = row+1;
    end
end
    
%% Swap oldColors for newColors %%
for branchInd = 1:length(d)
    currColor = get(d(branchInd), 'Color');
    row = find(sum(oldColors == repmat(currColor,size(oldColors,1),1),2) == 3, 1);          % color of current branch
    if ~isempty(row)                                                                        % if branch is not black
        set(d(branchInd), 'Color', newColors(row,:));
    end
end

%% Add a line cutting through the tree at the approrpriate level + finalize %%
xVal = xVals(Q==max(Q));
yLims = get(gca,'ylim');
plot([xVal, xVal], [0 yLims(2)], '--k');

set(gca, 'xlim', [0.8*min(tree(:,3)), 1.05*max(tree(:,3))]);
s = suptitle(['Dendrogram (cophenetic correlation = ', num2str(round(1000*cc)/1000), '), ', ...
    num2str(nClusters), ' clusters (Q = ', num2str(round(1000*max(Q))/1000), ')']);
set(s, 'fontname', 'calibri', 'fontsize', 14, 'fontweight', 'bold');