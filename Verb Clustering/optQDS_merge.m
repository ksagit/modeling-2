function cNew = optQDS_merge(A,c)

%%% Second stage of QDS optimization: merging communities %%%

%%% INPUT:
%%% A = a symmetric NxN matrix of weighted edges between N nodes
%%%     (or similarities between N items); must be between 0 and 1
%%% c = Nx1 vector, community / clusters assignments of the N nodes
%%%     (communities / clusters are denoted by integers)
%%% OUTPUT:
%%% cNew = new clustering solution after merging communities

%%% Idan Blank, Aug 09, 2017; EvLab rulz! 

cNames = unique(c);                 % names of communities (ideally, 1:cN)
cN = length(cNames);                % number of communities
QDS = computeQDS(A,c);              % cNX1 vector, with the contribution of each community to QDS
mergeMat = zeros(cN,cN);            % upper triangular; (i,j) = deltaQDS for merging communities i,j
for i = 1:(cN-1)
    c1 = cNames(i);
    for j = (i+1):cN
        c2 = cNames(j);
        Ac1c2 = A(c==c1,c==c2);
        if sum(Ac1c2(:)) > 0
            cMerge = c;
            cMerge(c==c1) = cN+1;               % merging communities i and j
            cMerge(c==c2) = cN+1;
            QDS_merge = computeQDS(A,cMerge);   % index of merged community c1 in the QDS_merge vector
            mergeInd = find(unique(cMerge)==cN+1,1);
            deltaQDS = QDS_merge(mergeInd) - QDS(i) - QDS(j);            
            if deltaQDS > 0
                mergeMat(i,j) = deltaQDS;
            end
        end
    end
end

mergePairs = zeros(cN,2);
row = 1;
isStop = sum(mergeMat(:)) == 0;
while ~isStop
    [i,j] = ind2sub([cN, cN], find(mergeMat == max(mergeMat(:)),1));    % community pair whose merge will result in the highest QDS increase
    if sum(mergePairs(:)==i)==0 && sum(mergePairs(:)==j)==0
        mergePairs(row,:) = [cNames(i),cNames(j)];                      % i and j will be merged
        row = row+1;
        mergeMat(i,:) = 0;                                              % i cannot be merged with anything else
        mergeMat(j,:) = 0;                                              % j cannot be merged with anything else
        mergeMat(:,i) = 0;
        mergeMat(:,j) = 0;
    end
    if sum(mergeMat(:)) == 0
        isStop = 1;
    end
end
ind = find(sum(mergePairs,2)==0,1);
mergePairs = mergePairs(1:(ind-1),:);

cNew = c;
for i = 1:size(mergePairs,1)
    c1 = mergePairs(i,1);
    c2 = mergePairs(i,2);
    cNew(c==c2) = c1;
end

cFinal = zeros(size(cNew));                             % renaming the communities to 1:cNewN
cNewNames = unique(cNew);
[~,inds] = sort(cNewNames,'ascend');
for i = 1:length(cNewNames)
    cFinal(cNew==cNewNames(inds(i))) = i;
end
cNew = cFinal;