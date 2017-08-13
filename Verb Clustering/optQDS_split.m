function cNew = optQDS_split(A,c)

%%% First stage of QDS optimization: splitting communities %%%

%%% INPUT:
%%% A = a symmetric NxN matrix of weighted edges between N nodes
%%%     (or similarities between N items); must be between 0 and 1
%%% c = Nx1 vector, community / clusters assignments of the N nodes
%%%     (communities / clusters are denoted by integers)
%%% OUTPUT:
%%% cNew = new clustering solution after splitting communities

%%% Idan Blank, Aug 09, 2017; EvLab rulz! 

cNames = unique(c);                 % names of communities (ideally, 1:cN)
cN = length(cNames);                % number of communities
QDS = computeQDS(A,c);              % cNX1 vector, with the contribution of each community to QDS
cNew = zeros(size(c));              % initialization of final solution

for i = 1:cN
    currA = A(c==cNames(i),c==cNames(i));  % current community   
    N = size(currA,1);
    nodes = find(c==cNames(i));     % node indices
    deltaBest = 0;                  % QDS for best split of community i
    if N > 1
        D = diag(sum(currA));       % degree matrix: diagonal matrix with degree of each node
        L = D-currA;                    % Laplacian; can also use the command: laplacian(graph(A))
        v0 = ones(N,1);                
        opts.v0 = v0/norm(v0);          % normalized initial Lancsoz vector for eigenvalue decomposition
        opts.issym = 1;
        opts.isreal = 1;
        [V,~] = eigs(L,2,'sa',opts);
        vFiedler = V(:,2);              % eigenvector of the second smallest eigenvalue (smallest eigenvalue is 0 for a graph laplacian)
        zeroTest = round((10^6)*vFiedler)/(10^6);
        ind = find(zeroTest == 0, 1);
        if isempty(ind)                 % I don't know how to handle Fiedler vectors with entries equal to 0   
            [~,inds] = sort(vFiedler,'descend');
            nodes = nodes(inds);
            nodeBest = 0;                   % node where best split is identified
            for j = 1:(N-1)                 % all possible ways of cutting current community in 2
                c1 = nodes(1:j);      % nodes in first sub-community (the rest go in the second subcommunity)
                cSplit = c;
                cSplit(c1) = cN+1;          % new label for sub-community c1 (sub-community c2 keeps the original label of community i)            
                QDS_new = computeQDS(A,cSplit);                
                c1_ind = find(unique(cSplit)==cN+1,1);          % index of sub-community c1 in the QDS vector
                c2_ind = find(unique(cSplit)==cNames(i),1);     % index of sub-community c2 in the QDS vector
                QDS_split = QDS_new(c1_ind) + QDS_new(c2_ind);
                deltaQDS = QDS_split - QDS(c2_ind);                  % difference in QDS after splitting community i 
                if deltaQDS > deltaBest
                    deltaBest = deltaQDS;
                    nodeBest = j;
                end
            end
        end
    end
    if deltaBest > 0
        cNew(nodes(1:nodeBest)) = i;                    % sub-community 1 is named i in the new community structure
        cNew(nodes((nodeBest+1):end)) = cN+i;           % sub-community 2 is named cN+i in the new community structure
    else
        cNew(nodes) = cNames(i);                        % no split is introduced for community i
    end
end

cFinal = zeros(size(cNew));                             % renaming the communities to 1:cNewN
cNewNames = unique(cNew);
[~,inds] = sort(cNewNames,'ascend');
for i = 1:length(cNewNames)
    cFinal(cNew==cNewNames(inds(i))) = i;
end
cNew = cFinal;