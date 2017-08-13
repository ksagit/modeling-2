function [tree, clusters, cc] = runGreedyOptQ(A,D)

%%% Hierarchical clustering based on modularity optimization
%%% INPUT:
%%% A = a symmetric NxN matrix of weighted edges between N nodes
%%%     (or similarities between N items)
%%% D = a symmetric NxN matrix of distances between the same N items
%%% OUTPUT:
%%% tree = a Nx3 matrix, containing information about the structure of 
%%%        the hierarchical tree (type "help linkage" for details)
%%% clusters = NxN matrix, each column is a different partition of the N nodes
%%%            obtained by cutting the tree at some level; all possible
%%%            partitions are computed (from each node in its own cluster, through 
%%%            to all nodes in a single cluster)
%%% cc = cophenetic correlation (how faithfully tree captures the original similarities)

%%% Idan Blank, August 8, 2017; EvLab rulz!

M = sum(sum(A))/2;                          % sum of all pairwise similarities 
N = size(A,1);
tree = zeros(N-1,3);                        % similar to what the linkage function returns
rowT = 1;                                   % for filling in tree

clusters = zeros(N);
clusters(:,1) = (1:N)';                     % initial clustering: each item is assigned to a different cluster
for i = 1:(N-1)
    maxC = max(clusters(:,i));              % for labeling new (i.e., merged) clusters;    
    cNames = unique(clusters(:,i));
    cN = length(cNames);                    % number of current clusters (should be N - i + 1)
    
    deltaQ = zeros(nchoosek(cN,2),3);       % potential change in modularity for merging each pair of clusters
    rowQ = 1;
    
    for C1_ind = 1:(cN-1)
        C1 = cNames(C1_ind);    
        currA = A(clusters(:,i)==C1,:);
        MC1 = sum(currA(:)); % total sum of similarities for cluster C1
        
        for C2_ind = (C1_ind+1):cN
            C2 = cNames(C2_ind);            
            currA = A(clusters(:,i)==C2,:);              
            MC2 = sum(currA(:)); % total sum of similarities for cluster C2            
            
            AC1C2 = A(clusters(:,i)==C1, clusters(:,i)==C2);
            MC1C2 = sum(AC1C2(:));              % sum of similarities between C1 and C2

            if MC1C2>0
                deltaQ(rowQ,1) = C1;                
                deltaQ(rowQ,2) = C2;            
                deltaQ(rowQ,3) = 2*(MC1C2/(2*M) - (MC1*MC2)/(4*(M^2)));
                rowQ = rowQ+1;
            end
        end
    end
    deltaQ = deltaQ(1:(rowQ-1),:);
    nextMerge = find(deltaQ(:,3) == max(deltaQ(:,3)),1);    
    C1 = deltaQ(nextMerge,1);
    C2 = deltaQ(nextMerge,2);   
    currD = D(clusters(:,i)==C1,clusters(:,i)==C2);
    tree(rowT,:) = [C1 C2 mean(currD(:))];       % two merged clusters and their average distance
    
    rowT = rowT + 1;
    
    cPrior = clusters(:,i);
    cNew = cPrior;
    cNew(cPrior==C1) = maxC+1;           % assign the items in C1 to the new (merged) cluster
    cNew(cPrior==C2) = maxC+1;           % assign the items in C2 to the new (merged) cluster
    clusters(:,i+1) = cNew;
end

clusters = fliplr(clusters);
for i = 1:N                              % re-labeling clusters to be between 1-N
    c = clusters(:,i);
    cNames = unique(c);
    cN = length(cNames);
    for j = 1:cN
        c(c==cNames(j)) = j;
    end
    clusters(:,i) = c;
end
    
tree(:,[1 2])=sort(tree(:,[1 2]),2);
distVec = zeros(nchoosek(N,2),1);                % for computing the cophenetic correlation
row = 1;
for i = 1:N
    distVec(row:(row+N-i-1)) = D((i+1):end, i);
    row = row+N-i;
end
distVec = distVec';
cc = cophenet(tree,distVec);                     % Matlab built-in function
tree(:,3) = 1:(N-1);                             % Edit tree for better visualization
                                                 % (tree is not used for any computations from this point)