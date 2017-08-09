function [tree, clusters, cc] = runGreedyOptQ(A,D)

%%% Hierarchical clustering based on modularity optimization
%%% INPUT:
%%% A = a symmetric NxN matrix of weighted edges between N nodes
%%%     (or similarities between N items)
%%% D = a symmetric NxN matrix of distances between the same N items
%%% OUTPUT:
%%% tree = a Nx3 matrix, containing information about the structure of 
%%%        the hierarchical tree (type "help linkage" for details)
%%% clusters = Nx(N-1) matrix, each column is a different partition of the N nodes
%%%            obtained by cutting the tree at some level; all possible
%%%            partitions are computed (from each node in its own cluster, through 
%%%            to all nodes in a single cluster)
%%% cc = cophenetic correlation (how faithfully tree captures the original similarities)

%%% Idan Blank, August 8, 2017; EvLab rulz!

A = triu(A,1) + tril(A,-1);                 % remove similarities between each node and itself 
                                            % (twoM below assumes there are no numbers in the diagonal)
M = sum(sum(triu(A,1)));                    % sum of all pairwise similarities 
N = size(A,1);
tree = zeros(N-1,3);                        % similar to what the linkage function returns
rowT = 1;                                   % for filling in tree

C = (1:N)';                                 % initial clustering: each item is assigned to a different cluster
for i = 1:(N-1)
    maxC = max(C);                          % for labeling new (i.e., merged) clusters;    
    cNames = unique(C);
    cN = length(cNames);                    % number of current clusters (should be N - i + 1)
    
    deltaQ = zeros(nchoosek(cN,2),3);       % potential change in modularity for merging each pair of clusters
    rowQ = 1;
    
    for C1_ind = 1:(cN-1)
        C1 = cNames(C1_ind);        
        A_in = A(C==C1,C==C1);              % similarities within cluster C1
        A_out = A(C==C1, ~(C==C1));         % similarities between items in C1 and items out of C1
        MC1 = sum(A_in(:)) + sum(A_out(:)); % total sum of similarities for cluster C1
        
        for C2_ind = (C1_ind+1):cN
            C2 = cNames(C2_ind);            
            A_in = A(C==C2,C==C2);              % similarities within cluster C2
            A_out = A(C==C2, ~(C==C2));         % similarities between items in C1 and items out of C2
            MC2 = sum(A_in(:)) + sum(A_out(:)); % total sum of similarities for cluster C2            
            
            AC1C2 = A(C==C1, C==C2);
            MC1C2 = sum(AC1C2(:));              % sum of similarities between C1 and C2

            deltaQ(rowQ,1) = C1;                
            deltaQ(rowQ,2) = C2;            
            deltaQ(rowQ,3) = 2*(MC1C2/(2*M) - (MC1*MC2)/(4*(M^2)));
            rowQ = rowQ+1;
        end
    end
  
    nextMerge = find(deltaQ(:,3) == max(deltaQ(:,3)),1);    
    C1 = deltaQ(nextMerge,1);
    C2 = deltaQ(nextMerge,2);   
    currD = D(C==C1,C==C2);
    tree(rowT,:) = [C1 C2 max(currD(:))];       % two merged clusters and their average distance
    
    rowT = rowT + 1;
    
    C(C==C1) = maxC+1;                           % assign the items in C1 to the new (merged) cluster
    C(C==C2) = maxC+1;                           % assign the items in C2 to the new (merged) cluster
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
clusters = cluster(tree, 'maxclust', 1:N);       % Matlab built-in function