function [tree, clusters, cc] = runHC(A)

%%% Hierarchical clustering based on average linkage
%%% INPUT:
%%% A = a symmetric NxN matrix of weighted edges between N nodes
%%%     (or similarities between N items)
%%% OUTPUT:
%%% tree = a Nx3 matrix, containing information about the structure of 
%%%        the hierarchical tree (type "help linkage" for details)
%%% clusters = Nx(N-1) matrix, each column is a different partition of the N nodes
%%%            obtained by cutting the tree at some level; all possible
%%%            partitions are computed (from each node in its own cluster, through 
%%%            to all nodes in a single cluster)
%%% cc = cophonetic correlation (how faithfully tree captures the original similarities)

%%% Idan Blank, August 8, 2017; EvLab rulz!

D = max(A(:))-A+min(A(:));                  % convert similarities to dissimilarities
N = size(D,1);
distVec = zeros(nchoosek(N,2),1);           % convert dissimilarities to a vector like those given by the pdist function
row = 1;
for i = 1:N
    distVec(row:(row+N-i-1)) = D((i+1):end, i);
    row = row+N-i;
end
distVec = distVec';
tree = linkage(distVec,'average');          % Matlab built-in function 
cc = cophenet(tree,distVec);                % Matlab built-in function
clusters = cluster(tree, 'maxclust', 1:N);  % Matlab built-in function