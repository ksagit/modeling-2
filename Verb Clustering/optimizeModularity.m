function [c,t] = optimizeModularity(A)

%%% Community detection in a similarity matrix

%%% INPUT:
%%% A = a symmetric NxN matrix of similarities between N items

%%% OUTPUT:
%%% c = Nx5 matrix, each column contains community assignments of the N items by a different method 
%%% t = 1x5 vector, modularity / modularity-density values for each method (higher = better)

%%% COMMUNITY DETECTION METHODS:
%%% HC = hierarchical clustering (hierarchy built regardless of modularity,
%%%      then the partition with the highest modularity is chosen)
%%% optQ = hierarchical clustering built by optimizing modularity
%%%        (then the partition with the highest modularity is chosen)
%%% Louvain = community detection using the Louvain method
%%% optQDS = a partition built in by optimizing density-modularity
%%%          using an iterative split-merge algorithm
%%% optQHEP = a partition obtained by optimizing modularity with a
%%%          HE-prime consolidation ratio heuristic to encourage balanced
%%%          community growth (instead of some huge, some small)

%%% REFERENCES:
%%% Original modularity measure: Newman, M. E., & Girvan, M. (2004). Finding and evaluating community structure in networks. Physical review E, 69(2), 026113.
%%% optQ: Newman, M. E. (2004). Fast algorithm for detecting community structure in networks. Physical review E, 69(6), 066133.
%%% Louvain: Mucha, P. J., Richardson, T., Macon, K., Porter, M. A., & Onnela, J. P. (2010). Community structure in time-dependent, multiscale, and multiplex networks. science, 328(5980), 876-878.
%%%         http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
%%% optQDS: Chen, M., Kuzmin, K., & Szymanski, B. K. (2014). Community detection via maximization of modularity and its variants. IEEE Transactions on Computational Social Systems, 1(1), 46-65.
%%% optQHEP: Wakita, K., & Tsurumi, T. (2007, May). Finding community structure in mega-scale social networks. In Proceedings of the 16th international conference on World Wide Web (pp. 1275-1276). ACM.

%%% Idan Blank, August 8 2017; EvLab Rulz!

N = size(A,1);      % number of nodes / items
D = 1-A;            % NxN matrix of dissimilarities
c = zeros(N,5);
t = zeros(1,5);

%% HC: hierarchical clustering with post-hoc modularity optimization %%
[tree, clusters, cc] = runHC(D);     % first column of clusters: single cluster; last column: singleton clusters
Q = zeros(N,1);
for cInd = 1:N
    Q(cInd) = computeQ(A,clusters(:,cInd));
end

cWin = clusters(:,find(Q==max(Q)));   % clustering with the highest modularity
Q = flipud(Q);                        % first value is now modularity for the singleton clustering; 
                                      % last value is modularity for a single cluster solution                                      
plotHC(tree,Q,cWin,cc);
c(:,1) = cWin;
t(1) = max(Q);

%% optQ: hierarchical clustering via by modularity optimization %%
[tree, clusters, cc] = runGreedyOptQ(A,D);
Q = zeros(N,1);
for cInd = 1:N
    Q(cInd) = computeQ(A,clusters(:,cInd));
end
cWin = clusters(:,find(Q==max(Q)));   % clustering with the highest modularity
Q = flipud(Q);                        % first value is now modularity for the singleton clustering; 
% plotHC(tree,Q,cWin,cc);             % this often gives a tangled tree, no need to plot
c(:,2) = cWin;
t(2) = max(Q);

%% Louvain method %%
addpath(genpath('/Users/iblank/Desktop/MIT/Experiments/DynamicNetworksTools/GenLouvain2.0/'));
nIter = 100;         % number of iterations (to take into account the random initalization of the method)
gamma = 1;
C = zeros(N,nIter);
Q = zeros(nIter,1);
for i = 1:nIter
    gamma = 1;                  % default value
    k = sum(A);
    twoM = sum(k);              % sum of all similarities in the A
    B = A - gamma*k'*k/twoM;	% Modularity matrix
    [currC, currQ] = genlouvain(B);
    C(:,i) = currC;
    Q(i) = currQ/twoM;
end

consensusMat = ones(N,N);      % consensus similarity matrix
for i = 1:(N-1)
    for j = (i+1):N
        consensusMat(i,j) = sum(C(i,:)==C(j,:))/nIter;
        consensusMat(j,i) = consensusMat(i,j);
    end
end
k = sum(consensusMat,1);
twoM = sum(k);                      % sum of all similarities in the A
B = consensusMat - gamma*k'*k/twoM;	% Modularity matrix
[C, ~] = genlouvain(B);
c(:,3) = C;
t(3) = computeQ(A,C);

% QDS = zeros(N,1);
% QDS(c) = computeQDS(A,clusters(:,c));
% QDS = flipud(QDS);