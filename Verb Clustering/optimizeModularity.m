function [c,t] = optimizeModularity(A)

%%% Community detection in a similarity matrix

%%% INPUT:
%%% A = a symmetric NxN matrix of similarities between N items

%%% OUTPUT:
%%% c = Nx4 matrix, each column contains community assignments of the N items by a different method 
%%% t = 1x4 vector, modularity / modularity-density values for each method (higher = better)

%%% COMMUNITY DETECTION METHODS:
%%% HC = hierarchical clustering (hierarchy built regardless of modularity,
%%%      then the partition with the highest modularity is chosen)
%%% optQ = hierarchical clustering built by optimizing modularity
%%%        (then the partition with the highest modularity is chosen)
%%% Louvain = community detection using the Louvain method
%%% optQDS = a partition built in by optimizing density-modularity
%%%          using an iterative split-merge algorithm

%%% REFERENCES:
%%% Original modularity measure: Newman, M. E., & Girvan, M. (2004). Finding and evaluating community structure in networks. Physical review E, 69(2), 026113.
%%% optQ: Newman, M. E. (2004). Fast algorithm for detecting community structure in networks. Physical review E, 69(6), 066133.
%%% Louvain: Mucha, P. J., Richardson, T., Macon, K., Porter, M. A., & Onnela, J. P. (2010). Community structure in time-dependent, multiscale, and multiplex networks. science, 328(5980), 876-878.
%%%         http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
%%% optQDS: Chen, M., Kuzmin, K., & Szymanski, B. K. (2014). Community detection via maximization of modularity and its variants. IEEE Transactions on Computational Social Systems, 1(1), 46-65.

%%% May be added in the future:
%%% optQHEP = a partition obtained by optimizing modularity with a
%%%          HE-prime consolidation ratio heuristic to encourage balanced
%%%          community growth (instead of some huge, some small)
%%% Wakita, K., & Tsurumi, T. (2007, May). Finding community structure in mega-scale social networks. In Proceedings of the 16th international conference on World Wide Web (pp. 1275-1276). ACM.

%%% Idan Blank, August 8 2017; EvLab Rulz!

N = size(A,1);                  % number of nodes / items
A = A-min(A(:));                % Shift A to have a min of 0
A = A/max(A(:));                % Scale A to have a max of 1
A0 = triu(A,1) + tril(A,-1);    % A with zeros along the diagonal (A as a graph with no self-loops)
A1 = A0 + eye(N);               % As with 1s along the diagonal (A as a similarity matrix)

D = 1-A1;                       % NxN matrix of dissimilarities
c = zeros(N,4);
t = zeros(1,4);

doHC = 1;
doQ = 1;
doLouvain = 1;
doQDS = 1;

%% HC: hierarchical clustering with post-hoc modularity optimization %%
if doHC
    [tree, clusters, cc] = runHC(D);      % first column of clusters: single cluster; last column: singleton clusters
    Q = zeros(N,1);
    for cInd = 1:N
        Q(cInd) = computeQ(A,clusters(:,cInd));
    end

    cWin = clusters(:,Q==max(Q));   % clustering with the highest modularity
    Q = flipud(Q);                        % first value is now modularity for the singleton clustering; 
                                          % last value is modularity for a single cluster solution                                      
    plotHC(tree,Q,cWin,cc);
    c(:,1) = cWin;
    t(1) = max(Q);
end

%% optQ: hierarchical clustering via by modularity optimization %%
if doQ
    [tree, clusters, cc] = runGreedyOptQ(A0,D);
    Q = zeros(N,1);
    for cInd = 1:N
        Q(cInd) = computeQ(A0,clusters(:,cInd));
    end
    cWin = clusters(:,Q==max(Q));       % clustering with the highest modularity
    Q = flipud(Q);                      % first value is now modularity for the singleton clustering; 
    plotHC(tree,Q,cWin,cc);             % this often gives a tangled tree, no need to plot
    c(:,2) = cWin;
    t(2) = max(Q);
end

%% Louvain method %%
if doLouvain
    addpath(genpath('/Users/iblank/Desktop/MIT/Experiments/DynamicNetworksTools/GenLouvain2.0/'));
    nIter = 100;         % number of iterations (to take into account the random initalization of the method)
    gamma = 1;
    C = zeros(N,nIter);
    Q = zeros(nIter,1);
    for i = 1:nIter
        gamma = 1;                  % default value
        k = sum(A0);
        twoM = sum(k);              % sum of all similarities in the A
        B = A0 - gamma*k'*k/twoM;	% Modularity matrix
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
end

%% optQDS: clustering via optimizing modularity-density %%
if doQDS
    C = ones(N,1);   % all items in a single cluster
    cN = 1;          % number of clusters
    cN_split = 0;
    cN_merge = 0;
    QDS = 0;
    while ~(cN==cN_split && cN==cN_merge)
        cN = length(unique(C));
        C = optQDS_split(A0,C);
        cN_split = length(unique(C));
        C = optQDS_merge(A0,C);
        cN_merge = length(unique(C));
        QDS_new = computeQDS(A0,C);
        QDS_new = sum(QDS_new);
        if QDS_new == QDS
            break
        else
            QDS = QDS_new;
        end
    end
    c(:,4) = C;
    t(4) = QDS_new;
end