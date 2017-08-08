function [c,t] = optimizeModularity(A)

%%% Community detection in a similarity matrix

%%% INPUT:
%%% A = a symmetric NxN matrix of similarities between N items

%%% OUTPUT:
%%% c = Nx3 matrix, each column contains community assignments of the N items by a different method 
%%% t = table of modularity and modularity-density values for each method
%%%     (higher = better)

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

%% Hierarchical clustering with post-hoc modularity optimization %%
[tree, clusters, cc] = runHC(A);
Q = zeros(N,1);
for c = 1:N
    Q(c) = computeQ(A,clusters(c,:));
end
Q = flipud(Q);            % first value is now modularity for the singleton clustering; 
                          % last value is modularity for a single cluster solution
plotHC(tree,Q,cc);

