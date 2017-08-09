function Q = computeQ(A,c)

%%% Compute modularity for a given clustering %%%
%%% INPUT:
%%% A = a symmetric NxN matrix of weighted edges between N nodes
%%%     (or similarities between N items)
%%% c = Nx1 vector, community / clusters assignments of the N nodes
%%%     (communities / clusters are denoted by integers)
%%% OUTPUT:
%%% Q = modularity (as defined by Newman & Girvan, 2004, Phys. Rev.)
%%%     (higher modularity = better clustering) 

%%% Idan Blank, Aug 08 2017; EvLab rulz! 

A = triu(A,1) + tril(A,-1);     % remove similarities between each node and itself 
                                % (twoM below assumes there are no numbers in the diagonal)
                                % also, the diagonal has a constant contribution to modularity across different partitions
k = sum(A,1);                   % degree of each node
twoM = sum(k);                  % Two M, twice the sum of all weights in the network
cNames = unique(c);             % cluster labels
nC = length(cNames);            % number of communities / clusters
Q = 0;
for cInd = 1:nC
    cCurr = cNames(cInd);
    cNodes = find(c==cCurr);                      % nodes in the current community / cluster    
    ACurr = sum(sum(A(cNodes,cNodes)));           % sum of all pairwise similarities between nodes within current community / cluster
    kCurr = k(cNodes);                            % degrees of nodes within current community / cluster
    kProducts = kCurr'*kCurr;                     % all products between the degrees of every two nodes
    kProducts = triu(kProducts,1) + tril(kProducts,-1);
    kProducts = sum(kProducts(:));                % sum of all products between the degrees of every two nodes
    Q = Q + ACurr - (kProducts/twoM);
end
Q = Q/twoM;