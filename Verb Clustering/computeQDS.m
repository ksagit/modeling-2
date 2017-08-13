function QDS = computeQDS(A,c)

%%% Compute modularity for a given clustering %%%

%%% INPUT:
%%% A = a symmetric NxN matrix of weighted edges between N nodes
%%%     (or similarities between N items); must be between 0 and 1
%%% c = Nx1 vector, community / clusters assignments of the N nodes
%%%     (communities / clusters are denoted by integers)

%%% OUTPUT:
%%% QDS = Nx1 vector of the contributions of each cluster to QDS 
%%%       the final QDS value is the sum of this vector (higher QDS = better clustering) 
%%%       (computed as defined by Chen et al., 2014, IEEE Trans. Comp. Soc. Sys.)

%%% Idan Blank, Aug 08 2017; EvLab rulz! 

M = sum(sum(triu(A,1)));        % sum of all weights in the network
cNames = unique(c);             % cluster labels
nC = length(cNames);            % number of communities / clusters

QDS = zeros(nC,1);
for C1_ind = 1:nC
    C1 = cNames(C1_ind);
    nC1 = sum(c==C1);                             % size of cluster C1
    if nC1 > 1                                    % if cluster is non-singleton
        M_out = sum(sum(A(c==C1,~(c==C1))));          % sum of edges from cluster C1 to other clusters        
        M_in = sum(sum(triu(A(c==C1,c==C1),1)));      % sum of edges within cluster C1
        dC1 = (2*M_in)/(nC1*(nC1-1));             % internal density of C1, original formula from Chen et al
        QDS(C1_ind) = QDS(C1_ind) + (M_in/M)*dC1 - (((2*M_in+M_out)/(2*M))*dC1)^2;
    end
    for C2_ind = 1:nC
        C2 = cNames(C2_ind);
        nC2 = sum(c==C2);
        if ~(C1 == C2)
            ACurr = A(c==C1,c==C2);              % edges between C1 & C2
            M_btw = sum(sum(ACurr));             % sum of edges between C1 & C2
            dC1C2 = M_btw/(nC1*nC2);             % pairwise density between clusters C1 & C2, , original formula from Chen et al
            QDS(C1_ind) = QDS(C1_ind) - ((M_btw)/(2*M))*dC1C2;            
        end
    end
end