function [S,Q] = voulain_clustering(gamma, sim_matrix) 
    num_trials = 10;
    sz = size(sim_matrix, 1);
    
    A = sim_matrix;
    k = sum(A);
    twom = sum(k);			% sum of all similarities in the A
    B = A - gamma*k'*k/twom;	% Modularity matrix
    S = cell(10,1);
    Q = cell(10,1);
    for i = 1:num_trials
        [s,q] = genlouvain(B);
        S{i} = s;
        Q{i} = q;
    end

    consensus = zeros([sz sz]);
    for i = 1:sz
        for j = i:sz
            total = 0;
            for k = 1:num_trials
                total = total + double(S{k}(i) == S{k}(j));
            end
            f = total / num_trials;
            consensus(i,j) = f;
            consensus(j,i) = f;
        end
    end
    
    A = consensus;
    gamma = 1;			% default value
    k = sum(A);
    twom = sum(k);			% sum of all similarities in the A
    B = A - gamma*k'*k/twom;	% Modularity matrix
    [S,Q] = genlouvain(B);
    Q = Q/twom;

    