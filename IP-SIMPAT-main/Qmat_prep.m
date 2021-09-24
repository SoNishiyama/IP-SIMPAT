function Q = Qmat_prep(nummar, chr, bp, len, qw)
% prepare a triangle matrix Q for the adjacency weight

Qt = triu(ones(nummar),1);
for i = 1:nummar-1
    for j = 1:nummar
        if i < j
            if chr(i) == chr(j)
               Qt(i,j) = (abs(bp(i) - bp(j)) / len(chr(i)))^(-1/qw);
            end
        end
    end
end
Q = sparse(Qt);