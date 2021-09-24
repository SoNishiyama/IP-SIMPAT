function Q = Qmat_prep_intq(nummar, chr, chrlen, bp, qw)
% prepare an adjacency matrix Q for QP, adopted for intlinprog

d1 = zeros(nummar,1);
d2 = ones(nummar*(nummar-1)/2,1);
Q = cat(1,d1,d2);
cnt = nummar + 1;
for i = 1:nummar-1
    for j = 1:nummar
        if i < j
            if chr(i) == chr(j)
                Q(cnt,1) = (abs(bp(i) - bp(j)) / chrlen.data(chr(i),2))^(-1/qw);
            end
            cnt = cnt + 1;
        end
    end
end