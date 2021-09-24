function [x, fval] = pat_intlinprog_quad(genmat, Q, hw, constthre)
%% implementation of quadratic programming
% main function for finding a small set of SNPs loci for paternity
% inference

% parameters 

% genmat: genotype matrix of a parental population, assuming [0 1 2], each
% individual in the column direction and each marker in the row direction.

% Q: a matrix of the adjacency weight q, assuming distance of pair of loci.
% A form needs to be adopted for intlinprog, and this can be done by
% Qmat_prep_intq.m

% hw: heterozygosity weight h, the power of discrimination for
% heterozygous genotype

% constthre: constraint threshold, how many distinct loci (homoz ref vs
% homoz alt, plus info from heterogyzous loci depending on hw) is needed
% for discriminating a pair of individual. Setting "1" worked in every
% optimization in our test. 
%% 

    [nummar, numind] = size(genmat);
    const_cnt = 1;
    a_id_cnt = 1;
    b_id_cnt = 1;
    M = [0 1 hw; 1 1/hw 1; hw 1 0];
    thre = hw*constthre;
    a1 = zeros(nummar*numind*(numind-1)/2+nummar*(nummar-1)/2*7,1);
    a2 = zeros(nummar*numind*(numind-1)/2+nummar*(nummar-1)/2*7,1);
    aval = zeros(nummar*numind*(numind-1)/2+nummar*(nummar-1)/2*7,1);
    b1 = zeros(numind*(numind-1)/2+nummar*(nummar-1)/2,1);
    b2 = zeros(numind*(numind-1)/2+nummar*(nummar-1)/2,1);
    bval = zeros(numind*(numind-1)/2+nummar*(nummar-1)/2,1);

    for i = 1:numind-1
        for j = 1:numind
            if i < j
                % missingness fill with zero for letting not them being used for
                % the uniqueness calculation
                vi = genmat(:,i);
                vj = genmat(:,j);
                vj(isnan(vi)) = 0;
                vi(isnan(vj)) = 0;
                vj(isnan(vj)) = 0;
                vi(isnan(vi)) = 0;
                
                a1(a_id_cnt:a_id_cnt+nummar-1) = ones(nummar,1)*const_cnt;
                a2(a_id_cnt:a_id_cnt+nummar-1) = 1:nummar;
                aval(a_id_cnt:a_id_cnt+nummar-1) = -hetmask(vi, vj, M);
                b1(b_id_cnt) = const_cnt;
                b2(b_id_cnt) = 1;
                bval(b_id_cnt) = -thre;
                a_id_cnt = a_id_cnt + nummar;
                b_id_cnt = b_id_cnt + 1;
                const_cnt = const_cnt + 1;
            end
        end
    end
    
    y_cnt = nummar + 1;
    for i = 1:length(vi)-1
        for j = 1:length(vi)
            if i < j               
                a1(a_id_cnt) = const_cnt;
                a2(a_id_cnt) = i;
                aval(a_id_cnt) = 1;
                a_id_cnt = a_id_cnt + 1;
                
                a1(a_id_cnt) = const_cnt;
                a2(a_id_cnt) = j;
                aval(a_id_cnt) = 1;
                a_id_cnt = a_id_cnt + 1;
                
                a1(a_id_cnt) = const_cnt;
                a2(a_id_cnt) = y_cnt;
                aval(a_id_cnt) = -1;
                a_id_cnt = a_id_cnt + 1;
                
                b1(b_id_cnt) = const_cnt;
                b2(b_id_cnt) = 1;
                bval(b_id_cnt) = 1;
                b_id_cnt = b_id_cnt + 1;
                const_cnt = const_cnt + 1;
                              
                a1(a_id_cnt) = const_cnt;
                a2(a_id_cnt) = i;
                aval(a_id_cnt) = -1;
                a_id_cnt = a_id_cnt + 1;
                
                a1(a_id_cnt) = const_cnt;
                a2(a_id_cnt) = y_cnt;
                aval(a_id_cnt) = 1;
                a_id_cnt = a_id_cnt + 1;
                
                const_cnt = const_cnt + 1;
                
                a1(a_id_cnt) = const_cnt;
                a2(a_id_cnt) = j;
                aval(a_id_cnt) = -1;
                a_id_cnt = a_id_cnt + 1;
                
                a1(a_id_cnt) = const_cnt;
                a2(a_id_cnt) = y_cnt;
                aval(a_id_cnt) = 1;
                a_id_cnt = a_id_cnt + 1;
                
                const_cnt = const_cnt + 1;
                y_cnt = y_cnt + 1;
            end
        end
    end
       
    intcon = 1:(nummar+nummar*(nummar-1)/2);
    lb = zeros((nummar+nummar*(nummar-1)/2),1);
    ub = ones((nummar+nummar*(nummar-1)/2),1);
    
    Aineq = sparse(a1,a2,aval,numind*(numind-1)/2+ nummar*(nummar-1)/2*3,nummar + nummar*(nummar-1)/2);
    bineq = sparse(b1,b2,bval,numind*(numind-1)/2+nummar*(nummar-1)/2*3,1);
    options = optimoptions("intlinprog","Maxtime",5e4);
    
    tic
    [res,fval,exitflag,output] = intlinprog(Q,intcon,Aineq,bineq,[],[],lb,ub,[],options);
    toc
    x = round(res(1:nummar));

end

function d = hetmask(vi, vj, M)
    d = zeros(length(vi),1);
    for i = 1:1:length(vi)
        d(i,1) = M(vi(i)+1, vj(i)+1);
    end
    %d = o.*sel;
end
