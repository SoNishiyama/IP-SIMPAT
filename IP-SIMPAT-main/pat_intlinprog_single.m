function [x, fval,tol] = pat_intlinprog_single(genmat, hw, constthre)
%% test function for finding a small set of SNPs loci for paternity inference 
% This is analogous to paternity_QPintlinprog.m but without distance weight
% parameters 

% genmat: genotype matrix of a parental population, assuming [0 1 2], each
% individual in the column direction and each marker in the row direction.

% hw: heterozygosity weight h, the power of discrimination for
% heterozygous genotype

% constthre: constraint threshold, how many distinct loci (homoz ref vs
% homoz alt, plus info from heterogyzous loci depending on hw) is needed
% for discriminating a pair of individual. Setting "1" worked in every
% optimization in our test. 

%% 
    [nummar, numind] = size(genmat);
    f = ones(nummar,1);
    const_cnt = 1;
    a_id_cnt = 1;
    b_id_cnt = 1;
    M = [0 1 hw; 1 1/hw 1; hw 1 0];
    thre = hw*constthre;
    a1 = zeros(nummar*numind*(numind-1)/2,1);
    a2 = zeros(nummar*numind*(numind-1)/2,1);
    aval = zeros(nummar*numind*(numind-1)/2,1);
    b1 = zeros(numind*(numind-1)/2,1);
    b2 = zeros(numind*(numind-1)/2,1);
    bval = zeros(numind*(numind-1)/2,1);

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
       
    intcon = 1:nummar;
    lb = zeros(nummar,1);
    ub = ones(nummar,1);
    
    Aineq = sparse(a1,a2,aval,numind*(numind-1)/2,nummar);
    bineq = sparse(b1,b2,bval,numind*(numind-1)/2,1);
    options = optimoptions("intlinprog","Maxtime",5e5);
    
    tic
    [res,fval,exitflag,output] = intlinprog(f,intcon,Aineq,bineq,[],[],lb,ub,[],options);
    toc
    x = round(res(1:nummar));
    tol = eval_tolerance(genmat, M, find(x));

end

function d = hetmask(vi, vj, M)
    d = zeros(length(vi),1);
    for i = 1:1:length(vi)
        d(i,1) = M(vi(i)+1, vj(i)+1);
    end
    %d = o.*sel;
end

function res = eval_tolerance(gen, M, target_mar_list)
    [nummar, numind] = size(gen); 
    effect = zeros(numind*(numind-1)/2,1);
    for m = 1:length(target_mar_list)
        c = 1;
        for k = 1:(numind - 1)
            for j = 1:numind
                if k < j
                    gen_k = gen(target_mar_list(m), k);
                    gen_j = gen(target_mar_list(m), j);
                    if sum(isnan([gen_k gen_j])) >= 1
                        gen_k = 0;
                        gen_j = 0;
                    end
                    tol = M(gen_k+1, gen_j+1);
                    effect(c) = effect(c) + tol;
                    c = c + 1;
                end
            end
        end
    end
    res = median(effect);
end