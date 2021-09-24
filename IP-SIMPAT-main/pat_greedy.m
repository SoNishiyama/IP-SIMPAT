function [x, fval, d, iter, fval_adj] = pat_greedy(genmat, Q, hw, constthre)
%% Implementation of greedy algorhithm
%implementation of greedy algorithm for finding an optimized set of SNPs
%loci for paternity inference.

% An individual in the first column of the genmat is considered as maternal
% parent.


% parameters 


% genmat: genotype matrix of a parental population, assuming [0 1 2], each
% individual in the column direction and each marker in the row direction.

% Q: a triangle matrix of the adjacency weight q, assuming distance of pair
% of loci 

% hw: heterozygosity weight h, the power of discrimination for
% heterozygous genotype

% constthre: constraint threshold, how many distinct loci (homoz ref vs
% homoz alt, plus info from heterogyzous loci depending on hw) is needed
% for discriminating a pair of individual. Setting "1" worked in every
% optimization in our test. 

%% 
    [nummar, numind] = size(genmat);
    curr = 0;
    
    M = [0 1 hw; 1 1/hw 1; hw 1 0];
    thre = hw*constthre;
    
    X = ones(numind)*thre;
    X = triu(X,1);
    
    res_array = [];
    
    while 1
        temp_array = zeros(nummar,1);
        for i = 1:nummar
            if sum(ismember(res_array,i)) == 0
                weight_curr = weight_d(Q, res_array, i);
                temp_array(i) = eval(genmat, M, X, weight_curr, i);
            end
        end
        
        eval_max = max(temp_array);
        eval_max_array = find(temp_array==eval_max);
        curr_res = eval_max_array(randsample(length(eval_max_array),1));
        
        res_array = [res_array curr_res];
        X = update_uniq_mat(genmat, M, X, curr_res);
        
        curr = curr + 1;
        disp("marker selection " + curr + ", " + sum(sum(X)) + " uniqueness remains")
        
        if sum(sum(X)) == 0
            break
        end
    end
    x = zeros(nummar,1);
    for i = 1:length(res_array)
        x(res_array(i)) = 1;
    end
    fval = x'*Q*x;
    fval_adj = x'*Q*x - ((sum(x)-1)*(sum(x)-2)/2-1);
    d = eval_depth(genmat, M, find(x));
    iter = curr;
end

function res = weight_d(Q, res_mar_array, target_mar_id)
    dist_sum = 0;
    for i = 1:length(res_mar_array)
        if res_mar_array(i) < target_mar_id
            d = Q(res_mar_array(i), target_mar_id);
        else
            d = Q(target_mar_id, res_mar_array(i));
        end
        if d-0 ~= 1
            dist_sum = dist_sum + d;
        end
    end
    if dist_sum == 0
        dist_sum = 1;
    end
    res = 1/dist_sum;
end

function effect = eval(gen, M, X, weight_d, target_mar_id)
    [nummar, numind] = size(gen);
    eval_vec = eval_vector_prep(gen, M, target_mar_id);
    effect = 0;
    c = 1;
    for k = 1:(numind - 1)
        for j = 1:numind
            if k < j
                if X(k, j) ~= 0
                    dist = min([eval_vec(c), X(k, j)]);
                    effect = effect + weight_d*dist;
                end
                c = c + 1;
            end
        end
    end
end

function res = eval_vector_prep(gen, M, target_mar_id)
    [nummar, numind] = size(gen);
    c2 = 1;
    res = zeros(numind*(numind-1)/2,1);
    for k = 1:(numind - 1)
        for j = 1:numind
            if k < j
                gen_k = gen(target_mar_id, k);
                gen_j = gen(target_mar_id, j);
                if sum(isnan([gen_k gen_j])) >= 1
                    gen_k = 0;
                    gen_j = 0;
                end
                res(c2) = M(gen_k+1,gen_j+1);
                c2 = c2 + 1;
            end
        end
    end
end

function X = update_uniq_mat(gen, M, X, curr_res)
    [nummar, numind] = size(gen);
    for k = 1:(numind - 1)
        for j = 1:numind
            if k < j
                if X(k, j) ~= 0
                    gen_k = gen(curr_res, k);
                    gen_j = gen(curr_res, j);
                    if sum(isnan([gen_k gen_j])) >= 1
                        gen_k = 0;
                        gen_j = 0;
                    end
                    dist = min([M(gen_k+1, gen_j+1) X(k, j)]);
                    X(k, j) = X(k, j) - dist;
                end
            end
        end
    end
end

function res = eval_depth(gen, M, target_mar_list)
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
