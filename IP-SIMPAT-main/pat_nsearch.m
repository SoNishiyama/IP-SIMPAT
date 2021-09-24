function [res_array, fval, d, iter, fval_adj] = pat_nsearch(genmat, Q, hw, constthre, init, flip_frac)
%% Implementation of neighborhood search for optimizing marker set for paternity inference.

% An individual in the first column of the input.genmat is considered as
% mother and thus not be used for the calculation of pairwise distance.

% In neighborhood search, pairs of the markers in initial list will be
% chosen, and the choice (select/unselect) will be flipped. Subsequently,
% the flipped solution will be evaluated again to have better set of
% markers. Here, for all the pairs, one would be selected from the
% incumbent solution and the other from the unselected set. Iteration will
% be terminated when any change is made on fval, num-markers, and depth. In
% order to reduce the computation time, the flip will be done only on
% "flip-frac" fraction of loci that are correlated to the genotype of a
% given marker locus. During the neighbor search, the set of markers from
% which one marker is excluded from incumbent solution will also be
% investigated.


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

% init: vector for the initial marker set
% flip_frac: fraction v of markers to be flipped in neighborhood search

%% 
    [nummar, numind] = size(genmat);    
    M = [0 1 hw; 1 1/hw 1; hw 1 0];
    
    check_mat = check_mat_prep(genmat, M, flip_frac);
    
    curr_fval_optim = eval_optim(init, Q);
    curr_fval_d = eval_depth(genmat, M, find(init));
    
    iter = 0;
    disp("iteration " + iter + ", fval " + curr_fval_optim + ", " + sum(init) + " markers, depth " + curr_fval_d)
    
    res_array = init;
    while 1
        [best_fval_optim, best_fval_d, best_mar] = neighbor_search(genmat, Q, M, res_array, hw, constthre, check_mat);
    
        if curr_fval_optim == best_fval_optim
            if curr_fval_d == best_fval_d
                iter = 1 + iter;
                disp("iteration " + iter + ", fval " + curr_fval_optim + ", " + sum(res_array) + " markers, depth " + curr_fval_d)
                break
            end
        end
        curr_fval_optim = best_fval_optim;
        curr_fval_d = best_fval_d;
        res_array = zeros(nummar,1);
        res_array(best_mar) = 1;
        iter = iter + 1;
        disp("iteration "+iter+", fval "+curr_fval_optim+", "+sum(res_array)+" markers, depth "+curr_fval_d)
    end
    
    fval = res_array'*Q*res_array;
    fval_adj = res_array'*Q*res_array - ((sum(res_array)-1)*(sum(res_array)-2)/2-1);
    d = best_fval_d;
end


function res = corr_mat_prep(gen, M)
    [nummar, numind] = size(gen); 
    evals_mat = zeros(numind*(numind-1)/2, nummar);
    for i = 1:nummar
        evals_mat(:,i) = eval_vector_prep(gen, M, i);
    end
    res = corr(evals_mat);
end


function res = check_mat_prep(gen, M, swap_frac)
    nummar = size(gen,1);
    corr_mat = corr_mat_prep(gen, M);
    corr_rank_mat = -tiedrank(corr_mat);
    swap_count = round(swap_frac*nummar);
    if swap_count >= nummar
        swap_count = nummar - 1;
    end
    res = zeros(swap_count, nummar);
    for i = 1:nummar
        [a,temp] = sort(corr_rank_mat(:,i));
        res(:,i) = temp(2:swap_count+1);
    end
end


function res = eval_vector_prep(gen, M, target_mar_id)
    numind = size(gen,2);
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

function res = eval_optim(x, Q)
    res = x'*Q*x;
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


function res = combination(curr_list_sparse, check_mat)
    % produce set of markers to be flipped
    res = zeros(length(curr_list_sparse)*size(check_mat,1),2);
    c = 1;
    for i = 1:length(curr_list_sparse)
        swap_list = check_mat(:,curr_list_sparse(i));
        swap_list = swap_list(ismember(swap_list,curr_list_sparse)==0);
        for j = 1:length(swap_list)
            res(c,:) = [curr_list_sparse(i) swap_list(j)];
            c = c + 1;
        end
    end
    res(sum(res,2)==0,:) = [];
end

function res = neighbor_gen(curr_list, comb)
    count1 = size(comb,1);
    res = zeros(count1+sum(curr_list),length(curr_list));
    for i = 1:size(comb,1)
        temp = curr_list;
        for j = 1:size(comb,2)
            gen_j = temp(comb(i,j));
            if gen_j == 0
               temp(comb(i,j)) = 1;
            elseif gen_j == 1
               temp(comb(i,j)) = 0;
            end
        end
        res(i,:) = temp;
    end
    curr_list_sparse = find(curr_list);
    for i = 1:sum(curr_list)
        temp = curr_list;
        temp(curr_list_sparse(i)) = 0;
        res(count1+i,:) = temp;
    end
end


function res = distance(vi, vj, M)
    if sum(isnan(vi)) + sum(isnan(vj)) >= 1
        vj(isnan(vi)) = 0;
        vi(isnan(vj)) = 0;
        vj(isnan(vj)) = 0;
        vi(isnan(vi)) = 0;
    end
    res = 0;
    for i = 1:length(vi)
        res = res + M(vi(i)+1,vj(i)+1);
    end
end
    

function flag = eval_uniqueness(gen, sparse, M, thre, numind)
    flag = 0;
    for i = 1:(numind-1)
        if flag == 1
            break
        end
        for j = 1:numind
            if i < j
                vi = gen(sparse,i);
                vj = gen(sparse,j);
                if distance(vj, vi, M) < thre
                    flag = 1;
                    break
                end
            end
        end
    end
end


function [best_fval_optim, best_fval_d, best_mar] = neighbor_search(gen, Q, M, curr_list, hw, constthre, check_mat)
    numind = size(gen,2);    
    curr_list_sparse = find(curr_list);

    comb = combination(curr_list_sparse, check_mat);

    best_fval_optim = eval_optim(curr_list, Q);
    best_fval_d = eval_depth(gen, M, curr_list_sparse);

    res_neighbors = curr_list_sparse;

    neighbors = neighbor_gen(curr_list, comb);
    num_neighbors = size(neighbors,1);
    c_uniq = 0;
    h = waitbar(0,'Please wait...'); 
    for n = 1:num_neighbors
        flag = eval_uniqueness(gen, find(neighbors(n,:)), M, hw*constthre, numind);
        if rem(n,100) == 0
            h = waitbar(n/num_neighbors,h, "neighbor search "+n+"/"+num_neighbors+", "+c_uniq+" fulfill constraints");
        end
        if flag == 0
            c_uniq = c_uniq + 1;
            curr_neighbor = neighbors(n,:)';
            temp_optim = eval_optim(curr_neighbor, Q);
            if temp_optim <= best_fval_optim
                if temp_optim < best_fval_optim
                    res_neighbors = find(curr_neighbor);
                    best_fval_optim = temp_optim;
                    best_fval_d = eval_depth(gen, M, find(curr_neighbor));
                    continue
                end
                temp_d = eval_depth(gen, M, find(curr_neighbor));
                if temp_d >= best_fval_d
                    if temp_d > best_fval_d
                        res_neighbors = find(curr_neighbor);
                        best_fval_optim = temp_optim;
                        best_fval_d = temp_d;
                        continue
                    end
                    res_neighbors = [res_neighbors find(curr_neighbor)];
                end
            end
        end
    end
    close(h);
    best_mar = res_neighbors(:,randsample(1:size(res_neighbors,2),1));
end                    
