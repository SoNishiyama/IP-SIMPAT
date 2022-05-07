function [x, fval, d] = pat_ga(genmat, D, hw, constthre, popsize, crossfrac, elitec, x0)
%% implemetation of genetic algorithm
% implementation of genetic algorithm for finding an optimized set of SNPs
% loci for paternity inference.

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

% popsize: population size for ga function
% crossfrac: crossover flaction for ga function
% elitc: elite count for ga function
%% 
    [nummar, numind] = size(genmat);
    c = 1;
    M = [0 1 hw; 1 1/hw 1; hw 1 0];
    thre = hw*constthre;
    Aineq = zeros(numind*(numind-1)/2,nummar);
    bineq = -ones(numind*(numind-1)/2,1)*thre;
    intcon = 1:nummar;
    lb = zeros(nummar,1);
    ub = ones(nummar,1);
    
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
                
                % adding hw weight for heterogyzous genotype 
                % i.e. 8 adds 1/8 weight for [0 1] pair compared with [0 2] pair 
                Aineq(c,:) = -hetmask(vi, vj, M);
                c = c + 1;
                if rem(c,10000) == 0
                    c
                end
            end
        end
    end
    
    %x0 = ones(1, nummar);
    options = optimoptions('ga','UseParallel',true, 'PlotFcn',{@gaplotbestf,@gaplotmaxconstr}, ...
        "PopulationSize", popsize, "CrossoverFraction", crossfrac, "EliteCount", elitec, "MaxGenerations",1000, ...
        'InitialPopulationRange', [0;1], 'InitialPopulationMatrix', x0);
    f = @(x)sum(sum(D.*(x'*x)));
    [x,fval] = ga(f,nummar, Aineq,bineq,[],[],lb, ub, [],intcon,options);
    d = eval_depth(genmat, M, find(x));
    
end

function d = hetmask(vi, vj, M)
    d = zeros(length(vi),1);
    for i = 1:1:length(vi)
        d(i,1) = M(vi(i)+1, vj(i)+1);
    end
    %d = o.*sel;
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