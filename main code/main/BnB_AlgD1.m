function [LB, UB, gap] = BnB_AlgD1(Sigma_all, lambda, GreedySet, k)

    % ============ Step 1: compute lower bound ====================
    F_greedy = MI_weighted(Sigma_all, lambda, GreedySet);
    LB = F_greedy;

    % ============ Step 2: compute singleton gains (submodular UB) ============
    n = size(Sigma_all{1},1);
    MG = zeros(1,n);

    for s = 1:n
        MG(s) = fast_singleton_MI(Sigma_all, lambda, s);
    end

    % Sort largest to smallest
    MG = sort(MG, 'descend');

    remaining = k - numel(GreedySet);
    if remaining < 0
        error("Greedy set larger than k.");
    elseif remaining == 0
        UB = LB;
        gap = UB - LB;
        return;
    end

    % ============ Step 3: Compute UB using Algorithm D.1 approximation ============
    alpha = 1 - 1/exp(1);   % (1 - 1/e)
    greedy_residual_gain = sum( MG(1:remaining) );  % submodular upper bound

    UB = LB + greedy_residual_gain / alpha;

    % ============ Step 4: compute optimality gap ================================
    gap = UB - LB;

    fprintf("Greedy MI = %.6f\n", LB);
    fprintf("Alg D.1 UB = %.6f\n", UB);
    fprintf("Optimality gap ≤ %.6f\n", gap);

end


% ========== Helper: singleton MI = upper bound on true marginal gain ========
function mg = fast_singleton_MI(Sigma_all, lambda, s)
    mg = 0;
    for i = 1:length(lambda)
        sigma_ss = Sigma_all{i}(s, s);
        mg = mg + lambda(i) * 0.5 * log(1 + sigma_ss);
    end
end
