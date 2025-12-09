function F = MI_weighted(Sigma_all, lambda, S)

    r = length(lambda);

    % Normalize eigenvalues to get weights
    w = lambda(:) / (1);

    F = 0;
    for i = 1:r
        Fi = MI_single_mode(Sigma_all{i}, S);
        F = F + w(i) * Fi;       % <--- use normalized weights
    end
    F = F/10000;
end
