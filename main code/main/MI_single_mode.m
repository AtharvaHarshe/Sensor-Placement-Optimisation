function Fi = MI_single_mode(Sigma_i, S)
    
    n = size(Sigma_i, 1);
    U = setdiff(1:n, S);

    if isempty(U)
        Fi = 0; 
        return;
    end

    Sigma_SS = Sigma_i(S, S);
    Sigma_UU = Sigma_i(U, U);
    Sigma_SU = Sigma_i(S, U);
    Sigma_US = Sigma_i(U, S);

    % Conditional covariance
    Sigma_U_given_S = Sigma_UU - Sigma_US * (Sigma_SS \ Sigma_SU);

    % Stable MI
    Fi = 0.5 * ( logdet(Sigma_UU) - logdet(Sigma_U_given_S) );
end

function y = logdet(A)
    % Stable log-determinant using Cholesky
    % Add tiny jitter for numerical stability
    A = A + 1e-12 * eye(size(A));
    U = chol(A);
    y = 2 * sum(log(diag(U)));
end
