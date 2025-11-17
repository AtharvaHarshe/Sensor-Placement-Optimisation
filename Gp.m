
%coords = readmatrix("C:\Users\madhu\OneDrive\Desktop\Textbooks\optimization\Final project\coordinates.txt");
X = coords(:, 2:3);                  % (54×2)
Xnorm = (X - mean(X)) ./ std(X);     % normalize coordinates

% --- Prepare training data ---
Ytrain = Yc;                         % from your POD residuals
Xtrain = Xnorm;

% --- Build quadratic basis (mean function for GP) ---
H = [ones(54,1), Xtrain(:,1), Xtrain(:,2), Xtrain(:,1).^2, ...
     Xtrain(:,1).*Xtrain(:,2), Xtrain(:,2).^2];   % 54×6

% --- Solve for mean coefficients (β) ---
beta = H \ Ytrain;

% --- Compute residuals (remove deterministic quadratic mean) ---
R = Ytrain - H * beta;   % 54×30

% --- Scale residuals for numerical stability ---
R = (R - mean(R(:))) / std(R(:));   % zero-mean, unit-variance

%% --- Define Negative Log Marginal Likelihood ---
nll = @(theta) sum_nll(theta, Xtrain, R);

%% --- Optimize GP hyperparameters (bounded) ---
theta0 = log([1, 1, 0.3, 0.1]);      % initial [ℓ_x, ℓ_y, σ_f, σ_n]
lb = log([1e-2, 1e-2, 1e-3, 1e-3]);  % lower bounds
ub = log([10, 10, 10, 1]);           % upper bounds

opts = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 300, ...
    'TolFun', 1e-4, ...
    'TolX', 1e-4);

fprintf('\n---------------------------------------------\n');
fprintf('  Optimizing GP hyperparameters (with bounds)\n');
fprintf('---------------------------------------------\n');

[theta_opt, fval] = fmincon(@(th) nll(th), theta0, ...
    [], [], [], [], lb, ub, [], opts);

% --- Extract optimized parameters ---
ell_x   = exp(theta_opt(1));
ell_y   = exp(theta_opt(2));
sigma_f = exp(theta_opt(3));
sigma_n = exp(theta_opt(4));

fprintf('\n Optimized GP hyperparameters:\n');
fprintf('   ell_x = %.3f\n   ell_y = %.3f\n   sigma_f = %.3f\n   sigma_n = %.3f\n', ...
    ell_x, ell_y, sigma_f, sigma_n);

% --- Build final covariance matrix ---
dX = pdist2(Xtrain(:,1), Xtrain(:,1));
dY = pdist2(Xtrain(:,2), Xtrain(:,2));

K = (sigma_f^2) * exp(-0.5 * ((dX./ell_x).^2 + (dY./ell_y).^2)) ...
  + (sigma_n^2) * eye(size(Xtrain,1));

%% --- Visualize covariance matrix ---
figure;
imagesc(K);
axis equal tight;
colorbar;
title('Final GP Covariance Matrix (Constrained Optimization)');
xlabel('Sensor index');
ylabel('Sensor index');


%clear beta dY dX fval H  lb nll opts theta_opt theta0 ub  ;

%% ============================================================
%  Local Function Definition (must be at file end)
% ============================================================
function val = sum_nll(theta, X, R)
    % Negative log marginal likelihood across all snapshots
    ell_x   = exp(theta(1));
    ell_y   = exp(theta(2));
    sigma_f = exp(theta(3));
    sigma_n = exp(theta(4));

    % Build kernel matrix
    K = (sigma_f^2) * exp(-0.5 * ( ...
        pdist2(X(:,1), X(:,1)).^2 / ell_x^2 + ...
        pdist2(X(:,2), X(:,2)).^2 / ell_y^2 )) ...
        + (sigma_n^2) * eye(size(X,1));

    % Cholesky decomposition (numerical stability)
    [L, p] = chol(K, 'lower');
    if p > 0
        val = Inf; 
        return;
    end

    N = size(X,1);
    T = size(R,2);
    logdetK = 2 * sum(log(diag(L)));

    % Sum NLL across all temporal snapshots
    val = 0;
    for t = 1:T
        y = R(:,t);
        alpha = L' \ (L \ y);
        val = val + 0.5 * (y' * alpha + logdetK + N * log(2*pi));
    end
end

