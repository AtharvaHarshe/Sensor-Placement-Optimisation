%% ============================================================
%   Weighted Mutual Information (WMI) + Greedy Sensor Selection
%   Uses: K or K_all, lambda, sigma_n (noise std)
% ============================================================

% ---- Prepare covariance set and weights
Ks = {};
if exist('K_all','var') && iscell(K_all) && ~isempty(K_all)
    Ks = K_all;                                % per-mode covariances (cell of 54x54)
elseif exist('K','var') && ~isempty(K)
    Ks = {K};                                  % single covariance fallback
else
    error('No covariance matrix found. Provide K or K_all.');
end

modes = numel(Ks);

% Weights from POD energies if available (normalize to sum=1)
if exist('lambda','var') && ~isempty(lambda)
    w = lambda(1:min(modes,numel(lambda)));
    w = w(:);  w = w / sum(w);
else
    w = ones(modes,1) / modes;                 % equal weights
end

% Noise variance (sigma_n^2). If not in workspace, use a modest default.
if exist('sigma_n','var') && ~isempty(sigma_n)
    noise_var = sigma_n^2;
else
    noise_var = 0.01;                          % fallback if sigma_n not available
end

% ---- Greedy selection parameters
k_select = 54;                                  % choose how many sensors to place
N = size(Ks{1},1);
cand = 1:N;
S = [];                                         % selected set
F_vals = zeros(k_select+1,1);                   % MI curve; F(0)=0

% Precompute tiny jitter for stability
jitter = 1e-10;

% ---- Helper: logdet(I + (1/noise_var)*K(S,S)) via Cholesky
logdet_gain = @(Ksub) local_logdet( eye(size(Ksub)) + (1/noise_var) * (Ksub + jitter*eye(size(Ksub))) );

% ---- Greedy loop
for t = 1:k_select
    best_gain = -inf;
    best_j = NaN;

    % Current base value (sum across modes)
    base_val = 0;
    if isempty(S)
        base_val = 0;
    else
        for i = 1:modes
            KSS = Ks{i}(S,S);
            base_val = base_val + 0.5 * w(i) * logdet_gain(KSS);
        end
    end

    for j = cand
        Sj = [S, j];
        cur_val = 0;
        for i = 1:modes
            KSSj = Ks{i}(Sj,Sj);
            cur_val = cur_val + 0.5 * w(i) * logdet_gain(KSSj);
        end
        mgain = cur_val - base_val;
        if mgain > best_gain
            best_gain = mgain;
            best_j = j;
        end
    end

    % Commit winner
    S = [S, best_j];
    cand(cand==best_j) = [];
    F_vals(t+1) = F_vals(t) + best_gain;

    fprintf('Step %2d: picked sensor %2d  |  ΔF=%.6f  |  F=%.6f\n', t, best_j, best_gain, F_vals(t+1));
end

fprintf('\nSelected sensors (k=%d):\n', k_select);
disp(S);

% ---- Plot MI curve
figure; plot(0:k_select, F_vals, 'o-','LineWidth',1.5);
xlabel('Number of selected sensors'); ylabel('Weighted MI F(S)');
title('Greedy WMI Sensor Selection'); grid on;

%% ===== Local function (kept at end of script) =====
function val = local_logdet(A)
    % Numerically-stable log(det(A)) using Cholesky (A must be SPD)
    [L,p] = chol(A,'lower');
    if p>0
        % Add a tad more jitter if needed
        epsj = 1e-12;
        [L,p] = chol(A + epsj*eye(size(A)),'lower');
        if p>0, error('Matrix not SPD even after jitter.'); end
    end
    val = 2*sum(log(diag(L)));
end
