%% ================== SECTION 5.4 — Reconstruction RMSE ==================
S_final = Pareto_Sets{knee_k};
% --- Ensure TestSet is T × 54 (snapshots × sensors) ---
if size(TestSet,1) == 54
    TestSet = TestSet';          % Now size = 30 × 54
end

% Training mean (must match POD center)
TrainMean = mean(TrainSet, 2);    % 54 × 1

% Subtract mean from TestSet
TestSet_centered = TestSet - TrainMean';   % 30 × 54

% Use only the POD modes you computed earlier
% Phi = 54 × r
r = size(Phi,2);

% Selected sensor set (from Pareto_Sets{knee_k})
S = S_final;                      % 1 × k

% Extract POD rows corresponding to selected sensors
Phi_S = Phi(S, :);                % k × r

% Compute pseudo-inverse (stable even if rank-deficient)
PhiS_pinv = pinv(Phi_S);          % r × k

% Prepare reconstructed matrix
T = size(TestSet_centered,1);     % number of test snapshots
Ypred_full = zeros(T, 54);        % reconstructions

%% ---- Reconstruct each snapshot using POD coefficients ----
for t = 1:T

    y_true = TestSet_centered(t, :);   % 1 × 54
    y_S    = y_true(S)';               % k × 1 measurements at selected sensors

    % POD coefficient estimate
    a_t = PhiS_pinv * y_S;             % r × 1

    % Reconstruct field:  y_hat = Φ a + mean
    y_hat = Phi * a_t;                 % 54 × 1
    Ypred_full(t,:) = y_hat' + TrainMean';   % add mean back
end

%% ---- Compute RMSE ----
sqErr = (Ypred_full - TestSet).^2;      % same orientation
RMSE = sqrt(mean(sqErr(:)));

fprintf("\n========== Section 5.4 Reconstruction Results ==========\n");
fprintf("Selected sensors: %d\n", length(S));
fprintf("Reconstruction RMSE = %.6f\n", RMSE);
fprintf("==========================================================\n");
%%
%% ================== SECTION 5.4 — RMSE vs k PLOT ==================

maxK = length(Pareto_Sets);   % usually 54
RMSE_k = zeros(maxK, 1);

% Training mean (must match POD centering)
TrainMean = mean(TrainSet, 2);     % 54 × 1

% Fix TestSet orientation
if size(TestSet,1) == 54
    TestSet = TestSet';            % 30 × 54
end

% Center TestSet using TRAIN mean
TestSet_centered = TestSet - TrainMean';


%% Precompute POD shapes
r = size(Phi,2);                   % number of POD modes


%% ---- LOOP: compute RMSE for k = 1 ... maxK ----
for k = 1:maxK
    
    S = Pareto_Sets{k};         % selected sensors for this k
    Phi_S = Phi(S, :);          % k × r

    % Stable pseudo-inverse
    PhiS_pinv = pinv(Phi_S);    % r × k

    Ypred_full = zeros(size(TestSet_centered));  % 30 × 54
    
    % ---- Reconstruct each test snapshot ----
    for t = 1:size(TestSet_centered,1)

        y_true = TestSet_centered(t,:);    % 1 × 54
        y_S    = y_true(S)';               % k × 1

        % POD coefficients
        a_t = PhiS_pinv * y_S;             % r × 1

        % Reconstruct full snapshot and add mean back
        Ypred_full(t,:) = (Phi * a_t)' + TrainMean';
    end

    % ---- Compute RMSE for this k ----
    sqErr = (Ypred_full - TestSet).^2;
    RMSE_k(k) = sqrt(mean(sqErr(:)));

    fprintf("k = %2d --> RMSE = %.4f\n", k, RMSE_k(k));
end


%% ===================== Plot RMSE vs k ========================

figure;
plot(1:maxK, RMSE_k, 'LineWidth', 2);
xlabel('Number of Selected Sensors (k)', 'FontSize', 12);
ylabel('RMSE of Reconstruction', 'FontSize', 12);
title('Section 5.4: Reconstruction Error vs Number of Sensors');
grid on;


