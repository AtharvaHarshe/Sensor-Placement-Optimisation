%% ============================ FIGURE 9 REPLICATION ============================
fprintf("Computing RMSE curves for LG / PV / UN / RA ...\n");

maxK = 23;                   % as in the paper
T = size(TestSet,1);
nSensors = 54;
numRandom = 20;              % number of random trials for RA curve

RMSE_LG = zeros(maxK,1);
RMSE_PV = zeros(maxK,1);
RMSE_UN = zeros(maxK,1);
RMSE_RA = zeros(maxK,1);

%% Loop over sensor count k
for k = 1:maxK
    fprintf("k = %d\n", k);

    % ---------- LG ----------
    S_LG = Pareto_Sets{k};


    % ---------- Predictive Variance (PV) ----------
    S_PV = pv_select(coords, k);      % predictive variance greedy
  

    % ---------- Uniform Placement (UN) ----------
    S_UN = uniform_select(coords, k); % spread-out sensors
   

    % ---------- Random placement (RA) ----------
    tmp = zeros(numRandom,1);
    for r = 1:numRandom
        S_r = randperm(54, k);
        tmp(r) = rmse_reconstruct(S_r , Phi, TrainSet, TestSet);
    end
    RMSE_RA(k) = mean(tmp);
   RMSE_LG(k) = rmse_reconstruct(S_LG, Phi, TrainSet, TestSet);
RMSE_PV(k) = rmse_reconstruct(S_PV, Phi, TrainSet, TestSet);
RMSE_UN(k) = rmse_reconstruct(S_UN, Phi, TrainSet, TestSet);


end




%% ================== PLOT (IDENTICAL STYLE TO FIGURE 9) ==================
figure; hold on;
plot(1:maxK, RMSE_LG, 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k', ...
     'DisplayName', 'LG');
plot(1:maxK, RMSE_PV, 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'PV');
plot(1:maxK, RMSE_UN, 'b-d', 'LineWidth', 2, 'MarkerFaceColor', 'b', ...
     'DisplayName', 'UN');
plot(1:maxK, RMSE_RA, 'g-^', 'LineWidth', 2, 'MarkerFaceColor', 'g', ...
     'DisplayName', 'RA');

xlabel('Number of sensors', 'FontSize', 12);
ylabel('RMSE', 'FontSize', 12);
title('RMSE vs Sensor Count (Reproduced Figure 9)');
legend('Location','northeast'); grid on;
