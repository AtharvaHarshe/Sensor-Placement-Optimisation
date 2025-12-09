

numSensors = size(Sigma_all{1}, 1);
maxSensors = 54;                
Pareto_MI   = zeros(maxSensors,1);
Pareto_Sets = cell(maxSensors,1);
prevMI = -inf;

for eps = 1:maxSensors
    
    S = [];                      % selected sensors start empty
    remaining = 1:numSensors;    % all candidates
    
    currentMI = 0;               % track MI for the greedy process
    
    % ============================================================
    %                TRUE LAZY–GREEDY SELECTION (C.1)
    % ============================================================
    for k = 1:eps
        
        % Step 1: initialize δ(s) = +∞ for all remaining sensors
        m = length(remaining);
        delta   = inf(1, m);        % cached marginal gains
        current = false(1, m);      % flags for whether δ(s) was recomputed
        
        % Step 2: LazyGreedy inner loop
        while true
            
            % Find s* = sensor with largest cached δ value
            [~, idx] = max(delta);
            s_star   = remaining(idx);
            
            % If already recomputed once, accept s*
            if current(idx)
                break;
            end
            
            % Otherwise recompute true marginal gain:
            F_before = MI_weighted(Sigma_all, lambda, S);
            F_after  = MI_weighted(Sigma_all, lambda, [S, s_star]);
            
            delta(idx) = F_after - F_before;   % true δ(s)
            current(idx) = true;               % mark as updated
        end
        
        % Step 3: Add selected sensor s*
        S = [S, s_star];
        
        % Remove it from "remaining"
        remaining(idx) = [];
        
        % Update current MI
        currentMI = MI_weighted(Sigma_all, lambda, S);
    end
    % ============================================================
    %                END OF LAZY GREEDY BLOCK
    % ============================================================
    
    
    % ---------- Non-decreasing MI stopping rule ----------
  % ---------- Combined stopping rule ----------
if eps > 2
    
    %deltaMI     = currentMI/10000 - prevMI/10000;   
    deltaMI     = currentMI - prevMI;   % current marginal gain
          % previous marginal gain

    % Conditions:
    cond_decreaseMI      = (currentMI < prevMI);
    cond_smallGain       = (deltaMI < 0.1);       % your threshold (adjust as needed)
 

    if cond_decreaseMI || cond_smallGain 
        fprintf("\nSTOP at ε = %d because: ", eps);

        if cond_decreaseMI
            fprintf("MI decreased; ");
        end
        if cond_smallGain
            fprintf("marginal gain < threshold; ");
        end
        

        fprintf("\n");

        Pareto_MI   = Pareto_MI(1:eps-1);
        Pareto_Sets = Pareto_Sets(1:eps-1);

        break;
    end
end

% Track MI values for next iteration
prevPrevMI = prevMI;
prevMI      = currentMI;

    
    % Store Pareto point

    Pareto_MI(eps)   = currentMI;
    Pareto_Sets{eps} = S;
    %normilised_MI = currentMI/10000;
    %Pareto_MI_Normalised(eps) = normilised_MI;
    fprintf('ε = %d sensors → MI = %.4f\n', eps, currentMI);
end

%%
x  = (1:eps-1);     % sensor counts
F  = Pareto_MI*10;          % MI values

[~, knee_k, knee_MI, MC] = findKneePoint(x, F);

fprintf('\n*** Knee point found at |S| = %d sensors ***\n', knee_k);
fprintf('Weighted MI at knee = %.4f\n', knee_MI);

clear cond_decreasingGain cond_decreaseMI cond_smallGain current delta;
clear deltaPrevMI deltaMi F_before F_after idx k m prevMI prevPrevMI S s_star;

%%
figure;
hold on;
plot(1:eps-1, Pareto_MI_Normalised, 'o-', 'LineWidth', 2);
stem(knee_k,knee_MI/10000,'s','LineWidth', 4);
xlabel('Number of Sensors (ε)');
ylabel('Weighted Mutual Information');
title('Pareto Frontier (LG–ε–constraint)');
grid on;
hold off;

%%

k = knee_k;
GreedySet = Pareto_Sets{k};
GreedyMI = Pareto_MI(k);

[LB, UB, gap] = BnB_AlgD1(Sigma_all, lambda, GreedySet, k);

%% === Interpretation & Validation Output ===

fprintf("\n==================== BnB VALIDATION ====================\n");
fprintf("Sensor budget (k): %d\n", k);
fprintf("Greedy MI (LB):    %.6f\n", LB);
fprintf("Upper Bound (UB):  %.6f\n", UB);
fprintf("Gap (UB - LB):     %.6f\n", gap);

% Threshold for "good" (you can adjust)
tol = 1e-3 * LB;      % relative threshold (0.1%)

if abs(gap) <= tol
    fprintf("\nSTATUS: ✔ Greedy solution is VALID and essentially OPTIMAL.\n");
    fprintf("        (UB and LB differ by less than %.2f%%)\n", 100 * tol / LB);
elseif gap <= 0.02 * LB
    fprintf("\nSTATUS: ✔ Greedy solution is NEAR-OPTIMAL.\n");
    fprintf("        (Gap ≤ 2%% of MI — acceptable for deployment)\n");
else
    fprintf("\nSTATUS: ⚠ Greedy solution may NOT be optimal.\n");
    fprintf("        Consider running tighter BnB or reducing candidate set.\n");
end

fprintf("=========================================================\n\n");

%%
%% === FIGURE 1: Greedy Pareto Frontier with Knee Point ===

% Assumes you already have these:
%   Pareto_MI   : vector of MI values for k = 1:maxSensors
%   Pareto_Sets : cell array of sets for each k
%   knee_k      : knee-point sensor count (integer)
%   knee_MI     : MI value at knee
%   maxSensors  : usually 54

figure('Color','white'); hold on;

% --- Plot Pareto frontier ---
plot(1:eps-1, Pareto_MI, '-o', ...
    'MarkerSize', 4, 'MarkerFaceColor', [0 0.45 0.74]);


% --- Highlight knee point ---
stem(knee_k, knee_MI/10, 'd', ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.85 0.33 0.10]);  % orange

% --- Labels and title ---
xlabel('Number of Sensors (k)', 'FontSize', 12);
ylabel('Weighted Mutual Information', 'FontSize', 12);
title('Pareto Frontier and Knee Point Selection', 'FontSize', 14);

% --- Grid and formatting ---
grid on;


% --- Annotation text near knee point ---
text(knee_k - 5, knee_MI/10+ 2, sprintf('  Knee point: k = %d', knee_k), ...
    'FontSize', 12, 'Color', [0.85 0.33 0.10]);

% --- Tight axis ---
xlim([1 eps-1]);

