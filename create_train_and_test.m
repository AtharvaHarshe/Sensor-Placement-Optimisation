% -------------------------------------------
% Create Train/Test snapshot sets (keep all 54 sensors)
% Assumes:
%   - Y is 54 × m snapshot matrix already in workspace
%   - m ≥ 50
% -------------------------------------------

% Pick 50 distinct snapshot indices
snap_idx = randperm(size(Y, 2), 100);

% Train = first 30 snapshots
train_idx = snap_idx(1:70);

% Test = last 20 snapshots
test_idx  = snap_idx(71:100);

% Create Train and Test matrices
TrainSet = Y(:, train_idx);   % 54 × 30
TestSet  = Y(:, test_idx);    % 54 × 20



clear snap_idx test_idx train_idx