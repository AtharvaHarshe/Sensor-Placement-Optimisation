load('train_test_coords.mat');

rowMean = mean(TrainSet, 2);
Yc = TrainSet - rowMean;

% SVD (economy)
[U, S, V] = svd(Yc, 'econ');

% Spatial POD modes
Phi = U;

% Temporal coefficients
A = S * V';

% Mode energies
lambda = diag(S).^2;

clear rowMean 


