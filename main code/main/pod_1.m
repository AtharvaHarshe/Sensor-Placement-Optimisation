%clear all; close all; clc;
%load('train_test_coords.mat');


% find mean and normalise it
meanTemp = mean(TrainSet, 2);       
Yc = TrainSet - meanTemp;        

%compute covariance matrix
R = Yc * Yc.';

% Eigen decomposition -----
[Phi, D] = eig(R, 'vector');      % D = vector of eigenvalues (unsorted)
[lambda, idx] = sort(D, 'descend');  % Sort by descending eigenvalue magnitude
Phi = Phi(:, idx);                % Reorder eigenvectors accordingly

% Normalize POD modes -----
for k = 1:length(lambda)
    Phi(:,k) = Phi(:,k) / norm(Phi(:,k));   % ensure unit norm
end

% Compute temporal coefficients (POD amplitudes) -----
A = Phi.' * Yc;                % r x 70 (each row = temporal evolution of mode k)

% Energy content (pseudo kinetic energy) -----
E_total = sum(lambda);
E_percent = (lambda / E_total) * 100;   % % energy per mode

%  Plot eigenvalue spectrum -----
figure;
stem(E_percent);
xlabel('Mode number'); ylabel('Energy (%)');
title('POD Mode Energy Distribution');

clear idx k meanTemp A D E_percent E_total Yc R 
%%
modes_7 = sum(E_percent(1:7))
modes_12 = sum(E_percent(1:12))
modes_23 = sum(E_percent(1:23))
modes_54 = sum(E_percent(1:54))



