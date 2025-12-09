

X = coords(:,2:3);          % only x and y coordinates
r = size(Phi, 2);           % number of POD modes
% Variable Inisitalisation
Sigma_all = cell(r,1);      % store Σ_i
mu_all    = cell(r,1);      % store μ_i
gprModels = cell(r,1);      % store GP models

for i = 1:r
    
    % ----- POD mode i becomes GP training target -----
    y_i = Phi(:,i);   % 54x1 vector
    
    % ----- Train GP model for mode i -----
    gprModels{i} = fitrgp(X, y_i, ...
        'KernelFunction','squaredexponential', ...
        'BasisFunction','none', ...
        'FitMethod','exact', ...
        'Standardize',true);
    
    % ----- Standardize coordinates exactly how GP did -----
    X_std = (X - gprModels{i}.Impl.StdMu') ./ gprModels{i}.Impl.StdSigma';
    
    % ----- Build kernel from learned hyperparameters -----
    theta = gprModels{i}.Impl.ThetaHat;          % [ell; sigma_f]
    kernelFcn = gprModels{i}.Impl.Kernel.makeKernelAsFunctionOfXNXM(theta);
    
    % ----- Compute covariance matrix Σ_i = K_i(X,X) + σ_n^2 I -----
    K_i = kernelFcn(X_std, X_std);
    Sigma_all{i} = K_i + gprModels{i}.Impl.SigmaHat^2 * eye(size(K_i));
    
    % ----- Mean (optional, not needed for MI) -----
    mu_all{i} = predict(gprModels{i}, X);
    
end

clear i k_i kernelFcn theta X X_std y_i K_i 