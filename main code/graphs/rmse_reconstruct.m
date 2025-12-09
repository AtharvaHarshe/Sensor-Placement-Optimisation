function RMSE = rmse_reconstruct(S, Phi, TrainSet, TestSet)

    % ensure TestSet is T × 54
    if size(TestSet,1)==54
        Test_full = TestSet';
    else
        Test_full = TestSet;
    end

    % mean-center
    TrainMean = mean(TrainSet,2);   % 54×1
    Test_center = Test_full - TrainMean';

    T = size(Test_full,1);
    n = size(Test_full,2);
    r = size(Phi,2);

    Phi_S = Phi(S,:);          % k×r
    PhiS_pinv = pinv(Phi_S);   % r×k

    Ypred = zeros(T,n);

    for t = 1:T
        yS  = Test_center(t,S)';   % k×1
        a   = PhiS_pinv * yS;      % r×1
        yrec = Phi*a + TrainMean;  % 54×1
        Ypred(t,:) = yrec';
    end

    RMSE = sqrt( mean( (Ypred - Test_full).^2 , 'all') );
end
