function S = pv_select(coords, k)
    n = size(coords,1);
    S = [];

    % start with random seed
    S(1) = randi(n);

    for i = 2:k
        varScore = zeros(n,1);
        for s = 1:n
            if ismember(s,S), varScore(s) = -Inf; continue; end
            d = min(vecnorm(coords(S,:) - coords(s,:),2,2));
            varScore(s) = d;  % crude proxy of predictive variance
        end
        [~,next] = max(varScore);
        S(end+1) = next;
    end
end
