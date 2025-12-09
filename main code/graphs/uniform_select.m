function S = uniform_select(coords, k)
    n = size(coords,1);
    S = zeros(1,k);

    S(1) = randi(n);

    for i = 2:k
        best_d = -Inf; best_s = -1;
        for s = 1:n
            if ismember(s,S), continue; end
            d = min(vecnorm(coords(S(1:i-1),:) - coords(s,:),2,2));
            if d > best_d
                best_d = d;
                best_s = s;
            end
        end
        S(i) = best_s;
    end
end
