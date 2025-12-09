function [knee_idx, knee_x, knee_y, MC] = findKneePoint(x, F)
% findKneePoint  Find knee point of a Pareto curve using Menger curvature.
%
% Implements Appendix A of:
%   "Multi-objective optimization for sensor placement: An integrated
%    combinatorial approach with reduced order model and Gaussian process"
%
% INPUTS:
%   x  - vector of x-coordinates (e.g., number of sensors)
%   F  - vector of y-coordinates (e.g., weighted MI at each x)
%
% OUTPUTS:
%   knee_idx - index of the knee point (2..n-1)
%   knee_x   - x(knee_idx)
%   knee_y   - F(knee_idx)
%   MC       - Menger curvature values for each point (NaN at 1 and n)
%
% The Menger curvature at point i (2 <= i <= n-1) is:
%   MC(x_i) = sqrt(A - B^2) / (|pq| * |qr| * |rp|)
% where:
%   p = (x_{i-1}, F_{i-1}), q = (x_i, F_i), r = (x_{i+1}, F_{i+1})
%   |pq| = distance between p and q, etc.
%   A = 4 |pq|^2 |qr|^2
%   B = |pq|^2 + |qr|^2 - |rp|^2

    % --- basic checks ---
    x = x(:);
    F = F(:);

    n = numel(x);
    if n < 3
        error('findKneePoint:NotEnoughPoints', ...
              'Need at least 3 points to define a knee.');
    end
    if numel(F) ~= n
        error('findKneePoint:SizeMismatch', ...
              'x and F must have the same length.');
    end

    % Preallocate curvature array
    MC = zeros(n, 1);

    % Loop over interior points (2..n-1)
    for i = 2:n-1
        % Points p, q, r
        xp = x(i-1); Fp = F(i-1);
        xq = x(i);   Fq = F(i);
        xr = x(i+1); Fr = F(i+1);

        % Distances |pq|, |qr|, |rp|
        pq = sqrt( (xq - xp)^2 + (Fq - Fp)^2 );
        qr = sqrt( (xr - xq)^2 + (Fr - Fq)^2 );
        rp = sqrt( (xr - xp)^2 + (Fr - Fp)^2 );

        % If any segment is (numerically) zero length, skip
        if pq == 0 || qr == 0 || rp == 0
            MC(i) = NaN;
            continue;
        end

        % A = 4|pq|^2|qr|^2, B = |pq|^2 + |qr|^2 - |rp|^2
        pq2 = pq^2;
        qr2 = qr^2;
        rp2 = rp^2;

        A = 4 * pq2 * qr2;
        B = pq2 + qr2 - rp2;

        % Numerical safety: clamp A - B^2 to >= 0
        disc = A - B^2;
        if disc < 0
            disc = 0;
        end

        MC(i) = sqrt(disc) / (pq * qr * rp);
    end

    % Knee point is index with maximum curvature among valid entries
    [~, knee_idx] = max(MC);

    knee_x = x(knee_idx);
    knee_y = F(knee_idx);
end

