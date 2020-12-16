function [W, A] = extract_sources(C, K)

% Extraction of the sources from the covariance matrix
% see section V

[V, L] = eig((C+C')/2);
L = diag(L);
[L, idx] = sort(L, 'descend');

PWR = zeros(K, 1);
V= V(:, idx);
V = V(:, 1:K);
L = L(1:K);
W = zeros(size(V));

VV = V;
for v = 1:K
    K = sum(abs(VV).^2, 2);
    [~, idx] = max(K);
    coeffs = VV(idx, :);
    coeffs = coeffs/norm(coeffs);
    W(:, v) = VV * coeffs';
    VV = VV - VV * (coeffs' * coeffs);
end
    
WW = zeros(size(W, 1)^2, size(W, 2));
for w = 1:size(W, 2)
    ww = W(:, w)*W(:, w)';
    WW(:, w) = ww(:);
end

CC = C(:);

A = WW\CC;



end