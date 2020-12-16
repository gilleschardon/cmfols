function P = imacs_cvx(G, A, Lh, xi, Niter)

% IMACS
%T. Yardibi, J. Li, P. Stoica, N. S. Zawodny, and L. N. Cattafesta,
%A covariance fitting approach for correlated acoustic source mapping,
%J. Acoust. Soc. Am.127(5), 2920–2931 (2010)


%Y. Li, M. Li, D. Yang, and C. Gao, 
%Research of the improved mapping ofacoustic correlated sources method,
%Appl. Acoust.145, 290–304 (2019)

% requires CVX


[N, L] = size(A);

[U, Lam] = eig(G);
Lam = diag(Lam);
[~,idx] = sort(abs(Lam), 'descend');
beta = sum(Lam);

if ~exist('xi', 'var')
   xi = sqrt(beta * L * Lh);
end

Q = eye(Lh);
Gt = U(:, idx(1:Lh)) * diag(sqrt(Lam(idx(1:Lh))));


 
% AAA = [];
% for v = 1:Lh
%     AAA = blkdiag(AAA, A);
% end
PP = zeros(L, L);
for u = 1:Niter
    gtq = Gt*Q';
 %   C = cvx_lasso(A,Gt*Q, xi);
cvx_begin
variable C(L, Lh)

minimize norm(gtq - A*C, 'fro')
subject to
norm(C(:), 1) <= xi
cvx_end

    %C = spg_lasso(AAA,gtq(:), xi);
    %C = reshape(C, L, Lh);
    [uu, ss, vv] = svd(Gt'*A*C);
    Q = vv*uu';
    %xi = norm(C(:), 1);
    PPp = C*C';
    
    norm(PP, 'fro')^2;
    norm(PPp-PP, 'fro')^2/norm(PP, 'fro')^2
    PP = PPp;
    
end
    
P = C;
    
    
end

