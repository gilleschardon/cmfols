%measurement cov matrix, freq., wavenumber, microphone positions
%% Localization of sparse and coherent sources by orthogonal least squares
%  Gilles Chardon, François Ollivier, and José Picheral
%  The Journal of the Acoustical Society of America146, 4873 (2019); doi: 10.1121/1.5138931

% approx. 1min
%% Experiments, 2D, CMF-OLS, uncorrelated sources

% Generates figures 9 and 10 (partial)

% approx 10sec.
%%

load data_unc
close all
Lx = 400;
Ly = 200;
xx = linspace(-2, 2, Lx)';
yy = linspace(-1, 1, Ly)';

[Xg Yg] = meshgrid(xx, yy);

Z = 4.3;

D = dictionary(Pmic, [Xg(:) Yg(:) ones(Lx*Ly, 1)*Z], k);
D0 = D;

inormsD = 1./sqrt(sum(abs(D.^2), 1));
D = D .* inormsD;
Dbf = D0./ (sum(abs(D0).^2, 1));


Cbf = sum(conj(Dbf) .* (Data*Dbf), 1);


%%
% number of iterations
K = 4;

% diagonal removal
s = svd(Data);
Datad = Data - eye(length(Data))*min(s);

% CMF-OLS
tic
[sel, Cols, nres] = correl_ols(Datad, D, K);

Cols = diag(inormsD(sel))*Cols*diag(inormsD(sel));
Tols = toc

[sel2, Cols2, nres2] = correl_ols(Data, D, 12);



[W, A] = extract_sources(Cols, 4);
[W2, A2] = extract_sources(Cols2, 4);


%%
close all


 figure

 img = (10*log10(abs(reshape(Cbf, Ly, Lx))));
 img(img < max(img(:)) - 20) = max(img(:))-20;
 imagesc(xx,yy,real(img))

 xlabel('x')
 ylabel('y')

 axis xy
 axis square
 colormap(hot)
colorbar
hold on


scatter(Xg(sel), Yg(sel), ((abs(W(:, 1, 1).^2*A(1, 1)))/12000+eps).^2, '+k', 'linewidth', 2)
scatter(Xg(sel), Yg(sel), ((abs(W(:, 2, 1).^2*A(2, 1)))/12000+eps).^2, 'xk', 'linewidth', 2)
scatter(Xg(sel), Yg(sel), ((abs(W(:, 3, 1).^2*A(3, 1)))/12000+eps).^2, 'ok', 'linewidth', 2)
scatter(Xg(sel), Yg(sel), ((abs(W(:, 4, 1).^2*A(4, 1)))/12000+eps).^2, 'sk', 'linewidth', 2)

xlim([-2, 2])
ylim([-1.5, 0])
legend('CMF-OLS group 1', 'CMF-OLS group 2', 'CMF-OLS group 3', 'CMF-OLS group 4')

xlabel('x')
ylabel('y')
axis xy
axis image



figure
subplot(2, 1, 1)
plot(10*log10(abs(nres2)), '-xk')
xlabel('Iteration')
ylabel('Captured energy')
subplot(2, 1, 2)
stem(10*log10(svd(Cols)), 'k')
xlim([0 K])
xlabel('Order')
ylabel('Singular value of Ĉ')

%%
% Uncomment to generate fig 10, with the results of FIG_manip_corr
% 
% figure
% subplot(1, 2, 1)
% 
% plot(10*log10(abs(nres2)), '-ok')
% hold on
% plot(10*log10(abs(nres1)), '-xk')
% xlabel('Iteration')
% ylabel(sprintf('Residual\nenergy (dB)'))
% ylim([63 69])
% legend('Uncorrelated sources');%, 'Correlated sources')
% 
% subplot(1, 2, 2)
% stem(1:4, 10*log10(svd(Cols)), 'ok', 'MarkerSize', 10)
% hold on
% stem(1:4, 10*log10(ss1), 'xk', 'MarkerSize', 10)
% legend('Uncorrelated sources')%, 'Correlated sources')
% 
% xlim([0.5 K+0.5])
% xlabel('Order')
% ylabel('Singular value of Ĉ')
% %ylim([20 50])
% ylim([20 70])
