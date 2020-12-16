clear variables



load data_corr

%measurement cov matrix, freq., wavenumber, microphone positions
%% Localization of sparse and coherent sources by orthogonal least squares
%  Gilles Chardon, François Ollivier, and José Picheral
%  The Journal of the Acoustical Society of America146, 4873 (2019); doi: 10.1121/1.5138931

% approx. 10sec.
%% Experiments, 2D, CMF-OLS, correlated sources

% Generates figures 10 (partial), 11 and 12

% approx 10sec.


%%
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
K = 4;
s = svd(Data);
Datad = Data - eye(length(Data))*min(s);


tic
[sel, Cols, nres] = correl_ols(Datad, D, K);

Cols = diag(inormsD(sel))*Cols*diag(inormsD(sel));
Tols = toc

[sel2, Cols2, nres2] = correl_ols(Data, D, 12);



[W, A] = extract_sources(Cols, 2);
[W2, A2] = extract_sources(Cols2, 2);


%%
close all


 figure

 img = (10*log10(abs(reshape(Cbf, Ly, Lx))));
 img(img < max(img(:)) - 20) = max(img(:))-20;
 imagesc(xx,yy,real(img))
 colormap(hot)
 
colorbar
hold on
scatter(Xg(sel), Yg(sel), ((abs(W(:, 1, 1).^2*A(1, 1)))/12000+eps).^2, '+k', 'linewidth', 2)
scatter(Xg(sel), Yg(sel), ((abs(W(:, 2, 1).^2*A(2, 1)))/12000+eps).^2, 'xk', 'linewidth', 2)
xlim([-2, 2])
ylim([-1.5, 0])
legend('CMF-OLS group 1', 'CMF-OLS group 2')

xlabel('x')
ylabel('y')
axis xy
axis image

 figure
 

 subplot(2, 1, 2)
 
 CRLCRL = Dbf'*Datad*Dbf(:, sel(4));

 img = (10*log10(abs(reshape(CRLCRL, Ly, Lx))));
 img(img < max(img(:)) - 20) = max(img(:))-20;
 imagesc(xx,yy,real(img))
 colormap(hot)
 
colorbar
hold on
scatter(Xg(sel), Yg(sel), (abs(Cols(:, 4))/10000+eps).^2, 'xk', 'linewidth', 2)
xlim([-2, 2])
ylim([-1.5, 0])
legend('CMF-OLS covariances with source 4')

xlabel('x')
ylabel('y')
axis xy
axis image

 
 subplot(2, 1,1)
 
 CRLCRL = Dbf'*Datad*Dbf(:, sel(1));

 img = (10*log10(abs(reshape(CRLCRL, Ly, Lx))));
 img(img < max(img(:)) - 20) = max(img(:))-20;
 imagesc(xx,yy,real(img))
 colormap(hot)
 
colorbar
hold on
scatter(Xg(sel), Yg(sel), (abs(Cols(:, 1))/10000+eps).^2, '+k', 'linewidth', 2)
xlim([-2, 2])
ylim([-1.5, 0])
legend('CMF-OLS covariances with source 1')

xlabel('x')
ylabel('y')
axis xy
axis image






figure
subplot(2, 1, 1)
plot(10*log10(abs(nres2)), '-xk')
xlabel('Iteration')
ylabel('Residual energy (d)')
subplot(2, 1, 2)
stem(10*log10(svd(Cols)), 'k')
xlim([0 K])
xlabel('Order')
ylabel('Singular value of Ĉ')
