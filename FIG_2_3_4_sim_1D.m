clear variables
close all

%% Localization of sparse and coherent sources by orthogonal least squares
%  Gilles Chardon, François Ollivier, and José Picheral
%  The Journal of the Acoustical Society of America146, 4873 (2019); doi: 10.1121/1.5138931

% approx. 1min

% IMACS needs CVX
%% Simulations, 1D, CMF-OLS and IMACS

% Generates 
% - fig. 2, beamforming vs. CMF-OLS
% - fig. 3, CMF-OLS vs. IMACS
% - fig. 4, energy and singular values decay


%% General parameters
% Number of microphones
Nm = 19;

% Micorphone positions
xx = linspace(-1, 1, Nm)';
Xm = xx;
Ym = zeros(size(xx));

% Frequency and wavenumber
F = 3500;
k = 2 * pi * F / 344;


% search grid
L = 400;
xx = linspace(-1, 1, L);
Xs = xx;
Ys = 0*Xs;

% Distance between sources and microphone array
Z = 5;

% Dictionary of sources, power = 1
D = dictionary([Xm(:) Ym(:) zeros(Nm, 1)], [Xs(:) Ys(:) ones(L, 1)*Z], k);
D0 = D;
% Normalized dictionary
inormsD = 1./sqrt(sum(abs(D.^2), 1));
D = D .* inormsD;

% Beamforming dictionary (output = power)
Dbf = D0./ (sum(abs(D0).^2, 1));



%% Simulation of the data
% Number of snapshots to estimate the covariance matrix of the measurements
Nsnap = 500;
% amplitudes
sig_source = randn(2, Nsnap);

% Definition of the sources
% positions of the block of sources
Cs1 = [0.9 0.0 ; 0 0 ; -0.8 0.0];
% amplitudes
A1 = [0.2 ; -1 ; 0.4]*100;
% steering vectors times amplitudes
S1 = dictionary([Xm(:) Ym(:) zeros(Nm, 1)], [Cs1 ones(3, 1)*Z], k) * A1;

% second block
Cs2 = [-0.3 0.0 ; 0.5 0.0];
A2 = [0.3 ; 0.8*i]*100;
S2 = dictionary([Xm(:) Ym(:) zeros(Nm, 1)], [Cs2 ones(2, 1)*Z], k) * A2;

sources = [S1 S2];

% measure signals
data = sources * sig_source;

% noise
pwr = norm(data, 'fro')^2;
SNR = 0;
noise = randn(size(data));
noise = noise * sqrt(1/ norm(noise, 'fro')^2 * pwr * 10^(-SNR/10));
data2 = data + noise;

%% Processing
% Number of iterations
K = 5;

% Covariance matrix of the measurements, partial removal of the diagonal
Data = data2*data2'/Nsnap;
s = svd(Data);
Datad = Data - eye(length(Data))*min(s);

% beamforming
CRLCRL = Dbf'*Datad*Dbf;

% CMF-OLS
tic
[sel, Cols, nres] = correl_ols(Datad, D, K);
Tols = toc



% run with more iterations to plot the energy decay
[sel2, C2, nres2] = correl_ols(Datad, D, 14);

% Scaling
Cols = diag(inormsD(sel))*Cols*diag(inormsD(sel));

% Complete matrix, with all grid points
Colstot = zeros(L, L);
Colstot(sel, sel) = Cols;

% Extraction of the sources
[W, A] = extract_sources(Cols, 2);

% Estimated positions
Cest = [Xs(sel)' Ys(sel)'];

% MACS
c1 = clock;
P = imacs_cvx(Datad, D, 2, 230, 50); 
c2 = clock;
P = diag(inormsD)*P;

Tmacs = etime(c2, c1)



%% Plots
% Beamforming vs. CMF-OLS
figure
subplot(1, 2, 1)
img = (10*log10(abs(CRLCRL)));
img(img < max(img(:)) - 20) = max(img(:))-20;
imagesc(xx,xx,real(img))
hold on
scatter(Cs1(:, 1), Cs1(:, 1), 100, 'xk', 'linewidth', 2)
scatter(Cs2(:, 1), Cs2(:, 1), 100, 'ok', 'linewidth', 2)

axis ij
axis square
colormap(hot)
xlabel('x')
ylabel('x')
title('DAS-C')

subplot(1, 2, 2)

hold on
for u = 1:K
    for v = 1:K
        scatter(xx(sel(u)), xx(sel(v)), (abs(Cols(u, v))/500)^2, 'k', 'filled')
    end
end
xlabel('x')
ylabel('x')
axis ij
axis image
title('CMF-OLS')
xlim([-1, 1]);
ylim([-1, 1]);
%% CMF-OLS


figure('Renderer', 'painters', 'Position', [10 10 600 200])
hold on

plot(xx, 10*log10(diag(abs(CRLCRL))), 'k')
stem([Cs1(:, 1); Cs2(:, 1)] , 10*log10([abs(A1).^2 ; abs(A2).^2]), '--+k', 'linewidth', 1, 'markersize', 20)


stem(xx(sel), 10*log10(abs(W(:, 2).^2*A(2))), 'xb', 'linewidth', 2, 'markersize', 15)

stem(xx(sel), 10*log10(abs(W(:, 1).^2*A(1))), 'or', 'linewidth', 2, 'markersize', 15)

legend('Beamforming', 'Sources', 'CMF-OLS group 1', 'CMF-OLS group 2')

xlabel('Position (m)')
ylabel('CMF-OLS Power (dB)')

xlim([-1 1])
ylim([20 45])

%% Energy and singular values decay

figure('Renderer', 'painters', 'Position', [10 10 600 200])


subplot(1, 2, 1)
plot(10*log10(abs(nres2)), '-xk')
xlabel('Iteration')
ylabel(sprintf('Residual\nenergy (dB)'))
subplot(1, 2, 2)
stem(10*log10(svd(Cols)), 'k')
xlabel('Order')
ylabel(sprintf('Singular values\nof Ĉ (dB)'))
xlim([0 K+1])
ylim([20 50])


[idx1] = abs(P(:, 1)) > max(abs(P(:, 1)))/20;
[idx2] = abs(P(:, 2)) > max(abs(P(:, 2)))/20;

%% IMACS
figure('Renderer', 'painters', 'Position', [10 10 600 200])
hold on

plot(xx, 10*log10(diag(abs(CRLCRL))), 'k')

stem([Cs1(:, 1); Cs2(:, 1)] , 10*log10([abs(A1).^2 ; abs(A2).^2]), '--+k', 'linewidth', 1, 'markersize', 20)
stem(xx(idx1), 10*log10(abs(P(idx1, 1).^2)), 'xb', 'linewidth', 2, 'markersize', 15);
stem(xx(idx2), 10*log10(abs(P(idx2, 2).^2)), 'or', 'linewidth', 2, 'markersize', 15);

xlabel('Position (m)')
ylabel('MACS Power (dB)')

hold on
xlim([-1 1])
ylim([20 45])

legend('Beamforming', 'Sources', 'IMACS group 1', 'IMACS group 2');


