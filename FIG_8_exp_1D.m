clear

load data1 %measurement cov matrix, freq., wavenumber, microphone positions
%% Localization of sparse and coherent sources by orthogonal least squares
%  Gilles Chardon, François Ollivier, and José Picheral
%  The Journal of the Acoustical Society of America146, 4873 (2019); doi: 10.1121/1.5138931

% approx. 1min
%% Experiments, 1D, CMF-OLS and IMACS

% Generates figure 8 (1D experiment)

K = 5;
close all

% search space
L = 400;
xx = linspace(-1, 1, L)';

% distance between array and sources
Z = 5.18;

% dictionary, power = 1
D = dictionary1D(Pmic, [xx ones(L, 1)*Z], k);
D0 = D;

% norm = 1
inormsD = 1./sqrt(sum(abs(D.^2), 1));
D = D .* inormsD;

% beamforming dictionary
Dbf = D0./ (sum(abs(D0).^2, 1));

% diagonal removal
s = svd(Data);
Datad = Data - eye(length(Data))*min(s);

% beamforming
CRLCRL = Dbf'*Datad*Dbf;



tic
[sel, C, nres] = correl_ols(Datad, D, K);
Tols = toc

C = diag(inormsD(sel))*C*diag(inormsD(sel));


% more iterations to plot the decay
[sel2, C2, nres2] = correl_ols(Data, D, 12);



[W, A] = extract_sources(C, 2);
[W2, A2] = extract_sources(C2, 2);

% IMACS
Lh = 2;
c1 = clock;
normsD = sqrt(sum(abs(D.^2), 1));
P = imacs_cvx(Data, D, Lh, 8000, 40);
c2 = clock;

P = diag(inormsD)*P;

Tmacs = etime(c2, c1)

%%
close all


figure('Renderer', 'painters', 'Position', [10 10 600 200])
hold on
plot(xx, 10*log10(diag(abs(CRLCRL))), 'k')
stem([-0.60 -0.18 -0.07 0.62]+0.05, ones(4, 1)*70, '--k', 'markersize', eps)

stem(xx(sel), 10*log10(abs(W(:, 1).^2*A(1))), 'or', 'linewidth', 2, 'markersize', 15)
stem(xx(sel), 10*log10(abs(W(:, 2).^2*A(2))), 'xb', 'linewidth', 2, 'markersize', 15)

xlabel('Position (m)')
ylabel('CMF-OLS Power (dB)')
xlim([-1, 1])
ylim([50 70])
legend('Beamforming', 'Actual sources', 'CMF-OLS group 1', 'CMF-OLS group 2')


[idx1] = abs(P(:, 1)) > max(abs(P(:, 1)))/20;
[idx2] = abs(P(:, 2)) > max(abs(P(:, 2)))/20;

figure('Renderer', 'painters', 'Position', [10 10 600 200])

%plot(xx, abs(P.^2), 'linewidth', 2)
hold on
plot(xx, 10*log10(diag(abs(CRLCRL))), 'k')

stem([-0.60 -0.18 -0.07 0.62]+0.05, ones(4, 1)*70, '--k', 'markersize', eps)

stem(xx(idx1), 10*log10(abs(P(idx1, 1).^2)), 'or', 'linewidth', 2, 'markersize', 15);
stem(xx(idx2), 10*log10(abs(P(idx2, 2).^2)), 'xb', 'linewidth', 2, 'markersize', 15);
%plot(xx, abs(P.^2), 'linewidth', 2)
xlabel('Position (m)')
ylabel('MACS Power (dB)')
legend('Beamforming', 'Actual sources', 'MACS group 1', 'MACS group 2')

hold on
xlim([-1 1])
ylim([50 70])









figure
subplot(2, 1, 1)
plot(10*log10(abs(nres2)), '-xk')
xlabel('Iteration')
ylabel('Captured energy')
subplot(2, 1, 2)
stem(10*log10(svd(C)), 'k')
xlabel('Order')
ylabel('Singular values of Ĉ')
xlim([0 K+1])
ylim([20 70])

