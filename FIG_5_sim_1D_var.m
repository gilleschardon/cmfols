clear variables
close all

%% Localization of sparse and coherent sources by orthogonal least squares
%  Gilles Chardon, François Ollivier, and José Picheral
%  The Journal of the Acoustical Society of America146, 4873 (2019); doi: 10.1121/1.5138931

%% CMF-OLS, 1D, CMF-OLS and IMACS, simulations, performances

% Generates fig. 4, performances of DAS-C, IMACS, CMF-OLS

% approx 20min.
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
Nsnap = 100;

% Coordinates
Cs2 = [-0.3 0.0 ; 0.5 0.0];
% Ampltiudes
A2 = [0.5; 1*i]*100;
% steering vectors time amplitudes
S2 = dictionary([Xm(:) Ym(:) zeros(Nm, 1)], [Cs2 ones(2, 1)*Z], k) * A2;

sources = S2;


% Estimated positions and covariances
NT =  20; % number of draws (50 in the paper)

NSNR = 9; % number of SNRs

POSBF = zeros(2, NSNR, NT);
POSOLS = zeros(2, NSNR, NT);
POSIMACS = zeros(2, NSNR, NT);

CORBF = zeros(2, 2, NSNR, NT);
COROLS = zeros(2, 2, NSNR, NT);
CORIMACS = zeros(2, 2, NSNR, NT);

SNRs = linspace(-20, 20, NSNR);


% loop over the SND
for nsnr = 1:NSNR
    waitbar(nsnr/NSNR)
    
% loop over draws
for nt = 1:NT
	sig_source = randn(1, Nsnap);
	data = sources * sig_source;

	% noise
	pwr = norm(data, 'fro')^2;
	SNR = SNRs(nsnr);%-20;
	noise = randn(size(data));
	noise = noise * sqrt(1 / norm(noise, 'fro')^2 * pwr * 10^(-SNR/10));
	data2 = data + noise;


% CMF-OLS iterations
K = 2;

% Covariance matrix, diagonal removal
Data = data2*data2'/Nsnap;
s = svd(Data);
Datad = Data - eye(length(Data))*min(s);

% CMS-OLS
[sel, Cols, nres] = correl_ols(Datad, D, K);

Cols = diag(inormsD(sel))*Cols*diag(inormsD(sel));
Colstot = zeros(L, L);
Colstot(sel, sel) = Cols;

% positions and covariances
xsel = (xx(sel));
[~, ix] = sort(xsel);
POSOLS(:, nsnr, nt) = xsel(ix);
COROLS(:, :, nsnr, nt) = Cols(ix, ix);

Cest = [Xs(sel)' Ys(sel)'];

% MACS
P = imacs_cvx(Datad, D, 1, 100, 20); 
P = diag(inormsD)*P;

 bfimacs = abs(diag(P*P'));
 bbb = P*P';
% positions and covariances
% we take the maxes for positive and negative positions (we add here a bit
% of informations, to all the methods)
 [~, idxm] = max(bfimacs .*(xx<0)');
 [~, idxp] = max(bfimacs.*(xx>=0)');
 POSIMACS(:, nsnr, nt) = xx([idxm idxp]);

 sel = [idxm idxp];
 CORIMACS(:, :, nsnr, nt) = bbb(sel, sel);




% beamforming
bf = sum((Dbf'*Datad) .* Dbf.', 2);

xxm = xx(xx<0);
xxp = xx(xx>=0);

[~, idxm] = max(bf .*(xx<0)');
[~, idxp] = max(bf .*(xx>=0)');

POSBF(:, nsnr,  nt) = xx([idxm idxp]);

CORBF(:, :, nsnr, nt) = Dbf(:, [idxm idxp])'*Datad * Dbf(:, [idxm idxp]);

end % draws
end % SNR

save simu1


%%
load simu1


% mean errors
mbf = mean((POSBF - repmat([-0.3 0.5]', 1, NSNR, NT)).^2, 3);
mols = mean((POSOLS - repmat([-0.3 0.5]', 1, NSNR, NT)).^2, 3);
mmacs = mean((POSIMACS - repmat([-0.3 0.5]', 1, NSNR, NT)).^2, 3);

% actual covariance matrix
C = A2*A2';

cpols = zeros(NSNR, 1);
cpbf = zeros(NSNR, 1);

cpimacs = zeros(NSNR, 1);
zzzols = zeros(NSNR, NT);
zzzbf = zeros(NSNR, NT);

% honestly I don't know
for vv = 1:NSNR
    for ww = 1:NT
        zzzols(vv, ww) = norm(reshape(COROLS(:, :, vv, ww) - C, 4, 1), 'inf');
        zzzbf(vv, ww) = norm(reshape(CORBF(:, :, vv, ww) - C, 4, 1), 'inf');
        
        relerr = (COROLS(:, :, vv, ww) - C)./C;
        cpols(vv) = cpols(vv) + sum(abs(relerr(:)).^2);
        
        relerr = (CORBF(:, :, vv, ww) - C)./C;
        cpbf(vv) = cpbf(vv) + sum(abs(relerr(:)).^2);
        relerr = (CORIMACS(:, :, vv, ww) - C)./C;
        cpimacs(vv) = cpimacs(vv) + sum(abs(relerr(:)).^2);
        

    end
end

% errors on the covariances (/4 because 4 terms)
cpols = cpols/NT/4;
cpbf = cpbf/NT/4;
cpimacs = cpimacs/NT/4;


figure

subplot(1, 2, 1)
semilogy(SNRs, mols(1, :)', 'k-', 'linewidth', 2)
xlabel('SNR')
ylabel({'MSE of the','position estimation'})
hold on
plot(SNRs, mbf(1, :)', 'k--', 'linewidth', 2)
plot(SNRs, mmacs(1, :)', 'k-.', 'linewidth', 2)
semilogy(SNRs, mols(2, :)', 'k-', 'linewidth', 2)
plot(SNRs, mbf(2, :)', 'k--', 'linewidth', 2)
plot(SNRs, mmacs(2, :)', 'k-.', 'linewidth', 2)
subplot(1, 2, 2)


semilogy(SNRs, cpols(:, 1), 'k-', 'linewidth', 2)
hold on
semilogy(SNRs, cpbf(:, 1), 'k--', 'linewidth', 2)
xlabel('SNR')
ylabel({'MSE of the power and','covariances estimation'})
semilogy(SNRs, cpimacs(:, 1), 'k-.', 'linewidth', 2)



legend('CMF-OLS', 'DAS-C', 'IMACS')

