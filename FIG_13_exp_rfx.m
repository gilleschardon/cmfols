clear variables



load data_rfx
%measurement cov matrix, freq., wavenumber, microphone positions
%% Localization of sparse and coherent sources by orthogonal least squares
%  Gilles Chardon, François Ollivier, and José Picheral
%  The Journal of the Acoustical Society of America146, 4873 (2019); doi: 10.1121/1.5138931

% approx. 10sec.
%% Experiments, 2D, CMF-OLS, reflections

% Generates figure 13

% approx 10sec.

close all
Lx = 600;
Ly = 150;
xx = linspace(-1, 5, Lx)';
yy = linspace(-1.5, 0, Ly)';

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

%CRLCRL = Dbf'*Datad*Dbf;



tic
[sel, Cols, nres] = correl_ols(Datad, D, K);

Cols = diag(inormsD(sel))*Cols*diag(inormsD(sel));
Tols = toc

[sel2, Cols2, nres2] = correl_ols(Data, D, 12);



[W, A] = extract_sources(Cols, 2);
[W2, A2] = extract_sources(Cols2, 2);


%save manip1
%%
%load manip1
close all


 figure
% 
% 
% 
% 
 img = (10*log10(abs(reshape(Cbf, Ly, Lx))));
 img(img < max(img(:)) - 20) = max(img(:))-20;
 imagesc(xx,yy,real(img))
% hold on
% scatter(Cs1(:, 1), Cs1(:, 2), 300, 'xw', 'linewidth', 5)
% scatter(Cs2(:, 1), Cs2(:, 2), 300, 'ow', 'linewidth', 5)
% 
 xlabel('x')
 ylabel('y')
% 
 axis xy
 axis square
 colormap(hot)
colorbar
%figure
hold on
%scatter([Cs1(:, 1) ; Cs2(:, 1)], [Cs1(:, 2) ; Cs2(:, 2)] ,200, '+k')


scatter(Xg(sel), Yg(sel), (abs(W(:, 1, 1).^2*A(1, 1))).^2/700000000+eps, '+g', 'linewidth', 2)
scatter(Xg(sel), Yg(sel), (abs(W(:, 2, 1).^2*A(2, 1))).^2/700000000+eps, 'xg', 'linewidth', 2)

%scatter(Xg(selomp), Yg(selomp), 20, 'ok', 'linewidth', 2)


xlim([-1, 5])
ylim([-1, 1])
legend('CMF-OLS group 1', 'CMF-OLS group 2')

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
%ylim([20 50])

% 
% figure
% hold on
% imglog = log(abs(reshape(sum(P.^2, 2), L, L)));
% imglog(imglog < max(imglog(:)) - 10) = max(imglog(:)) - 10;
% imagesc(xx, xx, imglog)
% scatter(Cs1(:, 1), Cs1(:, 2), 300, 'ok')
% scatter(Cs2(:, 1), Cs2(:, 2), 300, 'sk')
% axis square
% axis xy
% colormap(1-gray)
% xlim([-1, 1])
% ylim([-1, 1])
% 
% xlabel('x')
% ylabel('y')
