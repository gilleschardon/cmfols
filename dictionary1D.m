function D = dictionary1D(PX, PS, k)
% Green function (2D case)
dx = PX(:, 1) - PS(:, 1)';
dz = PX(:, 2) - PS(:, 2)';


d = sqrt(dx.^2 + dz.^2);

D = exp(- 1i * k * d)./d;

D = D;

end