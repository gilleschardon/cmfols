function D = dictionary(PX, PS, k)
% Green function
dx = PX(:, 1) - PS(:, 1)';
dz = PX(:, 3) - PS(:, 3)';
dy = PX(:, 2) - PS(:, 2)';


d = sqrt(dx.^2 + dz.^2 + dy.^2);

D = exp(- 1i * k * d)./d;

D = D;% ./ sqrt(sum(abs(D.^2), 1));

end