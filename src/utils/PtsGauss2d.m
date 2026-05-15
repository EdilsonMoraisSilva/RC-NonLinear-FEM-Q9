% =========================================================================
% PtsGauss2d
%
% Gauss-Legendre quadrature points and weights in two dimensions on the
% reference square [-1, 1] x [-1, 1], obtained as the tensor product of
% the 1D rule of given order.
%
% Universidade Federal do Pará (UFPA)
% Original author: Remo Magalhães de Souza  (remo@ufpa.br)
% =========================================================================

function [csi, eta, w] = PtsGauss2d(ordem)

[csi1, w1] = PtsGauss1d(ordem);
k = 0;

for i = 1:ordem
    for j = 1:ordem
        k = k + 1;
        csi(k) = csi1(j);
        eta(k) = csi1(i);
        w(k)   = w1(i) * w1(j);
    end
end

return
