% =========================================================================
% Quadratic_Q9_Membrane_element.m
%
% 9-node isoparametric membrane finite element with embedded smeared
% reinforcement, used for material-nonlinear analysis of reinforced
% concrete under plane stress.
%
% Two operating modes (selected by 'Mensagem'):
%   'RigidPeq'   : returns the tangent element stiffness matrix and the
%                  equivalent nodal force vector.
%   'Resultado'  : recovers Gauss-point stresses, strains and updates the
%                  constitutive history (concrete + steel + cracking).
%
% sp matrix (8 rows x 9 Gauss points) returned in 'Resultado' mode:
%   sp(1,:) = sigma_1       -- maximum principal stress
%   sp(2,:) = sigma_2       -- minimum principal stress
%   sp(3,:) = tau_max       -- maximum in-plane shear stress
%   sp(4,:) = theta         -- principal direction angle [rad]
%   sp(5,:) = crack flag s1 -- 0=intact | 1=cracked | 2=saturated
%   sp(6,:) = crack flag s2 -- same for direction s2
%   sp(7,:) = w_s1          -- crack opening along s1 [m]
%   sp(8,:) = w_s2          -- crack opening along s2 [m]
% =========================================================================

function [var] = Quadratic_Q9_Membrane_element(Mensagem, Dados, d, coordsnoselem, cont, Numnoselem, ResConc, DadosConc, Resultado, ResAco)

switch Mensagem
    case 'RigidPeq'
        [k, peq, coordPtGauss] = Mrigidezelem(Dados, coordsnoselem, cont, Numnoselem, ResConc, DadosConc, ResAco);
        var.k            = k;
        var.peq          = peq;
        var.coordPtGauss = coordPtGauss;

    case 'Resultado'
        [s, e, ResConc, sp, DadosConc, Pintelemento, ResAco] = obtemresultados(Dados, d, cont, coordsnoselem, Numnoselem, ResConc, DadosConc, Resultado, ResAco);
        var.s            = s;
        var.e            = e;
        var.ResConc      = ResConc;
        var.sp           = sp;
        var.DadosConc    = DadosConc;
        var.Pintelemento = Pintelemento;
        var.ResAco       = ResAco;
end

return

% =========================================================================
% Mrigidezelem -- element tangent stiffness and equivalent nodal forces.
% (Renamed from the original 'Mrigidezbilinear': despite the legacy name,
%  this is the 9-node quadratic Q9 element, not a bilinear one.)
% =========================================================================
function [k, peq, coordPtGauss] = Mrigidezelem(Dados, coordsnoselem, cont, Numnoselem, ResConc, DadosConc, ResAco)

ordem        = Dados.ele(cont, 4);
[csi, eta, w] = PtsGauss2d(ordem);
NumdePontos  = ordem * ordem;

k   = zeros(Numnoselem * Dados.NumGDLNo, Numnoselem * Dados.NumGDLNo);
peq = zeros(Numnoselem * Dados.NumGDLNo, 1);
t   = Dados.ele(cont, 3);                  % concrete thickness

for c = 1:NumdePontos
    [N, NG]         = funcoesdeforma(csi, eta, c);
    [dNdcsi, dNdeta] = derfuncformacsieta(csi, eta, c);
    [Jac]           = matrizjacobiana(dNdcsi, dNdeta, coordsnoselem);
    [dNdx, dNdy]    = derfunformaxy(dNdcsi, dNdeta, Jac);
    [Dc, Da]        = MatrizconstitutivaOrtho(Dados, ResConc, DadosConc, cont, c, ResAco);
    [B]             = matrizB(dNdx, dNdy);
    J               = det(Jac);

    % Smeared steel contribution: equivalent thickness in x + in y
    EspAco = ResAco.TA(cont, 1) + ResAco.TA(cont, 3);

    k = k + t  * (B' * Dc * B * J * w(c)) ...
          + EspAco * (B' * Da * B * J * w(c));

    % Body forces (kept for completeness; usually zero in these examples)
    b      = [Dados.ele(cont, 5); Dados.ele(cont, 6)];
    peq    = peq + N' * b * t * J * w(c);

    % Cartesian coordinates of the Gauss point (for postprocessing)
    x = NG * coordsnoselem(:, 1);
    y = NG * coordsnoselem(:, 2);
    coordPtGauss(c, :) = [x, y];
end

return

% =========================================================================
% obtemresultados -- updates Gauss-point stresses, strains, principal
% directions, cracking state and internal force vector.
% =========================================================================
function [s, e, ResConc, sp, DadosConc, Pintelemento, ResAco] = obtemresultados(Dados, d, cont, coordsnoselem, Numnoselem, ResConc, DadosConc, Resultado, ResAco)

ordem        = Dados.ele(cont, 4);
[csi, eta, w] = PtsGauss2d(ordem);
NumdePontos  = ordem * ordem;

s = Resultado.el(cont).Sigma;
e = Resultado.el(cont).Epslon;

Pintelemento = zeros(size(Dados.conect, 2) * Dados.NumGDLNo, 1);
t            = Dados.ele(cont, 3);

% sp has 8 rows: 4 principal-stress quantities + 4 cracking-state values
NumPtsTotal = (Dados.ele(cont, 4))^2;
sp = zeros(8, NumPtsTotal);

% Carry over previous principal stresses and tangent moduli (history)
if Dados.IteracoesTotais > 0
    TpAnterior     = Resultado.el(cont).TensoesPrinc;
    DadosConc.SpA  = TpAnterior(1:4, :);
    DadosConc.EpA  = DadosConc.el(cont).E;
else
    DadosConc.SpA = zeros(4, 9);
    Ea            = Dados.ele(cont, 1);
    DadosConc.EpA = Ea * ones(2, 9);
end

S    = zeros(3, NumdePontos);
Saco = zeros(3, NumdePontos);

for c = 1:NumdePontos

    [N]             = funcoesdeforma(csi, eta, c); %#ok<ASGLU>
    [dNdcsi, dNdeta] = derfuncformacsieta(csi, eta, c);
    [Jac]           = matrizjacobiana(dNdcsi, dNdeta, coordsnoselem);
    [dNdx, dNdy]    = derfunformaxy(dNdcsi, dNdeta, Jac);
    [Dc, Da]        = MatrizconstitutivaOrtho(Dados, ResConc, DadosConc, cont, c, ResAco); %#ok<ASGLU>
    [B]             = matrizB(dNdx, dNdy);

    % Incremental stress / strain update at this Gauss point
    DeltaEpslon = B * d;
    e(:, c)     = e(:, c) + DeltaEpslon;

    DeltaSigma  = Dc * DeltaEpslon;
    s(:, c)     = s(:, c) + DeltaSigma;

    % Principal stresses and direction
    [TensoesPrinc] = calculatensoesprinc(s, c);
    sp(1:4, c)     = TensoesPrinc;

    DadosConc.el(cont).Sp(:, c) = TensoesPrinc;

    % Steel and concrete constitutive updates
    [ResAco]  = constitutive_steel(Dados, cont, c, e, ResAco);
    [ResConc] = constitutive_concrete(Dados, DadosConc, cont, c, ResConc);

    DadosConc.el(cont).E(:, c) = [ResConc.E(1); ResConc.E(2)];
    DadosConc.eiu              = ResConc.eiu;

    % --- Store cracking information (rows 5-8 of sp) ---------------------
    sp(5, c) = ResConc.Fissura(1);
    sp(6, c) = ResConc.Fissura(2);
    sp(7, c) = ResConc.AberturFissura(1);
    sp(8, c) = ResConc.AberturFissura(2);

    % --- Rotate principal stresses back to the global x-y system ---------
    Sprinc = ResConc.Sp(1:2);
    theta  = sp(4, c);
    Raio   = (Sprinc(1) - Sprinc(2)) / 2;
    sx     =  Raio*cos(2*theta) + Raio + Sprinc(2);
    sy     = -Raio*cos(2*theta) + Raio + Sprinc(2);
    txy    =  Raio*sin(2*theta);
    S(:, c) = [sx; sy; txy];

    % Smeared steel contribution to internal forces
    EspAco       = ResAco.TA(cont, 1) + ResAco.TA(cont, 3);
    Saco(:, c)   = [ResAco.ST.x(cont, c); ResAco.ST.y(cont, c); 0];

    J = det(Jac);
    Pintelemento = Pintelemento + (B' * S(:, c)    * J * w(c)) * t ...
                                + (B' * Saco(:, c) * J * w(c)) * EspAco;
end

return

% =========================================================================
% Shape functions (9-node Lagrangian quadrilateral)
% =========================================================================
function [N, NG] = funcoesdeforma(csi, eta, c)

x = csi(c);
y = eta(c);

N1 = x*y*(x-1)*(y-1)/4;
N2 = x*y*(x+1)*(y-1)/4;
N3 = x*y*(x+1)*(y+1)/4;
N4 = x*y*(x-1)*(y+1)/4;
N5 = y*(1-x^2)*(y-1)/2;
N6 = x*(x+1)*(1-y^2)/2;
N7 = y*(1-x^2)*(y+1)/2;
N8 = x*(x-1)*(1-y^2)/2;
N9 = (1-x^2)*(1-y^2);

N  = [N1 0  N2 0  N3 0  N4 0  N5 0  N6 0  N7 0  N8 0  N9 0
      0  N1 0  N2 0  N3 0  N4 0  N5 0  N6 0  N7 0  N8 0  N9];
NG = [N1 N2 N3 N4 N5 N6 N7 N8 N9];

return

% =========================================================================
% Shape function derivatives w.r.t. natural coordinates (csi, eta)
% =========================================================================
function [dNdcsi, dNdeta] = derfuncformacsieta(csi, eta, c)

x = csi(c);
y = eta(c);

dNdcsi = [1/4*(y*(y+1)*(2*x+1)),  1/4*(y*(y+1)*(2*x-1)),  1/4*(y*(y-1)*(2*x-1)),  1/4*(y*(y-1)*(2*x+1)), ...
          -(x*y*(y+1)),           -1/2*((y-1)*(y+1)*(2*x-1)), -(x*y*(y-1)),       -1/2*((y-1)*(y+1)*(2*x+1)), 2*x*((y-1)*(y+1))];

dNdeta = [1/4*(x*(2*y+1)*(x+1)),  1/4*(x*(2*y+1)*(x-1)),  1/4*(x*(2*y-1)*(x-1)),  1/4*(x*(2*y-1)*(x+1)), ...
          -1/2*((x-1)*(x+1)*(2*y+1)), -(x*y*(x-1)),       -1/2*((x-1)*(x+1)*(2*y-1)), -(x*y*(x+1)),       2*y*((x-1)*(x+1))];

return

% =========================================================================
% Jacobian matrix of the isoparametric mapping
% =========================================================================
function [Jac] = matrizjacobiana(dNdcsi, dNdeta, coordsnoselem)
j   = [dNdcsi; dNdeta];
Jac = j * coordsnoselem;
return

% =========================================================================
% Shape function derivatives w.r.t. cartesian coordinates (x, y).
% Uses backslash instead of inv() for numerical stability and speed.
% =========================================================================
function [dNdx, dNdy] = derfunformaxy(dNdcsi, dNdeta, Jac)
dNcsieta = [dNdcsi; dNdeta];
dNxy     = Jac \ dNcsieta;     % equivalent to inv(Jac)*dNcsieta, faster
dNdx     = dNxy(1, :);
dNdy     = dNxy(2, :);
return

% =========================================================================
% Orthotropic constitutive matrix in global coordinates (concrete) and
% smeared-steel constitutive matrix.
% Implements equations (5)-(11) of the paper:
%   - first iteration of the first step: isotropic linear elastic state
%   - subsequent iterations: orthotropic state rotated by current theta
% =========================================================================
function [Dc, Da] = MatrizconstitutivaOrtho(Dados, ResConc, DadosConc, cont, c, ResAco)

if Dados.IteracoesTotais == 0
    nu    = Dados.ele(cont, 2);
    E1    = Dados.ele(cont, 1);
    E2    = Dados.ele(cont, 1);
    E1a   = ResAco.TA(cont, 2);
    E2a   = ResAco.TA(cont, 4);
    theta = 0;
else
    nu    = ResConc.nu;
    E1    = abs(DadosConc.el(cont).E(1, c));
    E2    = abs(DadosConc.el(cont).E(2, c));
    theta = DadosConc.el(cont).Sp(4, c);
    E1a   = ResAco.ET.x(cont, c);
    E2a   = ResAco.ET.y(cont, c);
end

% Coordinate-rotation matrix T (Eq. 6 in the paper)
T = [  cos(theta)^2            sin(theta)^2            sin(theta)*cos(theta)
       sin(theta)^2            cos(theta)^2           -sin(theta)*cos(theta)
      -2*sin(theta)*cos(theta) 2*sin(theta)*cos(theta) cos(theta)^2 - sin(theta)^2];

% Orthotropic constitutive matrix in principal axes (Eq. 11)
C  = 1 / (1 - nu^2);
G  = (1/4) * ((E1 + E2) - 2*nu*sqrt(E1*E2));
Dp = C * [E1            nu*sqrt(E1*E2)  0
          nu*sqrt(E1*E2) E2             0
          0              0              G];

% Rotate to global axes (Eq. 7)
Dc = T' * Dp * T;

% Smeared steel: uniaxial in x and y, no shear coupling
Da = [E1a 0   0
      0   E2a 0
      0   0   0];

return

% =========================================================================
% Strain-displacement matrix B (size 3 x 18)
% =========================================================================
function [B] = matrizB(dNdx, dNdy)
B = zeros(3, 18);
for n = 1:9
    B(1, 2*n - 1) = dNdx(n);
    B(2, 2*n    ) = dNdy(n);
    B(3, 2*n - 1) = dNdy(n);
    B(3, 2*n    ) = dNdx(n);
end
return

% =========================================================================
% Principal stresses + direction from in-plane stress state
%   returns [s1; s2; tau_max; theta]
% =========================================================================
function [TensoesPrinc] = calculatensoesprinc(s, c)
sxx    = s(1, c);
syy    = s(2, c);
txy    = s(3, c);
centro = (sxx + syy) / 2;
a      = (sxx - syy) / 2;
R      = sqrt(a^2 + txy^2);
TensoesPrinc = [centro + R
                centro - R
                R
                atan2(txy, a) / 2];
return
