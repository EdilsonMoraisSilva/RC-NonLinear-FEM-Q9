% =========================================================================
% reinforcement_setup
%
% Computes the smeared reinforcement equivalent thickness per element,
% used by the distributed-reinforcement representation described in
% Section 2.6 of the paper (Figure 5).
%
% For each reinforcement bar group i = 1 (x direction) and i = 2 (y
% direction), the routine loops over all bars defined in Dados.Arm.Barras{i}
% and over all finite elements. When an element falls inside the rectangular
% region defined by the bar group coordinates, its equivalent steel
% thickness is computed as:
%
%       espessuraequivalente(elem, i) = V_steel(bar) / A_element
%
% where V_steel = n_bars * (pi * D^2 / 4) * L  (L is the element edge
% length along the bar direction), and A_element is the element area.
%
% Output:
%   TaxArmElem  - NumElem x 4 matrix:
%       column 1 : equivalent steel thickness in x  [m]
%       column 2 : initial steel modulus in x       [kPa]
%       column 3 : equivalent steel thickness in y  [m]
%       column 4 : initial steel modulus in y       [kPa]
% =========================================================================

function [TaxArmElem] = reinforcement_setup(Dados)

N = size(Dados.conect, 1);
espessuraequivalente{1} = zeros(N, 2);
espessuraequivalente{2} = zeros(N, 2);

for i = 1:2

    for j = 1:size(Dados.Arm.Barras{i}, 1)

        coordAreaAco{i}(j, :) = Dados.Arm.Barras{i}(j, 4:7);

        for k = 1:N
            Nosdoelem      = Dados.conect(k, 1:4);
            coordnosdoelem = Dados.coords(Nosdoelem, :);

            % Bar volume inside this element
            n  = Dados.Arm.Barras{i}(j, 1);
            D  = Dados.Arm.Barras{i}(j, 2);
            Lx = coordnosdoelem(1, 1) - coordnosdoelem(2, 1);
            Ly = coordnosdoelem(2, 2) - coordnosdoelem(3, 2);
            L  = [Lx; Ly];

            VolArm{i}(j) = n * (pi * D^2 / 4) * L(i);

            % Element / bar bounding box
            MaxCoordNoelem.x  = max(coordnosdoelem(:, 1));
            MaxCoordNoelem.y  = max(coordnosdoelem(:, 2));
            MaxCoordAreaaco.x = max(coordAreaAco{i}(j, 1:2));
            MaxCoordAreaaco.y = max(coordAreaAco{i}(j, 3:4));

            MinCoordNoelem.x  = min(coordnosdoelem(:, 1));
            MinCoordNoelem.y  = min(coordnosdoelem(:, 2));
            MinCoordAreaaco.x = min(coordAreaAco{i}(j, 1:2));
            MinCoordAreaaco.y = min(coordAreaAco{i}(j, 3:4));

            EX = MaxCoordNoelem.x;   AX = MaxCoordAreaaco.x;
            EY = MaxCoordNoelem.y;   AY = MaxCoordAreaaco.y;
            Ex = MinCoordNoelem.x;   Ax = MinCoordAreaaco.x;
            Ey = MinCoordNoelem.y;   Ay = MinCoordAreaaco.y;

            % Element is inside the bar group's coverage area
            if (EX <= AX) && (EY >= AY) && (Ex >= Ax) && (Ey <= Ay)

                Areaelemento = (coordnosdoelem(2, 1) - coordnosdoelem(1, 1)) * ...
                               (coordnosdoelem(3, 2) - coordnosdoelem(2, 2));

                espessuraequivalente{i}(k, 1) = VolArm{i}(j) / Areaelemento;
                espessuraequivalente{i}(k, 2) = Dados.Arm.Ei;
            end
        end
    end
end

x = espessuraequivalente{1};
y = espessuraequivalente{2};
TaxArmElem = [x, y];

return
