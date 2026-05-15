% =========================================================================
% cracking_model.m
%
% Computes the ultimate tensile strain etu used in the post-peak (softening)
% branch of the smeared crack model, based on the fracture-energy criterion
% described in Section 2.5 of the paper (Eq. 20):
%
%       etu = 2 * Gf / (sic_i * h)
%
% where Gf is the fracture energy associated with the area under the
% uniaxial tensile stress-strain curve, h is the characteristic length of
% the integration point (taken as the element projection on the crack
% plane), and sic_i is the peak tensile stress in the current principal
% direction.
%
% This regularises the post-peak response and reduces mesh-size dependence.
% =========================================================================

function [etu] = cracking_model(Dados, DadosConc, cont, c, sic, eic, eiu, i) 

% Element edge lengths
B = Dados.coords(Dados.conect(cont, 2), 1) - Dados.coords(Dados.conect(cont, 1), 1);
H = Dados.coords(Dados.conect(cont, 3), 2) - Dados.coords(Dados.conect(cont, 2), 2);

% Characteristic crack length projected on the principal direction
h = B*cos(DadosConc.el(cont).Sp(4, c)) + H*sin(DadosConc.el(cont).Sp(4, c));

% Area under the uniaxial tension curve (energy per unit volume)
gf = sic(i) * eiu;

% Fracture energy per unit area
Gf = h * gf;

% Ultimate strain
etu = (2 * Gf) / (sic(i) * h);

return
