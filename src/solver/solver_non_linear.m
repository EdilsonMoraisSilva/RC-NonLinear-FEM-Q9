% =========================================================================
% solver_non_linear  -  Main nonlinear FE driver
%
% Universidade Federal do Pará (UFPA)
% Instituto de Tecnologia (ITEC)
% Faculdade de Engenharia Civil
% Author: Edilson Morais Lima e Silva
% Co-author: Remo Magalhães de Souza
%
% Top-level routine that performs incremental-iterative nonlinear analysis
% of reinforced concrete plane-stress structures with smeared reinforcement.
%
% Algorithm (per load step):
%   1. Assemble the global tangent stiffness matrix K (concrete + steel)
%      from element contributions integrated over Gauss points.
%   2. Compute the external load increment.
%   3. Newton-Raphson loop on the unbalanced force vector P_un = P_ext - P_int.
%   4. At each iteration, update equivalent uniaxial strain, principal
%      directions, tangent moduli (concrete and steel) and cracking state
%      at every Gauss point.
%
% NOTE: Only the distributed (smeared) reinforcement representation is
%       included in this release. The discrete-bar variant is reserved
%       for future extensions (steel-concrete bond modelling).
%
% INPUT:
%   Dados  - struct with mesh, materials, loads, BCs and solver settings
%
% OUTPUT:
%   Resultado - struct with displacements, stresses, strains, principal
%               stresses and cracking data at all integration points
%   Pr        - reaction force vector (restrained DOFs)
%   D         - full nodal displacement vector
%   GDL       - DOF numbering table
% =========================================================================

function [Resultado, Pr, D, GDL,Dados] = solver_non_linear(Dados)

% --- Mesh sizes -----------------------------------------------------------
NumElem    = size(Dados.conect, 1);
Numnoselem = size(Dados.conect, 2);

% --- DOF numbering --------------------------------------------------------
[GDL, NGDLR, GDLll] = obtgdlest(Dados);

% --- Initialise per-element state variables -------------------------------
% ResConc.eiu{e}(i, gp) : accumulated equivalent uniaxial strain at
%                        Gauss point 'gp' along principal direction 'i'
%                        for element 'e'.
ResConc = [];
ResAco  = [];
for h = 1:NumElem
    ResConc.eiu{h} = zeros(2, 9);
    Resultado.el(h).Sigma  = zeros(3, 9);
    Resultado.el(h).Epslon = zeros(3, 9);
end

DadosConc = [];
Resultado.S = zeros(3, 9);
Dados.IteracoesTotais = 0;

% Free / restrained displacement vectors
Dl = zeros(GDLll, 1);
Pr = zeros(NGDLR, 1);

% Element-level incremental displacement vector
d = zeros(Numnoselem * Dados.NumGDLNo, 1);

% Reinforcement ratio per element (smeared model)
[TaxArmElem] = reinforcement_setup(Dados);
ResAco.TA    = TaxArmElem;

% =========================================================================
% Load step loop (incremental-iterative Newton-Raphson)
% =========================================================================
for np = 1:Dados.NumPassos

    erro = 1;
    Dados.iteracoes = 0;

    while (erro > Dados.tol)

        % --- Allocate global arrays for this iteration --------------------
        % Sparse assembly is used to keep memory footprint low even for
        % refined meshes. The element stiffness matrices are dense (18x18).
        K        = sparse(Dados.NumGdlEst, Dados.NumGdlEst);
        Peq      = zeros(Dados.NumGdlEst, 1);
        Pinterno = zeros(Dados.NumGdlEst, 1);

        % --- Element loop: assemble K and equivalent nodal forces --------
        for i = 1:NumElem

            Nosdoelem = Dados.conect(i, :);
            coordsnoselem = Dados.coords(Nosdoelem, :);
            cont = i;

            % Distributed (smeared) reinforcement representation
            [var] = Quadratic_Q9_Membrane_element( ...
                'RigidPeq', Dados, d, coordsnoselem, cont, Numnoselem, ...
                ResConc, DadosConc, Resultado, ResAco);

            k_e   = var.k;
            peq_e = var.peq;

            % Element DOF map
            [GDLelem] = obtemgdlelem(Dados, Nosdoelem, GDL, Numnoselem);
            NumGdlElem = length(GDLelem);

            % Scatter element contribution into global arrays
            % (vectorized form replaces the original double for-loop)
            K(GDLelem, GDLelem) = K(GDLelem, GDLelem) + k_e;
            Peq(GDLelem)        = Peq(GDLelem)        + peq_e;

            % Store cartesian coordinates of Gauss points for postprocessing
            Resultado.el(i).CoordsPtsGauus = var.coordPtGauss;
        end

        % --- Partition K into free / restrained blocks --------------------
        Kll = K(1:GDLll,            1:GDLll);
        Krr = K(GDLll+1:Dados.NumGdlEst, GDLll+1:Dados.NumGdlEst);
        Klr = K(1:GDLll,            GDLll+1:Dados.NumGdlEst);
        Krl = K(GDLll+1:Dados.NumGdlEst, 1:GDLll);

        % Free nodal forces and restrained displacements
        [pl, dr] = ObtVetorpLdr(Dados, GDL, NGDLR);

        Peql = Peq(1:GDLll, 1);

        PL = Peql + pl;
        PextlivreMax = PL;

        % External load vector at current step (proportional loading)
        Pexternolivre = PextlivreMax * np / Dados.NumPassos;

        % Pinternolivre only exists after the first global iteration
        if Dados.IteracoesTotais > 0
            Pdesbl = Pexternolivre - Pinternolivre;
        else
            Pdesbl = Pexternolivre;
            Pinternolivre = zeros(GDLll, 1);
        end

        % --- Solve for free DOF displacement increment -------------------
        deltaDl = Kll \ (Pdesbl - Klr*dr);
        deltaPr = Krl*deltaDl + Krr*dr;

        Dl = Dl + deltaDl;
        Pr = Pr + deltaPr;

        % Assemble full displacement vector (free DOFs then restrained)
        deltaD = [deltaDl; dr];
        D      = [Dl;      dr];
        Resultado.D = D;

        % Assemble full force vector
        P = [PL; Pr];

        % =================================================================
        % Recover stresses, strains and update internal state at every
        % Gauss point. This also feeds back the constitutive history
        % (ResConc, ResAco) needed by the next iteration's stiffness.
        % =================================================================
        Pinternolivre = zeros(GDLll, 1);

        for i = 1:NumElem

            Nosdoelem = Dados.conect(i, :);
            coordsnoselem = Dados.coords(Nosdoelem, :);
            [GDLelem] = obtemgdlelem(Dados, Nosdoelem, GDL, Numnoselem);
            NumGdlElem = length(GDLelem);

            % Element incremental displacement vector
            for j = 1:NumGdlElem
                d(j) = deltaD(GDLelem(j));
            end

            cont = i;

            [var] = Quadratic_Q9_Membrane_element( ...
                'Resultado', Dados, d, coordsnoselem, cont, Numnoselem, ...
                ResConc, DadosConc, Resultado, ResAco);

            ResAco    = var.ResAco;
            DadosConc = var.DadosConc;
            ResConc   = var.ResConc;

            % Stress / strain fields for postprocessing
            Resultado.el(i).Sigma        = var.s;
            Resultado.el(i).Epslon       = var.e;
            Resultado.el(i).TensoesPrinc = var.sp;
            Resultado.DadosConc          = DadosConc;

            % Scatter internal forces into global vector
            for g = 1:NumGdlElem
                kk = GDLelem(g);
                Pinterno(kk, 1) = Pinterno(kk, 1) + var.Pintelemento(g);
            end
        end

        % Free part of the internal force vector
        for i = 1:Dados.NumNos
            for j = 1:Dados.NumGDLNo
                if Dados.restr(i, j) == 0
                    g = GDL(i, j);
                    Pinternolivre(g, 1) = Pinterno(g, 1);
                end
            end
        end

        % --- Convergence check --------------------------------------------
        Pdesblivre = Pexternolivre - Pinternolivre;
        erro = norm(Pdesblivre / norm(Pexternolivre));

        Dados.iteracoes       = Dados.iteracoes       + 1;
        Dados.IteracoesTotais = Dados.IteracoesTotais + 1;

        if Dados.iteracoes > 1000
            disp('Iteration limit exceeded (>1000)');
            break
        end
    end

% Monitored DOF: current load and displacement (in kN and mm)
Dados.loadStep  = sum(Pexternolivre);
Dados.dispStep  = D( GDL(Dados.No, Dados.GDLNo) );
itersStep = Dados.iteracoes;

fprintf(' %4d | %6d | %10d |   %.3e  | %12.3f  | %10.4f\n', ...
        np, itersStep, Dados.IteracoesTotais, erro, ...
        abs(Dados.loadStep), abs(Dados.dispStep)*1000);

    % --- Per-step displacement and load at monitored DOF ------------------
    for i = 1:Dados.NumNos
        for j = 1:Dados.NumGDLNo
            Resultado.DeslocNosC(i, j) = D(GDL(i, j));
        end
    end

    DeslocNoC(1, 1)      = 0;
    Dados.CargaNoC(1, 1)       = 0;
    Resultado.DeslocNoC(np+1, 1)   = Resultado.DeslocNosC(Dados.No, Dados.GDLNo);
    Dados.CargaNoC(np+1, 1)    = sum(Pexternolivre);
end

% =========================================================================
% Final post-processing: pack nodal displacements and reactions
% =========================================================================
for i = 1:Dados.NumNos
    for j = 1:Dados.NumGDLNo
        DeslocNos(i, j) = D(GDL(i, j));
        ForcasNos(i, j) = P(GDL(i, j));
    end
end

Resultado.DeslocNos = DeslocNos;
Resultado.ForcasNos = ForcasNos;

% Load-displacement curve at the monitored DOF
figure(7);
plot(abs(Resultado.DeslocNoC), abs(Dados.CargaNoC), 'bs-', 'LineWidth', 2, 'MarkerSize', 10);
grid on;
xlabel('Displacement (m)');
ylabel('Load (kN)');
title('Load-displacement curve at monitored node');




return


% =========================================================================
% Internal helpers
% =========================================================================

% -------------------------------------------------------------------------
% obtgdlest -- builds the DOF numbering table:
%   free DOFs are numbered first (1 .. GDLll),
%   restrained DOFs come last  (GDLll+1 .. NumGdlEst).
% -------------------------------------------------------------------------
function [GDL, NGDLR, GDLll] = obtgdlest(Dados)

con = 0;
LL  = 0;
for no = 1:Dados.NumNos
    for i = 1:Dados.NumGDLNo
        if Dados.restr(no, i) == 0
            con = con + 1;
            GDL(no, i) = con;
            LL = LL + 1;
        end
    end
end
GDLll = LL;

RR = 0;
for no = 1:Dados.NumNos
    for i = 1:Dados.NumGDLNo
        if Dados.restr(no, i) == 1
            con = con + 1;
            GDL(no, i) = con;
            RR = RR + 1;
        end
    end
end
NGDLR = RR;

return

% -------------------------------------------------------------------------
% ObtVetorpLdr -- assembles the free nodal-force vector and the prescribed
% displacement vector (zero by default, kept for future extensions).
% -------------------------------------------------------------------------
function [pl, dr] = ObtVetorpLdr(Dados, GDL, NGDLR)

pl = zeros(0, 0);
for m = 1:Dados.NumNos
    for n = 1:Dados.NumGDLNo
        if Dados.restr(m, n) == 0
            pl(GDL(m, n), 1) = Dados.CargaNo(m, n);
        end
    end
end

dr = zeros(NGDLR, 1);

return
