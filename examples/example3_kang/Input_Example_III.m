%==========================================================================
% Universidade Federal do Pará (UFPA)
% Instituto de Tecnologia (ITEC)
% Faculdade de Engenharia Civil
% Author: Edilson Morais Lima e Silva
%
% Input script for the FE nonlinear analysis program (plane stress, Q9
% membrane element, smeared reinforcement).
%
% Example III -- Kang beam (1977)
%   Simply supported RC beam, single concentrated load at mid-span.
%   Geometry : L = 6.40 m, h = 0.55 m, b = 0.228 m
%   Loading  : total P = 311.36 kN, applied in 14 monotonic load steps
%   Reinforcement: three longitudinal layers (one compressive at top,
%                  two tensile layers at the bottom)
%
% Units: kN, m, MPa
%==========================================================================

function Input_Example_III

clc;format short g;close all;

% =========================================================================
% Console banner -- program identification
% =========================================================================
fprintf('\n');
fprintf('==========================================================================\n');
fprintf('=                 Federal University of Para (UFPA)                      =\n');              
fprintf('=                    Technology Institute (ITEC)                         =\n');
fprintf('=               Civil Engineering Graduate Program (PPGEC)               =\n');
fprintf('= Record of revisions                                                    =\n');
fprintf('= Date            Authors                     Discription of Change      =\n');
fprintf('= -----------  ------------------     -----------------------------------=\n');
fprintf('= 14/05/2026   Edilson Morais                                            =\n');
fprintf('==========================================================================\n');
fprintf('=DEVELOPED BY                                                            =\n');
fprintf('=                  Experimental Group on Dynamics,                       =\n');
fprintf('=             Instrumentatio and Singal Processing (GEDIPS)              =\n');
fprintf('==========================================================================\n');
fprintf('=Input script for the FE nonlinear analysis program                      =\n');
fprintf('=(plane stress, Q9 membrane element, smeared reinforcement).             =\n');
fprintf('= Example III -- Kang beam (1977), total load 311.36 kN       =\n');
fprintf('=Units: kN, m, MPa                                                       =\n');
fprintf('==========================================================================\n');
fprintf('\n');


%% ========================================================================
%  1. Solver and problem parameters
%  ========================================================================
% Newton-Raphson controls (number of load steps, tolerance) and concrete
% compressive strength fc.

    % Number of load steps (monotonic, equal-size increments)
    Dados.NumPassos = 14;

    % Convergence tolerance on the relative unbalanced force norm
    Dados.tol = 1e-4;

    % Concrete compressive strength
    Dados.fc = 38.8;   % MPa


%% ========================================================================
%  2. Geometry and mesh generation
%  ========================================================================
% Beam dimensions and inline structured mesh of Q9 elements. The mesh is
% built in seven horizontal strips so that the three reinforcement
% layers can be aligned with element boundaries (bottom tensile group at
% strips 2-4, top compressive group at strip 6).

    % Number of DOFs per node (in-plane membrane: ux, uy)
    Dados.NumGDLNo = 2;

    % Slab edge lengths
    Lx = 6.40;     % beam length      (m)
    Ly = 0.55;     % beam depth       (m)

    % Mesh density (Q9 element divisions)
    nx = 15;       % divisions along x
    ny = 13;       % divisions along y

    % Segment table:
    %               Lx       Ly       Nx   Ny
    Segment =  [   6.4   0.06400    nx    2     % bottom cover
                   6.4   0.02850    nx    1     % lower tensile layer (3 phi 28.5 mm)
                   6.4   0.03200    nx    1     % spacing between bottom layers
                   6.4   0.02850    nx    1     % upper tensile layer (2 phi 28.5 mm)
                   6.4   0.33350    nx    5     % web
                   6.4   0.01270    nx    1     % top compressive layer (2 phi 12.7 mm)
                   6.4   0.05080    nx    2 ];  % top cover

[Dados] = mesh_generator(Dados, Segment);

%% ========================================================================
%  3. Element properties
%  ========================================================================
% Per-element properties are assigned uniformly. The connectivity matrix
% follows the 9-node ordering required by the element routines (see
% diagram below).

    % Total nodes and DOFs
    Dados.NumNos    = size(Dados.coords, 1);
    Dados.NumGdlEst = Dados.NumNos * Dados.NumGDLNo;


%% ========================================================================
%  4. Boundary conditions and applied loads
%  ========================================================================
% Identification of restrained nodes (supports) and node where the
% external load is applied. Single concentrated load at mid-span top.

    % Nodal restraint flags: 0 = free, 1 = restrained
    Dados.restr = zeros(Dados.NumNos, Dados.NumGDLNo);

    % Position of restraints (simply supported beam)
    %           x       y
    PosRestr = [ 0       0     % pinned support  (left)
                 6.4     0 ];  % roller support  (right)

    cont = 1;
    for i = 1:Dados.NumNos
        if (PosRestr(cont,1) == Dados.coords(i,1)) && ...
           (PosRestr(cont,2) == Dados.coords(i,2))
            no(cont) = i;
            cont = cont + 1;
        end
        if cont > size(PosRestr, 1)
            break
        end
    end

    % Apply restraints: pinned (ux, uy) on the left, roller (uy) on the right
    Dados.restr(no(1), :) = [1 1];
    Dados.restr(no(2), 2) =  1;

    % ---------------------------------------------------------------------
    % Build the Q9 connectivity matrix.
    %
    % The 9-node ordering required by the element routines is shown below.
    % Throughout the program, the shape functions are evaluated according
    % to this pattern, so the order MUST be preserved.
    %
    %      2-----------5-----------1
    %      |                       |
    %      |                       |
    %      |                       |
    %      6           9           8
    %      |                       |
    %      |                       |
    %      |                       |
    %      3-----------7-----------4
    %
    % The indices below are derived from the structured nodal grid
    % (2*nx+1 nodes per row).
    % ---------------------------------------------------------------------
    k = 1;
    for i = 1:ny
        for j = 1:nx
            Noi = 3 + (j-1)*2 + 2*i*(2*nx + 1);
            Noj = Noi - 2;
            Nok = Noj - (4*nx + 2);
            Nol = Nok + 2;
            Nom = Noi - 1;
            Non = Noj - (2*nx + 1);
            Noo = Nok + 1;
            Nop = Nol + (2*nx + 1);
            Noq = Nop - 1;

            Dados.conect(k, :) = [Noi Noj Nok Nol Nom Non Noo Nop Noq];
            k = k + 1;
        end
    end

    % Initial concrete tangent modulus (reference value; the value used
    % per element is set explicitly in Dados.ele below).
    E = (5600 * sqrt(Dados.fc)) * 1000;

    % Number of elements and per-element properties
    Dados.NumElem = size(Dados.conect, 1);
    for i = 1:length(Dados.conect)
        % Columns:   E0[kPa]      nu     t[m]    GaussOrd   bz   Mx   My
        Dados.ele(i, :)       = [3355.8e4  0.2  0.228   3        0   0   0];
        Dados.TipoElemento(i, :) = 'ElemMembrana';
    end

    % --- Applied load (single concentrated load at mid-span) ------------
    Dados.CargaNo = zeros(Dados.NumNos, Dados.NumGDLNo);
    Dados.P       = -311.36;        % kN, downward
    PosCarga      = [3.20 0.55];    % mid-span, top face

    for i = 1:Dados.NumNos
        if (PosCarga(1,1) == Dados.coords(i,1)) && ...
           (PosCarga(1,2) == Dados.coords(i,2))
            noCarga = i;
            break
        end
    end
    Dados.CargaNo(noCarga, 2) = Dados.P;

    % --- Monitored node for the load-displacement curve -----------------
    % Bottom face, mid-span: vertical displacement compared with the
    % experimental curve of Kang (1977).
    PosDeslNo = [3.20 0];
    for i = 1:Dados.NumNos
        if (PosDeslNo(1,1) == Dados.coords(i,1)) && ...
           (PosDeslNo(1,2) == Dados.coords(i,2))
            NoDesl = i;
            break
        end
    end
    Dados.No     = NoDesl;
    Dados.GDLNo  = 2;       % vertical direction


%% ========================================================================
%  5. Reinforcement (smeared representation)
%  ========================================================================
% Bar groups for the smeared steel model. Each row of Dados.Arm.Bars{i}
% defines a group of bars in direction i (1 = x, 2 = y) with:
%   [ n_bars, diameter, length, x1, x2, y1, y2, z1, z2 ]
% The smeared steel equivalent thickness per element is computed later
% by Obtemtaxaarmadura.
%
% NOTE: this example uses the field name 'Bars' (in English). The
% Examples I and II use 'Barras'. If Obtemtaxaarmadura reads the
% Portuguese name, the field below must be renamed to 'Barras'
% (or the routine adapted to read either spelling).

    Dados.TipoArm = 'Distributed';

    switch Dados.TipoArm
        case 'Distributed'

            % Three longitudinal layers in x:
            %   - 2 phi 12.7 mm at the top  (compressive reinforcement)
            %   - 2 phi 28.5 mm  upper tensile layer
            %   - 3 phi 28.5 mm  lower tensile layer
            %
            %                       n.bars   diameter   L     x1   x2    y1      y2      z1   z2
            Dados.Arm.Barras =  { [    2       0.01270   Lx   0.0  Lx   0.4866  0.4991  0.0  0.0
                                     2       0.02850   Lx   0.0  Lx   0.1246  0.1520  0.0  0.0
                                     3       0.02850   Lx   0.0  Lx   0.0650  0.0924  0.0  0.0 ]

                                % Reinforcement in y direction (none in this example)
                                [    0       0.00      0    0.0  0.0   0.0    0.0    0.0  0.0
                                     0       0.00      0    0.0  0.0   0.0    0.0    0.0  0.0 ] };

        case 'Discretizada'
            % Reserved for future extension (discrete bars + bond-slip).
    end

    % ---------------------------------------------------------------------
    % Plot the undeformed mesh
    % ---------------------------------------------------------------------
    escala               = 0.06;
    Dados.InicioSeta     = 'Top-down';   % 'Bottom-up' | 'Top-down'
    Dados.Desenho.LabelNo = 'no';
    Dados.Desenho.LabelEl = 'no';

    DesenhaMalhaIndefQ9(Dados, escala);


%% ========================================================================
%  6. Steel constitutive parameters (Giuffre-Menegotto-Pinto)
%  ========================================================================
% Uniaxial monotonic envelope parameters used by ParametrosAco.

    Dados.Arm.fpy = 55.2e4;        % Steel yield stress       (kPa)
    Dados.Arm.Ei  = 21167.6e4;     % Initial elastic modulus  (kPa)
    Dados.Arm.fpu = 68.9e4;        % Ultimate stress          (kPa)
    Dados.Arm.Ey  = 288.2e4;       % Hardening tangent modulus (kPa)

    % Hardening ratio b = Ey / E0
    Dados.Arm.b   = Dados.Arm.Ey / Dados.Arm.Ei;

    % Curvature parameter R defining the transition between the elastic
    % and hardening branches of the GMP envelope (Eq. 22 of the paper).
    Dados.Arm.r   = 100;

    % Yield strain and ultimate strain
    Dados.Arm.epy = Dados.Arm.fpy / Dados.Arm.Ei;
    Dados.Arm.epf = ((Dados.Arm.fpu - Dados.Arm.fpy) / Dados.Arm.Ei) + Dados.Arm.epy;


%% ========================================================================
%  7. Console banner -- analysis parameters summary
%  ========================================================================
% Print mesh size, materials and solver settings before launching the
% Newton-Raphson driver.

    NumElem    = size(Dados.conect, 1);
    Numnoselem = size(Dados.conect, 2);

    fprintf('==========================================================================\n');
    fprintf('  RC Nonlinear FE Analysis  --  Plane stress, Q9 element, smeared steel\n');
    fprintf('==========================================================================\n');
    fprintf('  Mesh\n');
    fprintf('    Number of nodes       : %d\n', Dados.NumNos);
    fprintf('    Number of elements    : %d\n', NumElem);
    fprintf('    Nodes per element     : %d\n', Numnoselem);
    fprintf('    DOFs per node         : %d\n', Dados.NumGDLNo);
    fprintf('    Total DOFs            : %d\n', Dados.NumGdlEst);
    fprintf('  Concrete\n');
    fprintf('    fc                    : %8.2f MPa\n', Dados.fc);
    fprintf('    E0                    : %8.0f MPa\n', Dados.ele(1,1)/1e3);
    fprintf('    nu                    : %8.3f\n',     Dados.ele(1,2));
    fprintf('    Thickness             : %8.3f m\n',   Dados.ele(1,3));
    fprintf('  Steel (Giuffre-Menegotto-Pinto)\n');
    fprintf('    fy                    : %8.2f MPa\n', Dados.Arm.fpy/1e3);
    fprintf('    E0                    : %8.0f MPa\n', Dados.Arm.Ei /1e3);
    fprintf('    Hardening ratio b     : %8.4f\n',     Dados.Arm.b);
    fprintf('    Curvature R           : %8.1f\n',     Dados.Arm.r);
    fprintf('  Solver\n');
    fprintf('    Load steps            : %d\n',   Dados.NumPassos);
    fprintf('    Convergence tolerance : %.1e\n', Dados.tol);
    fprintf('    Monitored DOF (node)  : %d (dir %d)\n', Dados.No, Dados.GDLNo);
    fprintf('==========================================================================\n');
    fprintf('\n');
    fprintf(' Step |  Iters |  Cum.iters |  Final error |   Load (kN)   |  Disp (mm) \n');
    fprintf('------+--------+------------+--------------+---------------+------------\n');


%% ========================================================================
%  8. Run nonlinear FE analysis
%  ========================================================================
% Calls the Newton-Raphson driver. Returns nodal displacements,
% reactions, DOF table, and the updated Dados struct (which carries
% iteration counters and history vectors used in the report below).

    tAnalysis = tic;

    [Resultado, Pr, D, GDL, Dados] = solver_non_linear(Dados);


%% ========================================================================
%  9. Console summary -- analysis report
%  ========================================================================
    elapsed = toc(tAnalysis);

    fprintf('------+--------+------------+--------------+---------------+------------\n');
    fprintf('\n');
    fprintf('==========================================================================\n');
    fprintf('  Analysis completed\n');
    fprintf('==========================================================================\n');
    fprintf('    Total load steps          : %d\n',          Dados.NumPassos);
    fprintf('    Total Newton iterations   : %d\n',          Dados.IteracoesTotais);
    fprintf('    Avg iterations per step   : %.2f\n',        Dados.IteracoesTotais / Dados.NumPassos);
    fprintf('    Final monitored load      : %12.3f kN\n',   abs(Dados.loadStep));
    fprintf('    Final monitored disp.     : %12.4f mm\n',   abs(Dados.dispStep)*1000);
    fprintf('    Stiffness ratio (sec/init): %12.4f\n',      abs(Dados.loadStep)/max(abs(Dados.dispStep),eps) / ...
                                                            ( abs(Dados.CargaNoC(2,1)) / max(abs(Resultado.DeslocNoC(2,1)),eps) ));
    fprintf('    Elapsed wall-clock time   : %12.2f s\n',    elapsed);
    fprintf('==========================================================================\n\n');

   
    
    % --- Cracking statistics across all Gauss points --------------------
    nCracked   = 0;
    nSaturated = 0;
    nGPTotal   = 0;
    for e = 1:NumElem
        if isfield(Resultado.el(e), 'TensoesPrinc') && ~isempty(Resultado.el(e).TensoesPrinc)
            sp = Resultado.el(e).TensoesPrinc;
            nGPTotal   = nGPTotal   + size(sp, 2);
            nCracked   = nCracked   + sum(sp(5,:) == 1 | sp(6,:) == 1);
            nSaturated = nSaturated + sum(sp(5,:) == 2 | sp(6,:) == 2);
        end
    end
    if nGPTotal > 0
        fprintf('  Cracking state at end of analysis (per Gauss point):\n');
        fprintf('    Intact      : %5d  (%5.1f%%)\n', nGPTotal - nCracked - nSaturated, 100*(nGPTotal - nCracked - nSaturated)/nGPTotal);
        fprintf('    Cracked     : %5d  (%5.1f%%)\n', nCracked, 100*nCracked/nGPTotal);
        fprintf('    Total GPs   : %5d\n\n',          nGPTotal);
    end


%% ========================================================================
% 10. Post-processing
%  ========================================================================
% Plot the deformed configuration and the principal stress field.
% Optional results printout to a text file is left commented for
% reference.

    % --- Deformed configuration -----------------------------------------
    EscalaDesloc = 30;
    ElementMembranaDeformada(Dados, Resultado, EscalaDesloc);

    % --- Optional: write results to a file ------------------------------
    % ImprResAnaliseTensoes(Dados, 'ResultadoElemMembrana.m', Resultado, Pr, D);

    % --- Principal stress field on the deformed mesh --------------------
    escala = 0.03;
    DesenhaTensoesprincipais(Dados, Resultado, escala);

end