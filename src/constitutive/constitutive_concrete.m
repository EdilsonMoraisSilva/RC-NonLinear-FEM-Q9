% =========================================================================
% constitutive_concrete.m
%
% Concrete constitutive update at a single Gauss point.
%
% Implements the framework described in Section 2 of the paper:
%   - Kupfer biaxial failure envelope (Section 2.4, Tables 1-2)
%   - Equivalent uniaxial strain concept (Section 2.2-2.3)
%   - Smeared rotating crack model with fracture-energy regularisation
%     (Section 2.5, Eq. 21)
%
% For each principal direction i = 1, 2 the routine:
%   1. Computes the stress ratio alpha = sigma_1/sigma_2 and selects the
%      corresponding region (a, b, c, d) of the Kupfer envelope.
%   2. Extracts the peak stress sic and peak strain eic for that direction.
%   3. Increments the accumulated equivalent uniaxial strain eiu(i).
%   4. Selects the active branch of the equivalent uniaxial stress-strain
%      curve (Eqs. 14-19) and updates the principal stress Sp(i) and the
%      tangent modulus E(i).
%   5. If tension exceeds the peak, activates a smeared crack and records
%      the crack opening w_i and crack state (1=cracking, 2=saturated).
%
% Output struct ResConc carries:
%   ResConc.Sp(1:2)        : updated principal stresses
%   ResConc.E(1:2)         : updated tangent moduli
%   ResConc.eiu{cont}(i,c) : accumulated equivalent uniaxial strain history
%   ResConc.Fissura(i)     : 0=intact | 1=cracking | 2=saturated
%   ResConc.AberturFissura : crack opening [m]
%   ResConc.h_caract       : characteristic length per direction [m]
%   ResConc.nu             : effective Poisson coupling parameter
% =========================================================================

function [ResConc] = constitutive_concrete(Dados, DadosConc, cont, c, ResConc)

% --- Current principal stresses at this Gauss point (from previous it.) --
sigma1 = DadosConc.el(cont).Sp(1, c);
sigma2 = DadosConc.el(cont).Sp(2, c);

a     = sigma1 / sigma2;
sigma = [sigma1, sigma2];

% --- Material parameters --------------------------------------------------
fc = Dados.fc * 10^3;     % Concrete compressive strength [kPa]
ft = 0.1 * fc;            % Concrete tensile strength    [kPa]
E0 = Dados.ele(cont, 1);  % Initial tangent modulus      [kPa]
ec = 0.002;               % Strain at peak compression
eu = 0.0038;              % Ultimate compressive strain

eiu = ResConc.eiu{cont}(:, c);

% =========================================================================
% Kupfer biaxial envelope (Tables 1 and 2 of the paper)
%   Region (a): 0     <= alpha <= 1     -- biaxial compression
%   Region (b): -0.17 <= alpha <= 0     -- compression / small tension
%   Region (c): -inf  <= alpha <= -0.17 -- compression / large tension
%   Region (d): 1     <= alpha <= inf   -- biaxial tension
% =========================================================================
if 0 <= a && a <= 1
    s2c = ((1 + 3.65*a) / ((1 + a)^2)) * fc;
    s1c = a * s2c;
    p1  = s1c / fc;
    p2  = s2c / fc;
    e1c = ec * (-1.6*p1^3 + 2.25*p1^2 + 0.35*p1);
    e2c = ec * (3*p2 - 2);
    sic = [s1c, s2c];
    eic = [e1c, e2c];

elseif -0.17 <= a && a <= 0
    s2c = ((1 + 3.28*a) / (1 + a)^2) * fc;
    s1t = abs(a) * s2c;
    p2  = s2c / fc;
    e1t = s1t / E0;
    e2c = ec * (-2.58*p2^3 + 7.54*p2^2 - 8.38*p2 + 4.42);
    sic = [s1t, s2c];
    eic = [e1t, e2c];

elseif -inf <= a && a <= -0.17
    s2c = 0.65 * fc;
    s1t = ft;
    p2  = s2c / fc;                                                  
    e1t = s1t / E0;
    e2c = ec * (-2.58*p2^3 + 7.54*p2^2 - 8.38*p2 + 4.42);
    sic = [s1t, s2c];
    eic = [e1t, e2c];

else
    s2t = ft;
    s1t = s2t;
    e1t = ft / E0;
    e2t = e1t;
    sic = [s1t, s2t];
    eic = [e1t, e2t];
end

% --- Element dimensions for the characteristic length --------------------
B_elem    = Dados.coords(Dados.conect(cont, 2), 1) - Dados.coords(Dados.conect(cont, 1), 1);
H_elem    = Dados.coords(Dados.conect(cont, 3), 2) - Dados.coords(Dados.conect(cont, 2), 2);
theta_pt  = DadosConc.el(cont).Sp(4, c);

% Reset cracking flags for this Gauss point
ResConc.Fissura        = [0 0];
ResConc.AberturFissura = [0 0];
ResConc.h_caract       = [0 0];

% =========================================================================
% Loop over the two principal directions
% =========================================================================
for i = 1:2

    Satual    = sigma(i);
    Santerior = DadosConc.SpA(i, c);
    Eanterior = DadosConc.EpA(i, c);

    % Equivalent uniaxial strain increment (Eq. 12-13)
    deltaeu(i, 1) = (Satual - Santerior) / Eanterior;
    eiu(i, 1)     = eiu(i, 1) + deltaeu(i, 1);
    ResConc.eiu{cont}(i, c) = eiu(i, 1);

    Es = abs(sic(i) / eic(i));

    if sigma(i) < 0
        % ---------- COMPRESSION -----------------------------------------
        if sigma(i) <= 0.3 * fc
            % Linear elastic stage
            eiua          = abs(eiu(i, 1));
            ResConc.Sp(i) = -(E0 * eiua);
            ResConc.E(i)  = E0;
        else
            if abs(eiu(i, 1)) < abs(eic(i))
                % Case 1 -- ascending branch (Eq. 14, 15)
                eiua          = abs(eiu(i, 1));
                Sp_val        = (eiua * E0) / (1 + ((E0/Es) - 2)*(eiua/eic(i)) + (eiua/eic(i))^2);
                ResConc.Sp(i) = -Sp_val;
                m             = E0 * (1 - (eiua/eic(i)));
                b             = 1 + ((E0/Es) - 2)*(eiua/eic(i)) + (eiua/eic(i))^2;
                ResConc.E(i)  = m / (b^2);
            else
                % Case 2 -- descending branch (Eq. 16, 17)
                eiua          = abs(eiu(i, 1));
                Sp_val        = sic(i) + ((0.2*fc - sic(i)) / (eu - eic(i))) * (eiua - eic(i));
                ResConc.Sp(i) = -Sp_val;
                ResConc.E(i)  = abs((sic(i) - 0.2*fc) / (eic(i) - eu));
            end
        end

    else
        % ---------- TENSION ---------------------------------------------
        if abs(eiu(i, 1)) < abs(eic(i))
            % Case 3 -- linear tension below peak (Eq. 18)
            eiua          = abs(eiu(i, 1));
            ResConc.Sp(i) = E0 * eiua;
            ResConc.E(i)  = E0;

        else
            % Case 4 -- tension above peak: active crack (Eq. 19, 20)
            eiua = abs(eiu(i, 1));

            % Characteristic length projected on the crack direction.
            % Crack is normal to sigma_i, so direction i=1 uses theta+pi/2;
            % direction i=2 uses theta.
            if i == 1
                h_i = abs(B_elem*cos(theta_pt + pi/2)) + abs(H_elem*sin(theta_pt + pi/2));
            else
                h_i = abs(B_elem*cos(theta_pt))        + abs(H_elem*sin(theta_pt));
            end
            if h_i < 1e-12
                h_i = sqrt(B_elem^2 + H_elem^2);
            end

            % Ultimate strain via fracture energy (Eq. 21)
            [etu] = cracking_model(Dados, DadosConc, cont, c, sic, eic, eiua, i);

            % Post-peak stress and tangent modulus
            Sp_val        = (sic(i) / (etu - eic(i))) * (etu - eiua);
            Sp_val        = max(Sp_val, 0);  % no negative post-peak stress
            ResConc.Sp(i) = Sp_val;
            ResConc.E(i)  = abs(sic(i) / (eic(i) - etu));

            % --- Crack opening and crack-state flags ------------------
            w_i = (eiua - eic(i)) * h_i;     % crack opening [m]

            if eiua >= etu
                ResConc.Fissura(i) = 2;      % saturated crack (stress = 0)
            else
                ResConc.Fissura(i) = 1;      % active softening
            end
            ResConc.AberturFissura(i) = max(w_i, 0);
            ResConc.h_caract(i)       = h_i;
        end
    end
end

% --- Effective Poisson-like coupling (Section 2.1) ----------------------
% Darwin & Pecknold (1977): under tension-compression states the
% effective Poisson coupling is enhanced by the compressive direction.
% Each term uses its own principal-stress component.
nu0 = 0.2;
if (ResConc.Sp(1) < 0 && ResConc.Sp(2) > 0) || (ResConc.Sp(1) > 0 && ResConc.Sp(2) < 0)

    % Pick the compressive principal stress for the enhancement term
    if ResConc.Sp(1) < 0
        sc = ResConc.Sp(1);
        st = ResConc.Sp(2);
    else
        sc = ResConc.Sp(2);
        st = ResConc.Sp(1);
    end

    % Normalised ratios, clamped to avoid the (sigma/fc)^4 term blowing up
    % when the local stress overshoots the strength envelope inside a
    % Newton-Raphson iteration.
    rc = min(abs(sc) / fc, 1.0);
    rt = min(abs(st) / fc, 1.0);

    nu = nu0 + 0.6 * rc^4 + 0.4 * rt^4;

    % Cap the effective Poisson coupling to keep (1 - nu^2) safely
    % positive in the orthotropic constitutive matrix (Eq. 11 of the paper).
    nu = min(nu, 0.495);

    ResConc.nu = nu;
else
    ResConc.nu = nu0;
end

return
