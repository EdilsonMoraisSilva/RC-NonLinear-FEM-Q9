% =========================================================================
% constitutive_steel.m
%
% Steel constitutive update at a single Gauss point.
%
% Implements the uniaxial Giuffrè-Menegotto-Pinto monotonic envelope
% (Eq. 22 of the paper) to obtain the tangent modulus and current stress
% along the x and y reinforcement directions.
%
% Inputs:
%   Dados  - global problem data (steel parameters in Dados.Arm)
%   cont   - element index
%   c      - Gauss point index within the element
%   e      - 3 x 9 strain matrix of the element (current Gauss point in
%            column 'c'); e(1,:) = exx, e(2,:) = eyy, e(3,:) = gamma_xy
%   ResAco - struct accumulating steel state; ResAco.TA(:,1) and
%            ResAco.TA(:,3) hold the smeared steel thickness in x and y
%
% Outputs (fields of ResAco):
%   ResAco.ET.x(cont,c) - tangent modulus along x reinforcement
%   ResAco.ET.y(cont,c) - tangent modulus along y reinforcement
%   ResAco.ST.x(cont,c) - signed stress along x reinforcement
%   ResAco.ST.y(cont,c) - signed stress along y reinforcement
%
% Notes on bug fixes (relative to the original Master-thesis code):
%   * The original code overwrote ST1 instead of computing ST2 in the
%     biaxial-reinforcement branch; ST2 is now properly assigned.
%   * The original code used the local Gauss-point index 'c' (lowercase)
%     instead of the auxiliary variable 'C' inside the GMP expressions
%     for the uniaxial-reinforcement branches; this is corrected.
%   * The undefined variable 'ST' in the y-reinforcement-only branch
%     (raised a runtime error if e(1,c) > 0) is replaced by ST2.
% =========================================================================

function [ResAco] = constitutive_steel(Dados, cont, c, e, ResAco)

% --- Steel parameters (Giuffrè-Menegotto-Pinto) --------------------------
epy = Dados.Arm.epy;   % yield strain
fpy = Dados.Arm.fpy;   % yield stress
b   = Dados.Arm.b;     % hardening ratio  E_y / E_0
r   = Dados.Arm.r;     % curvature parameter of the elastic-hardening transition
% epf = Dados.Arm.epf;   % ultimate strain (kept for future ductility checks)

% --- Current strain magnitudes at this Gauss point -----------------------
epxx = abs(e(1, c));
epyy = abs(e(2, c));

% Auxiliary common factor in GMP (uses |eps/eps_y|, not Gauss-point index)
Cx = (epxx/epy)^r + 1;
Cy = (epyy/epy)^r + 1;

has_x = ResAco.TA(cont, 1) > 0;
has_y = ResAco.TA(cont, 3) > 0;

if has_x && has_y
    % ----- Reinforcement in both directions -----------------------------
    ET1 = fpy * ( b/epy - (b-1)/(epy*(Cx^(1/r))) ...
                 + epxx*((epxx/epy)^(r-1))*(b-1) / ((epy^2)*(Cx^(1+1/r))) );
    ST1 = fpy * ( b*(epxx/epy) + ((1-b)*(epxx/epy)) / ((1 + (epxx/epy)^r)^(1/r)) );

    ET2 = fpy * ( b/epy - (b-1)/(epy*(Cy^(1/r))) ...
                 + epyy*((epyy/epy)^(r-1))*(b-1) / ((epy^2)*(Cy^(1+1/r))) );
    ST2 = fpy * ( b*(epyy/epy) + ((1-b)*(epyy/epy)) / ((1 + (epyy/epy)^r)^(1/r)) );

    ResAco.ET.x(cont, c) = ET1;
    ResAco.ET.y(cont, c) = ET2;

    if e(1, c) < 0
        ResAco.ST.x(cont, c) = -ST1;
    else
        ResAco.ST.x(cont, c) =  ST1;
    end

    if e(2, c) < 0
        ResAco.ST.y(cont, c) = -ST2;
    else
        ResAco.ST.y(cont, c) =  ST2;
    end

elseif has_x
    % ----- Reinforcement only in x --------------------------------------
    ET1 = fpy * ( b/epy - (b-1)/(epy*(Cx^(1/r))) ...
                 + epxx*((epxx/epy)^(r-1))*(b-1) / ((epy^2)*(Cx^(1+1/r))) );
    ST1 = fpy * ( b*(epxx/epy) + ((1-b)*(epxx/epy)) / ((1 + (epxx/epy)^r)^(1/r)) );

    ResAco.ET.x(cont, c) = ET1;
    ResAco.ET.y(cont, c) = 0;

    if e(1, c) < 0
        ResAco.ST.x(cont, c) = -ST1;
    else
        ResAco.ST.x(cont, c) =  ST1;
    end
    ResAco.ST.y(cont, c) = 0;

elseif has_y
    % ----- Reinforcement only in y --------------------------------------
    ET2 = fpy * ( b/epy - (b-1)/(epy*(Cy^(1/r))) ...
                 + epyy*((epyy/epy)^(r-1))*(b-1) / ((epy^2)*(Cy^(1+1/r))) );
    ST2 = fpy * ( b*(epyy/epy) + ((1-b)*(epyy/epy)) / ((1 + (epyy/epy)^r)^(1/r)) );

    ResAco.ET.x(cont, c) = 0;
    ResAco.ET.y(cont, c) = ET2;

    if e(2, c) < 0
        ResAco.ST.y(cont, c) = -ST2;
    else
        ResAco.ST.y(cont, c) =  ST2;
    end
    ResAco.ST.x(cont, c) = 0;

else
    % ----- No reinforcement in this element -----------------------------
    ResAco.ET.x(cont, c) = 0;
    ResAco.ET.y(cont, c) = 0;
    ResAco.ST.x(cont, c) = 0;
    ResAco.ST.y(cont, c) = 0;
end

return
