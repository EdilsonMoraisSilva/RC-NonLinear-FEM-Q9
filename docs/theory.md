# Theory summary

This document is a short companion to the paper's Section 2. It maps
each component of the formulation to the file in the source tree where
it is implemented.

## 1. Orthotropic incremental constitutive matrix

The incremental orthotropic stress-strain relationship of Darwin and
Pecknold (1977) is written, in the rotated principal-direction frame, as:

```
{dσ_p} = [Dp] {dε_p}
```

with the constitutive matrix in principal directions (Eq. 11 of the paper):

```
[Dp] = 1/(1 - ν²) · [ E1            ν·sqrt(E1·E2)   0
                      ν·sqrt(E1·E2) E2              0
                      0             0               G_p ]
```

where the transverse modulus is taken as

```
G_p = (1/4) · ( (E1 + E2) - 2·ν·sqrt(E1·E2) )
```

A coordinate transformation `[T]` (rotation by `θ`) brings the
constitutive matrix into the global axes:

```
[D] = [T]^T · [Dp] · [T]
```

**File:** `src/element/ElemMembranaNaoLinearFisicaArmDistr.m`,
function `MatrizconstitutivaOrtho`.

## 2. Equivalent uniaxial strain

The equivalent uniaxial strain `eiu` is accumulated incrementally
along the loading path:

```
eiu_i ← eiu_i + (σ_i - σ_i^prev) / E_i^prev
```

For each principal direction `i = 1, 2`, `eiu_i` is the variable that
is fed into the uniaxial-equivalent stress-strain curve (Section 2.3 of
the paper) to obtain the current principal stress `Sp(i)` and tangent
modulus `E(i)`.

**File:** `src/constitutive/ParametrosMatrizconstitutivaMembranaNLF.m`.

## 3. Kupfer biaxial envelope

Four regions are distinguished according to the stress ratio
`α = σ1 / σ2`:

| Region | Range                | Description                  |
|--------|----------------------|------------------------------|
| (a)    | `0 ≤ α ≤ 1`          | Biaxial compression          |
| (b)    | `-0.17 ≤ α ≤ 0`      | Compression / small tension  |
| (c)    | `-∞ ≤ α ≤ -0.17`     | Compression / large tension  |
| (d)    | `1 ≤ α ≤ ∞`          | Biaxial tension              |

For each region the peak stress `σ_ic` and the corresponding peak strain
`ε_ic` are extracted from the envelope formulas of Kupfer et al. (1969)
(Tables 1 and 2 of the paper).

**File:** `src/constitutive/ParametrosMatrizconstitutivaMembranaNLF.m`.

## 4. Equivalent uniaxial stress-strain curves

Four branches of the equivalent uniaxial response are implemented:

| Branch   | Range                       | Equation in paper |
|----------|-----------------------------|-------------------|
| Compression below peak | `\|eiu\| < \|ε_ic\|` | Eq. 14 (Desayi-Krishnan) |
| Compression above peak | `\|eiu\| ≥ \|ε_ic\|` | Eq. 16 (linear softening) |
| Tension below peak     | `\|eiu\| < \|ε_it\|` | Eq. 18 (linear elastic) |
| Tension above peak     | `\|eiu\| ≥ \|ε_it\|` | Eq. 19 (smeared crack)  |

The corresponding tangent moduli are obtained analytically (Eqs. 15, 17,
20 of the paper).

## 5. Smeared rotating crack model

When the tensile strain exceeds the peak (`eiu_i > ε_it`) a smeared
crack is opened perpendicular to the current `σ_i` direction. The
post-peak (softening) branch is regularised using the fracture-energy
criterion:

```
ε_tu = 2 · G_f / (σ_it · h)
```

where `h` is the characteristic length of the integration point,
projected on the crack plane.

The crack state is stored per Gauss point and per principal direction:

- `0` -- intact
- `1` -- active crack (softening branch)
- `2` -- saturated crack (zero stress)

**File:** `src/constitutive/ParametrosFissuracao.m`.

## 6. Smeared reinforcement and steel constitutive law

The steel reinforcement is represented in a smeared (distributed) form
via an equivalent thickness per element and per direction (x, y),
computed from the user-defined bar groups (Section 2.6 of the paper).

The uniaxial Giuffrè-Menegotto-Pinto monotonic envelope is used to
update the steel stress and tangent modulus along each direction:

```
σ = f_py · [ b · ε/ε_py + (1 - b) · ε/ε_py / (1 + (ε/ε_py)^R)^(1/R) ]
```

**Files:**
- `src/preprocessing/Obtemtaxaarmadura.m` (smeared steel thickness)
- `src/constitutive/ParametrosAco.m`        (G-M-P update)

## 7. 9-node isoparametric Q9 element

Plane-stress membrane element with 9 nodes (4 corners + 4 mid-sides +
1 centre) and 2 DOFs per node (u_x, u_y). The full 3 × 3 Gauss rule is
used by default.

The element stiffness matrix sums the concrete and smeared-steel
contributions at each Gauss point (Eq. 32 of the paper):

```
[k] = sum_gp { t_c · B^T · D_c · B · |J| · w_gp
             + t_s · B^T · D_s · B · |J| · w_gp }
```

**File:** `src/element/ElemMembranaNaoLinearFisicaArmDistr.m`.

## 8. Global solution: Newton-Raphson

The incremental-iterative driver in `src/solver/AnaliseElemFinitosNaoLinear.m`
performs, for each load step `np`:

1. Assemble `K` and `P_ext` using current state `(ResConc, ResAco)`.
2. Solve `K · ΔD = P_ext - P_int` for the displacement increment.
3. Update Gauss-point stresses, strains, principal directions and
   tangent moduli at every integration point.
4. Recompute `P_int` and the unbalanced residual; iterate until
   `‖P_ext - P_int‖ / ‖P_ext‖ < tol`.

A sparse global stiffness matrix and the `backslash` direct solver are
used. Quadratic convergence is not strictly observed, due to the
non-smooth envelope-based concrete update (Sec. 6 of the paper).

## References (used in this document)

- Darwin, D., Pecknold, D.A. (1977). *Nonlinear biaxial stress-strain law
  for concrete.* J. Eng. Mech. Div., 103(2), 229-241.
- Kupfer, H., Hilsdorf, H.K., Rüsch, H. (1969). *Behavior of concrete
  under biaxial stresses.* ACI Journal, 66(8), 656-666.
- Desayi, P., Krishnan, S. (1964). *Equation for the stress-strain curve
  of concrete.* ACI Journal, 61, 345-350.
- Giuffrè, A., Pinto, P.E. (1970), and Menegotto, M., Pinto, P.E. (1973).
- Bresler, B., Scordelis, A.C. (1963). *Shear strength of reinforced
  concrete beams.* ACI Journal, 60, 51-74.
