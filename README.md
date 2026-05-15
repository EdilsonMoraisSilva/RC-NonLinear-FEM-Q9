# Orthotropic FE Framework for Nonlinear Analysis of Reinforced Concrete

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2018a%2B-blue)](https://www.mathworks.com/products/matlab.html)

MATLAB implementation of an orthotropic finite element framework for
monotonic nonlinear analysis of reinforced concrete structures under
plane stress. This repository accompanies the paper:

> Silva, E., Lobo, L., Barreto, D., de Souza, R. (2026).
> **Formulation and implementation of an orthotropic finite element for
> material nonlinear analysis of reinforced concrete structures considering
> cracking under plane stress state.**
> *Computers and Concrete* (under review).
> DOI: *to be added upon acceptance*.

The framework integrates well-established constitutive ingredients (Darwin
& Pecknold incremental orthotropy, Kupfer biaxial envelope, smeared
rotating cracks with fracture-energy regularisation, and a smeared
Giuffrè-Menegotto-Pinto steel law) into a transparent, reusable MATLAB
codebase. The primary contribution is the implementation itself: a
documented, modular environment that can be inspected, reused, and
extended for further studies.

---

## Features

- **Concrete constitutive model:** incremental orthotropic with continuous
  update of principal stress directions (Darwin & Pecknold, 1977).
- **Equivalent uniaxial strain** as an accumulated history variable
  governing tangent-modulus degradation.
- **Kupfer biaxial failure envelope** (1969) for the extraction of peak
  stress / strain in each principal direction.
- **Smeared rotating crack model** with fracture-energy regularisation
  (Eq. 21 of the paper) to reduce mesh-size dependence.
- **Smeared steel reinforcement** with the uniaxial Giuffrè-Menegotto-Pinto
  monotonic envelope.
- **9-node isoparametric quadrilateral membrane element** (Q9) with
  full 3 × 3 Gauss integration.
- **Newton-Raphson incremental-iterative solver** with sparse global
  stiffness assembly and direct backslash factorisation.

---

## Requirements

- MATLAB R2018a or later
- No additional toolboxes required (uses only core MATLAB functions)

---

## Quick start

Clone the repository and add `src/` to the MATLAB path:

```matlab
addpath(genpath('src'));
cd examples/example1_bresler_scordelis;
inputExemploITrechos;
```

Each example script under `examples/` is self-contained: it builds the
mesh, defines materials, applies loads, calls the solver and produces
the load-displacement curve plus a principal-stress visualisation.

---

## Repository structure

```
rc-nonlinear-fem/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── src/
│   ├── solver/
│   │   └── AnaliseElemFinitosNaoLinear.m   # Newton-Raphson driver
│   ├── element/
│   │   └── ElemMembranaNaoLinearFisicaArmDistr.m   # Q9 membrane element
│   ├── constitutive/
│   │   ├── ParametrosMatrizconstitutivaMembranaNLF.m  # concrete model
│   │   ├── ParametrosAco.m                  # Giuffrè-Menegotto-Pinto steel
│   │   └── ParametrosFissuracao.m           # fracture-energy criterion
│   ├── preprocessing/
│   │   └── Obtemtaxaarmadura.m              # smeared reinforcement ratio
│   ├── postprocessing/
│   │   ├── DesenhaMalhaIndef.m              # undeformed mesh plot
│   │   ├── DesenhaMalhaIndefQ9.m            # undeformed Q9 mesh plot
│   │   ├── DesenhaTensoesprincipais.m       # principal stress + crack map
│   │   ├── ElementMembranaDeformada.m       # deformed configuration
│   │   └── criaSeta.m                       # arrow drawing helper
│   └── utils/
│       ├── obtemgdlelem.m                   # element DOF map
│       ├── PtsGauss1d.m                     # 1D Gauss-Legendre rule
│       └── PtsGauss2d.m                     # 2D Gauss-Legendre rule
├── examples/
│   ├── example1_bresler_scordelis/          # Bresler & Scordelis beam (140 kN)
│   ├── example2_aurich/                     # Aurich beam (50 kN)
│   └── example3_kang/                       # Kang beam (311.36 kN)
└── docs/
    └── theory.md                            # short theory summary
```

---

## Reproducing the paper results

The three benchmark beams used in the paper are reproduced by the input
scripts in `examples/`:

| Example | Reference            | Total load | File |
|--------:|----------------------|-----------:|------|
| I       | Bresler & Scordelis (1963) | 140.00 kN  | `inputExemploITrechos.m` |
| II      | Aurich & Campos Filho (2003) |  50.00 kN  | `inputExemploIIQ9Trechos.m` |
| III     | Kang (1977)          | 311.36 kN  | `inputExemploIIITrechos.m` |

Each script writes to MATLAB figures the load-displacement curve at the
mid-span monitored node and the principal-stress field on the deformed
mesh.

---

## Naming conventions

Function and struct names in the source code follow the original
Portuguese terminology used during the model development. The most
relevant names are listed below for non-Portuguese readers:

| Code name      | Meaning                                         |
|----------------|-------------------------------------------------|
| `Dados`        | global problem data struct                      |
| `Dados.conect` | element connectivity table                      |
| `Dados.coords` | nodal coordinates                               |
| `Dados.restr`  | restraint flags (1 = restrained, 0 = free)      |
| `Dados.Arm`    | reinforcement (smeared steel) parameters        |
| `Dados.ele`    | per-element properties (E0, nu, thickness, ...) |
| `ResConc`      | concrete state (eiu, Sp, E, crack flags)        |
| `ResAco`       | steel state (tangent modulus, stress)           |
| `DadosConc`    | per-element concrete history (Sp, E, eiu)       |
| `Resultado`    | output struct (Sigma, Epslon, TensoesPrinc)     |
| `eiu`          | equivalent uniaxial strain                      |
| `Sp`           | principal stresses                              |
| `Fissura`      | crack flag (0 = intact, 1 = active, 2 = saturated) |
| `aco`          | steel                                           |
| `armadura`     | reinforcement                                   |
| `tensao`       | stress                                          |
| `deformacao`   | strain                                          |
| `passo`        | load step                                       |

All comments, headers, and user-facing messages are in English.

---

## Scope and limitations

The implementation has a deliberately restricted scope, consistent with
the paper:

- **Monotonic loading only.** No cyclic unload/reload rules are included.
- **Smeared (distributed) reinforcement only.** Discrete-bar
  representation with bond-slip is reserved for future extensions.
- **Validation on simply supported beams.** Hyperstatic benchmarks are
  identified as a priority for further work.
- **Newton-Raphson convergence is not strictly quadratic** due to the
  non-smooth envelope-based concrete update.

These items are explicitly discussed in the paper's conclusion and form
the natural roadmap for future development of this codebase.

---

## Citing this work

If you use this code in academic work, please cite both the paper and
the software archive:

```bibtex
@article{silva2026orthotropic,
  author  = {Silva, Edilson and Lobo, Leon and Barreto, Diego and de Souza, Remo},
  title   = {Formulation and implementation of an orthotropic finite element
             for material nonlinear analysis of reinforced concrete
             structures considering cracking under plane stress state},
  journal = {Computers and Concrete},
  year    = {2026},
  note    = {Under review}
}
```

See `CITATION.cff` for the full citation metadata.

---

## License

This code is released under the MIT License — see [`LICENSE`](LICENSE).

---

## Authors

- **Edilson Silva** — Federal University of Pará (UFPA) — `edilson_morais@ufpa.br`
- **Leon Lobo** — University of São Paulo (USP)
- **Diego Barreto** — Federal University of Pará (UFPA)
- **Remo de Souza** — Federal University of Pará (UFPA)

The original implementation was developed at the Graduate Program in
Civil Engineering, UFPA, with contributions from members of the
*Núcleo de Instrumentação e Computação Aplicada à Engenharia* (NICAE).

---

## Acknowledgements

The authors gratefully acknowledge financial support from FAPESPA
(Fundação Amazônia de Amparo a Estudos e Pesquisas).
