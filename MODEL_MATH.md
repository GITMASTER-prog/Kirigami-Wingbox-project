# Mathematical Model — Morphable Cellular Kirigami Wingbox

**Project:** Parametric Simulation of Shear-Center Migration and Bending–Torsion
Coupling in a Morphable Cellular Kirigami Wingbox
**Units:** SI everywhere (m, N, Pa, N/m, N·m, rad).
**Companion code:** `main.m` plus the functions listed in `README.md`.
Every equation below is implemented exactly once in the MATLAB code file noted.

---

## 1. Physical model and motivation

The model is a **simplified cellular kirigami wingbox cross-section** inspired by
2023 papers on cellular kirigami morphing aerostructures (e.g. *“Morphing
aerostructures based on cellular kirigami”*). Those papers use a two-dimensional
lattice of cells whose walls rotate about discrete hinges so that the skin can
change shape; the driving kinematics is a **fold/deploy angle θ**.

We idealise one wingbox cross-section as:

- a flat **bottom skin** along the full chord,
- **vertical webs** at each cell boundary,
- a **top skin** made, in each cell, of **two flat panels hinged at the web
  tops and joined at a raised seam** (a “peak”).

The seam height is the single morphing degree of freedom:

```
Δ(θ) = (c/N) · sin(θ)
```

- θ = 0° → Δ = 0 → the top skin is flat → **closed rectangular N-cell box**
- θ = 90° → Δ = dc → **fully deployed**, deepest (and most asymmetric) section

The seam is offset from the cell centre (`γ = P.seam_pos = 0.40`, not 0.5), so
the deployed sections are chordwise-asymmetric. This is what produces
Ixy ≠ 0, x_SC ≠ x_centroid, and a genuine shear-center **migration** — exactly
the physics the project wants to demonstrate.

This is a *course-project simplification* of the kirigami mechanism: the
real rotating-flap kinematics is reduced to an equivalent seam-lift. No
research-paper equation is reproduced; only classical thin-walled-beam theory
is used (Sections 4–7).

## 2. Assumptions (explicit list)

1. **Thin-walled idealisation** — every wall is a straight strip, t ≪ L;
   thickness is accounted for in the wall's axial direction only.
   Terms of order (t/L)² are dropped.
2. **Linear elasticity**, small deflections; E and G constant.
3. **Euler–Bernoulli/Vlasov beam behaviour**: cross-sections do not deform in
   their own plane (rigid contour); shear-lag and wall buckling are ignored.
4. **Uniform torsion** (Saint-Venant) only — warping torsion is **not**
   included in Kt (stated limitation; the open top skin makes warping
   significant at large θ).
5. Shear flows are computed from **transverse forces applied at the
   centroid**; the shear center follows from the zero-twist condition
   (Section 6). Superposition holds (linear theory).
6. Cells are **prismatic** along the span; the sweep varies the section
   geometry, not the spanwise load.
7. Thermal effects, residual stresses, joints/hinge flexibility: ignored.

## 3. Parameters (all in `default_params.m`)

### Geometric (m)
| symbol | field | meaning |
|---|---|---|
| c | `P.chord` | wingbox chord |
| h | `P.height` | web height (wingbox depth) |
| N | `P.n_cells` | number of cells |
| dc = c/N | — | cell width |
| γ | `P.seam_pos` | seam offset fraction within a cell (0 < γ < 1) |
| t_up | `P.t_up` | upper skin thickness |
| t_lo | `P.t_lo` | lower skin thickness |
| t_w | `P.t_web` | web (side-wall) thickness |

### Morphing sweep
`theta_min = 0°, theta_step = 15°, theta_max = 90°` — the deployed geometry
is set by Δ(θ) above.

### Material (Pa)
E = 72×10⁹ (2024-T3), G = 27×10⁹, ν = 0.33 (informational).

### Loading
V_y (default 1000 N), V_x (default 0), load line x = `P.load_pt(1)`,
reference beam length L_ref = 0.5 m for twist estimates.

## 4. Cross-section properties — `section_properties.m`

Wall midline data per segment e: endpoints (x1,y1)-(x2,y2), thickness t_e,
length L_e = √(dx²+dy²), midpoint (xm,ym).

**Area:**
```
A = Σ_e  t_e · L_e
```

**Centroid** (first moments Qx, Qy about the axes, then divide):
```
Qx = Σ t_e L_e ym ,   Qy = Σ t_e L_e xm
x̄ = Qy/A ,  ȳ = Qx/A
```

**Second moments and product of inertia about the centroid** — for a thin
strip, the exact integral splits into a parallel-axis part plus the strip's
own (b h³/12) part; for inclined strips that self-part is (t_e L_e)·(d²/12)
where d is the direction vector component:
```
Ix  = Σ t_e L_e [ (ym − ȳ)²  + dy²/12 ]
Iy  = Σ t_e L_e [ (xm − x̄)²  + dx²/12 ]
Ixy = Σ t_e L_e [ (xm − x̄)(ym − ȳ) + dx·dy/12 ]
```
These are **exact** for thin rectangles (verified to machine precision in
`self_test.m` Check 1) and reduce to the familiar parallel-axis theorem when
walls are axis-aligned.

**Validity check** (also in self_test): at θ = 0, A = c(t_up+t_lo) + (N+1)h t_w,
ȳ = h/2, Ix = 2c·t_up(h/2)² + (N+1)t_w h³/12, Ixy = 0.

## 5. Shear flow — `shear_flow.m`

### 5.1 Governing equation (unsymmetric bending + shear)

For an arbitrary (possibly unsymmetric) section under Vx, Vy the
axial-stress equilibrium of a wall element gives the classic result
(e.g. Megson, *Aircraft Structures*, §16–20; Budynas, *Shigley*, Ch. 8):

```
q(s) = a·Qx(s) − b·Qy(s)          (open section)
Δ0 = Ix·Iy − Ixy²
a = −(Vy·Iy − Vx·Ixy)/Δ0
b = +(Vx·Ix − Vy·Ixy)/Δ0
```

with the **centroidal** first moments of the cut portion
```
Qx(s) = ∫₀ˢ (y−ȳ) t ds ,    Qy(s) = ∫₀ˢ (x−x̄) t ds
```
When Ixy = 0 and Vx = 0 this reduces to the familiar q = −Vy·Qx/Ix.

**Sign convention (important):** q > 0 means flow in the direction of
increasing s *of the stored node direction* n1→n2; the constants a,b above
carry the sign so that the **resultant of the flow equals the applied
force**. Two pitfalls, both caught by `self_test.m` Check 2 (which tests
BOTH load directions at θ = 0° and 60°):

1. the naive a = +(Vy·Iy−…)/Δ0 gives the *reaction* instead of the applied
   force;
2. the **Ixy cross-terms carry MINUS signs** — an earlier draft had
   +Vx·Ixy and +Vy·Ixy, which satisfies equilibrium only when Ixy = 0 and
   silently breaks at deployed θ. This was found by an independent
   truth-solve (solve for the coefficients whose flow resultant equals the
   applied force) and confirmed against the exact unsymmetric-bending
   derivation; the same fix is applied in `shear_flow.m`, `verify_core.py`
   and the interactive demo.

### 5.2 Numerical construction (transparent, no black boxes)

1. **Slit the section.** Kirigami geometry stores the top skin with a
   duplicated node at the cut (`slit_side='right'|'left'`), making the
   analysis graph a **tree**. The closed-section answer is independent of
   the cut location (self_test Check 6).
2. **Tree walk.** BFS from the free end accumulates, for each segment,
   the accumulated Qx, Qy at its far end. Each segment contributes
   exactly `t·L·(midpoint offset)` because Q varies linearly along it.
3. **Endpoint q and exact integral.** q is **quadratic in s** along a
   straight wall; the code evaluates q at both ends and the integral
   ```
   ∫ q ds = a[Qx(s1)·L + t((y1−ȳ)L²/2 + ty L³/6)]
          − b[Qy(s1)·L + t((x1−x̄)L²/2 + tx L³/6)]
   ```
   in closed form — no quadrature error.
4. **Cell redundancies.** The tree killed N independent circulation
   paths; one constant q0k per cell k restores continuity. Each cell must
   have **zero twist**:
   ```
   ∮_k  q ds/t = 0      (compatibility / no ring shear strain)
   ```
   which, with q = q_open + Σ q0m, is an N×N linear system
   `A·q0 = b` with A the loop-incidence matrix (`±L_e/t_e`), b from the
   open-section integrals.
5. **Superposition.** q_final = q_open + Σ q0k on every wall; the resultant
   (F_x, F_y) = Σ Iq_e·(tx,ty) is checked against the applied load
   (should match to ~1e-12, printed by the verbose mode).

## 6. Shear center — `shear_center.m`

**Definition.** The shear center S is the point where the line of action of
the resultant shear force must pass so that the section twist is zero.
Equivalently: the moment of the internal shear flows about S vanishes.

**Zero-twist condition about any reference point P:**
```
0 = M_P + ((P − S) × R)_z
```
where M_P is the moment of the internal flow about P and R = (Vx,Vy) the
resultant. Taking P = centroid c and S = c + (ex, ey):
```
0 = M_c + (−ex·Vy + ey·Vx)          [z-component of the cross product]
```
**Two independent unit load cases** (linear superposition):
```
pure Vy:  ex = +M_y0 / Vy
pure Vx:  ey = −M_x0 / Vx
```
with M_y0, M_x0 the flow moments about the centroid for each case.

**Moment of the flow.** Along a straight wall the moment arm integrand
(r × t̂)_z is *constant* (the s-dependence cancels — derive it once, it is
two lines), so exactly:
```
M_c = Σ_e Iq_e · [ (x1_e − x̄)·ty_e − (y1_e − ȳ)·tx_e ]
```
with Iq_e = ∫ q ds per segment (exact from Section 5.2).

**Absolute position:** x_SC = x̄ + ex, y_SC = ȳ + ey.
**Do not** confuse S with the centroid: for this kirigami box they differ at
every deployed θ (that difference is the eccentricity the project studies).

**Verification:** for the C-channel (open section, no cells) the code
reproduces the textbook `e = b²h²t/(4Ix)` to 1e-15 (self_test Check 4), and
for the symmetric closed box S = centroid exactly (Check 3).

## 7. Torsion and bending–torsion coupling

### 7.1 Multi-cell Bredt–Batho — `bredt_torsion.m`

For each cell k with enclosed area A0_k (shoelace formula on the CCW loop
polygon) and constant cell flow q0_k:
```
Equilibrium:  T = 2 Σ_k q0_k A0_k
Twist rate:   ∮_k q ds/t = 2 A0_k G (dφ/dz)     (same q as in 5)
```
Solving the combined (N+1)×(N+1) system for unit torque gives the twist
rate φ′; then
```
J  = T / (G·φ′)                  [m⁴]
Kt = G·J / L_ref                 [N·m/rad]   (uniform-torsion stiffness)
```

### 7.2 Eccentricity and coupling — `morphing_simulation.m`

The load acts along the vertical line x = P.load_pt(1). The signed
eccentricity and torque are
```
e(θ) = x_loadline − x_SC(θ)          [m]
T(θ) = Vy · e(θ)                     [N·m]      (CCW positive)
```
and the **estimate** of the torsional response (uniform torsion, no
warping, static):
```
θ_t(θ) = T(θ) / Kt(θ)                [rad]      (reported in deg)
```
This is exactly the project formula T = V·e, θ_t = T/Kt, with all
assumptions from Section 2 applying. When the load line passes through S,
T = 0 — that is the *definition* of the shear center, so the coupling study
is self-consistent with Part 4.

## 8. Validation checklist (what `self_test.m` proves)

| # | check | expected | verified |
|---|---|---|---|
| 1 | θ=0 box A, ȳ, Ix, Ixy vs closed forms | machine ε | 1e-16 |
| 2 | flow resultant = applied force, θ=0°,60°, Vy and Vx | machine ε | 1e-13 |
| 3 | symmetric box → S = centroid | 0 | 1e-17 |
| 4 | C-channel → e = b²h²t/(4Ix) | analytic | 7e-16 |
| 5 | single-cell Bredt J = 4A0²/∮ds/t, q0 = T/2A0 | analytic | 1e-15 |
| 6 | cut-independence (right vs left slit) | identical q | 9e-16 |

(Errors quoted from the development harness; `self_test.m` recomputes them
in MATLAB with tolerance 1e-9…1e-12.)

## 9. Mapping to the 10 project parts

| Part | where |
|---|---|
| 1 parametric geometry | `kirigami_geometry.m` |
| 2 section properties | `section_properties.m` |
| 3 shear flow | `shear_flow.m`, Fig 8 |
| 4 shear center | `shear_center.m`, Figs 3–4, 9 |
| 5 morphing sweep + table | `morphing_simulation.m` (returns MATLAB `table`) |
| 6 visualisation | `visualization.m` (Figs 1–10) |
| 7 coupling | `bredt_torsion.m` + Fig 11 |
| 8 animation | `kirigami_animation.m` |
| 9 UI | `kirigami_app.m` |
| 10 ANSYS validation | `validation_ansys.m` |

## 10. Known simplifications vs the real kirigami paper

- The rotating-flap/sheet kirigami mechanism is reduced to an equivalent
  seam-lift; the morphing amplitude Δ(θ) is a kinematic idealisation, not a
  paper result.
- Joints are ideal hinges (zero radius, zero mass).
- No skin wrinkling, no hinge stiffness, no pre-stress.
- Warping torsion neglected in Kt (Section 2.4).
These are exactly the kind of simplifications expected in a structures
course project — state them in the report.
