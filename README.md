# Kirigami Wingbox — Shear-Center Migration Simulator (MATLAB)

Course project: **Parametric Simulation of Shear-Center Migration and
Bending–Torsion Coupling in a Morphable Cellular Kirigami Wingbox.**

The section is a cellular kirigami wingbox whose top skin lifts at the cell
seams as the morphing angle θ goes 0° (closed) → 90° (deployed). Because the
geometry changes, the centroid, Ix, Iy, Ixy, shear flow, shear center and
torsional stiffness all change — and the shear center **migrates**.

> **New here?** Open `demo/kirigami_demo.html` in a browser (no MATLAB
> needed) and read **[GUIDE.md](GUIDE.md)** — a plain-language walkthrough
> of the physics, the app, and how to present it.

## Quick start

```matlab
>> cd <this folder>
>> main            % self-test + sweep + 11 figures + animation + CSV export
>> self_test       % analytic validation only (run this first on a new machine)
>> kirigami_app    % interactive UI (Part 9)
```

Requirements: MATLAB R2020b+ (uses `table`, `sgtitle`, `exportgraphics`,
`uifigure`, `turbo`). No toolboxes required.

## Files

| file | role |
|---|---|
| `main.m` | orchestrator: self-test → sweep → figures → animation → export |
| `default_params.m` | all geometric/material/loading parameters (SI) |
| `kirigami_geometry.m` | Part 1 — parametric cross-section for angle θ |
| `section_properties.m` | Part 2 — A, centroid, Ix, Iy, Ixy |
| `shear_flow.m` | Part 3 — multi-cell q(s) incl. cell redundancies |
| `shear_center.m` | Part 4 — SC from the zero-twist condition |
| `bredt_torsion.m` | Part 7 — multi-cell Bredt J, Kt |
| `morphing_simulation.m` | Part 5 — θ-sweep + results table |
| `visualization.m` | Part 6 — report Figures 1–10 (+Fig 11 coupling) |
| `kirigami_animation.m` | Part 8 — morphing animation |
| `kirigami_app.m` | Part 9 — interactive UI |
| `validation_ansys.m` | Part 10 — ANSYS export + comparison table |
| `self_test.m` | analytic validation suite |
| `MODEL_MATH.md` | **every equation, assumption, unit and derivation** |
| `verify_core.py` | Python mirror of the solver used during verification |
| `demo/kirigami_demo.html` | **interactive browser demo** — live visuals, guided tour, self-test; no MATLAB required |
| `GUIDE.md` | **plain-language guide** — how to run, understand and demo everything |

## Method in one paragraph

Slit the section → open-section tree walk accumulates centroidal first
moments Qx, Qy → q_open = a·Qx − b·Qy (a,b from unsymmetric-bending shear
theory) → one compatibility equation per cell (`∮q ds/t = 0`) solved for the
N cell redundancies → superposition → moment of the flow about the centroid
→ zero-twist condition gives the SC → Bredt multi-cell solve gives J, Kt →
T = V·e, θ_t = T/Kt close the loop. Full details with derivations:
`MODEL_MATH.md`.

## Workflow for the report

1. `self_test` — quote the PASS lines as your code-validation section.
2. `main` — the console prints a full intermediate-calculations dump for the
   three validation configurations (closed / intermediate / deployed).
3. Figures 1–11 are report-ready; set `P.save_figs = true` in
   `default_params.m` to also export PNGs to `figures/`.
4. The results table `T` (also written to `kirigami_results.csv`) is the
   Part 5 deliverable.
5. For Part 10: run `validation_ansys(P, R)` to export node/wall CSVs of the
   three validation configurations, rebuild them in ANSYS, then fill
   `ansys_validation/sc_ansys_template.csv` and run
   `validation_ansys(P, R, true)` for the error-percentage comparison table.

## UI notes (Part 9)

`kirigami_app` exposes all inputs, the six required buttons, and live
outputs. It calls the *same* solver functions as the batch run, so UI and
script results are identical.

## Troubleshooting

- **"no free edge - slit the section"** — the geometry builder always slits
  the top skin; you should only see this if you modified `kirigami_geometry`.
- **Self-test failures** — do not ignore them; they compare against closed
  forms to 1e-9…1e-12. Check that you did not change units of any default.
- **Animation slow** — reduce `P.anim_frames` (default 90).
- **Older MATLAB** — replace `exportgraphics` with `print -dpng`, and
  `sgtitle` with `suptitle`-style workarounds.

## Verification statement

The solver exists in three independent implementations — MATLAB (`shear_flow.m`),
Python (`verify_core.py`) and JavaScript (the browser demo) — which agree and
all pass the same 22 analytic checks to ≤ 1e-9 (most to ~1e-15): closed-form
section properties, equilibrium for **both** load directions (Vy and Vx) at
θ = 0° and 60°, symmetric-box SC = centroid, the textbook C-channel shear
center, single-cell Bredt formulas, and cut-independence.

During this verification an error was found and fixed in the unsymmetric
bending coefficients of `shear_flow.m`: the Ixy cross-terms carried the wrong
sign (`+Vx·Ixy`, `+Vy·Ixy` instead of `−Vx·Ixy`, `−Vy·Ixy`), which violated
equilibrium whenever Ixy ≠ 0; the original tests (pure-Vy load at θ = 0°, where
Ixy = 0) could not catch it. The corrected coefficients, an independent
truth-solve confirmation, and the extended both-direction checks are documented
in `MODEL_MATH.md` §5.1 and `GUIDE.md` §6.
