# GUIDE — How to Run, Understand, and Demo This Project

*Read this once and everything else in the folder will make sense.*

---

## 0. The 30-second version

This project answers one question:

> **When the kirigami skin of a morphing wing deploys, where does the shear
> center go — and why does that matter?**

The **shear center** is the point where you can push on a beam's
cross-section sideways without twisting it. Push anywhere else and the beam
twists as well as bends. For a normal aircraft wing you want the load line to
pass through the shear center, so you need to know where it is — and in a
*morphing* wing it **moves** as the structure changes shape. That movement is
called **shear-center migration**, and it is exactly what this code
simulates, angle by angle.

To do that, the code must compute, for every deploy angle θ:

1. the shape of the cross-section (parametric geometry),
2. its area, centroid and inertias (section properties),
3. how shear force flows through the walls (shear flow),
4. where the shear center sits (zero-twist condition),
5. how stiff the section is against twisting (Bredt torsion),
6. and therefore how much unwanted twist a given load causes (coupling).

Everything in this folder — MATLAB files, the Python verifier, the browser
demo, this guide — exists to do those six things transparently and to *prove*
that the answers are right.

---

## 1. Run the interactive demo (no MATLAB needed)

Open this file in any modern browser (double-click it, or drag it into a
browser window):

```
demo/kirigami_demo.html
```

That's it. It is one self-contained file: the full solver has been ported
1:1 from the MATLAB code to JavaScript, so **the browser demo and the MATLAB
code produce the same numbers** (the Self-test tab proves it — see §7).

If you have MATLAB instead: run `main` for the full report pipeline, or
`kirigami_app` for the MATLAB UI. The two UIs show the same physics.

---

## 2. Tour of the app, panel by panel

### Left panel — INPUTS (drag & watch)

| slider | what it does |
|---|---|
| **θ (morphing angle)** | 0° = flat closed box, 90° = fully deployed. The top skin rises at each seam by Δ = (c/N)·sin θ. This is *the* morphing knob. |
| **N (cells)** | number of box cells along the chord. More cells = narrower cells, smaller seam lift per cell. |
| **γ (seam position)** | where inside each cell the peak sits. γ = 0.5 would make the deployed section symmetric (nothing interesting happens); the default 0.40 makes it asymmetric — that is what creates Ixy ≠ 0 and the migration. |
| **c, h** | chord and height of the box (metres). |
| **t_up, t_lo, t_web** | skin and web thicknesses (mm). |
| **V_y** | the transverse shear force applied to the section (N). |
| **x_load** | where that force acts along the chord — the load line. |

Every slider re-runs the whole solver instantly (a few ms — it's tiny
matrices).

### Middle panel — CROSS-SECTION (live)

The drawing is the actual cross-section at the current θ:

- **wall colour = local shear-flow magnitude q(s)**: blue ≈ 0 → red = q_max.
  Watch the flow redistribute as you drag θ.
- **+ grey cross** = centroid (the geometric "middle").
- **◆ cyan diamond** = the shear center. **Watch it slide as you deploy.**
- **→ orange arrow** = the applied load line at x_load. Its horizontal gap
  to the ◆ is the eccentricity *e* that produces torque.

### Right panel — OUTPUTS

Every number is computed live at the current θ. Click any **Ⓘ** for a
one-sentence plain-words explanation of that quantity. At the bottom there is
a live **equilibrium check**: the solver adds up the force carried by all the
wall flows and compares it to the applied (Vx, Vy). It should always read
`✔ … residual ~1e-13 N` — if it ever doesn't, something is broken.

### Bottom tabs

| tab | what you see |
|---|---|
| **📈 SC migration** | x_SC(θ) vs centroid x̄(θ) over the full sweep — the two curves separating IS the project result. Second chart: torque T = V·e and twist θ_t = T/K_t. |
| **📋 Results table** | the full θ-sweep table (same columns as the MATLAB table / `kirigami_results.csv`), with a CSV download button. |
| **✅ Self-test** | re-runs 22 analytic checks against textbook closed forms, in your browser, right now. All green = the demo's solver is verified. |
| **📖 Guide** | a short version of this document inside the app. |

### ▶ Deploy animation & 🎬 Guided tour (top bar)

- **Deploy animation** sweeps θ from 0° to 90° and back so you can *see* the
  section breathe and the ◆ migrate.
- **Guided tour** is a 10-step popup walkthrough that points at the numbers
  it talks about. Start there if this is your first open.

---

## 3. The physics, in plain words

**Shear flow = water in pipes.** When you shear the section with a vertical
force V_y, the force has to travel from the load point into the walls, like
water finding its way through a pipe network. The wall thickness is the pipe
diameter: thin walls carry less, so the flow *speed* (N/m) is higher there.
The solver computes that distribution q(s) along every wall.

**Centroid ≠ shear center.** The centroid is the balance point of the area.
The shear center is the point where a sideways push causes bending but *no
twist*. They coincide only for doubly-symmetric sections. For this deployed
kirigami box they differ — at θ = 45° the centroid sits at ȳ = 57.4 mm while
the shear center has dropped to y_SC = 28.2 mm, nearly 30 mm below it.

**Why does the shear center move?** Two mechanisms, both visible in the app:

1. The seam lifts the top skin, so the *area distribution* changes — the
   centroid itself rises (ȳ: 40 mm → 67.6 mm from 0° to 90°).
2. The deployed section is asymmetric fore-aft, so the inertia matrix gains
   an off-diagonal term **Ixy**. Ixy is what "mixes" horizontal and vertical
   bending; it also tilts the shear-flow pattern so its resultant no longer
   passes through the centroid. The x-shear-center moves out to
   x_SC = 151.4 mm (vs centroid 148.6 mm) at 90°.
   **Try it:** set γ = 0.50 → Ixy → 0, the two SC curves flatten onto the
   centroid, migration stops. Set γ back to 0.40.

**The twist you care about.** The load line sits at x_load (default 300 mm =
trailing edge). The horizontal gap between the load line and the shear
center is the eccentricity *e*, and it produces a torque **T = V·e**. The
section twists by θ_t = T/K_t, where K_t is the torsional stiffness from the
multi-cell Bredt solve. That bend-twist coupling is the engineering
consequence of the migration — in a real wing it interacts with aeroelastic
divergence/flutter margins, which is why you must know where S goes.

**The "cells" trick (why the method is clever).** A closed box is
statically indeterminate: infinitely many shear-flow patterns carry the same
total force. The solver cuts (slits) the box open — then the flow is unique
and follows from simple equilibrium along a tree of walls — and finally adds
one constant circulation per cell, chosen so each cell experiences **zero
twist** (compatibility ∮q ds/t = 0). Cutting differently (left vs right)
gives the identical answer — the Self-test proves it. This is the classical
Multi-cell analysis, done with nothing but linear algebra you can read.

---

## 4. A 3-minute demo script (use this when presenting)

1. **Open the demo, θ = 45°.** "This is a kirigami wingbox cross-section.
   The colours are the shear force flowing through the walls; the diamond is
   the shear center, the cross is the centroid. The project question: where
   does the diamond go as the wing morphs?"
2. **Drag θ slowly to 90°.** "The skin deploys, the section gets deeper and
   asymmetric — and the shear center drops ~30 mm below the centroid and
   drifts aft. That drift is the migration we quantify."
3. **Drag γ to 0.50.** "If I make the mechanism symmetric, Ixy vanishes and
   the migration disappears — confirming *why* it happens. Back to 0.40."
4. **Click the 📈 SC migration tab.** "Over the full sweep the SC moves
   smoothly; the gap to the centroid is the eccentricity that turns the
   lift into a torque." Point at the second chart: "That torque twists the
   section by up to a twentieth of a degree per metre with our loads —
   small, but it couples bending into torsion."
5. **Click the ✅ Self-test tab.** "And all 22 checks against textbook
   closed forms pass in your browser right now — same solver, same answers
   as the MATLAB code."

---

## 5. Which MATLAB file does what

| MATLAB file | what it does (plain words) | where you see it in the demo |
|---|---|---|
| `default_params.m` | every number (sizes, thicknesses, loads) in one place | the sliders' default positions |
| `kirigami_geometry.m` | builds the wall layout for angle θ, slits the top skin | the drawing's shape |
| `section_properties.m` | A, centroid, Ix, Iy, Ixy by exact strip integrals | the first output rows |
| `shear_flow.m` | the q(s) "water in pipes" solve (slit → tree → cells) | the wall colours + q_max |
| `shear_center.m` | two unit-load cases + zero-twist condition → S | the ◆ diamond |
| `bredt_torsion.m` | multi-cell torque solve → J, K_t | K_t row |
| `morphing_simulation.m` | sweeps θ, assembles the results table | the Results table tab |
| `visualization.m` | the report figures | the two chart tabs |
| `kirigami_animation.m` | the morphing animation | ▶ Deploy animation |
| `kirigami_app.m` | the MATLAB interactive UI | the whole demo (browser twin) |
| `validation_ansys.m` | ANSYS export + comparison | (Part 10, offline) |
| `self_test.m` | the 22-check validation suite | the ✅ Self-test tab |

The math behind each, with derivations and references, is in
`MODEL_MATH.md` (§4 properties, §5 shear flow, §6 shear center, §7 torsion).

---

## 6. How we know the answers are right

Three independent implementations of the same solver:

1. **MATLAB** (`shear_flow.m` + `self_test.m`) — the deliverable.
2. **Python** (`verify_core.py`) — an independent mirror; `py verify_core.py`
   re-runs all checks and prints PASS lines.
3. **JavaScript** (inside `demo/kirigami_demo.html`) — the Self-test tab.

All three now agree, and all pass the same 22 checks to ≤ 1e-9 (most to
machine precision ~1e-15): closed-form properties at θ = 0, equilibrium for
**both** load directions at θ = 0° and 60°, symmetric-box SC = centroid,
textbook C-channel shear center, single-cell Bredt formulas, and
cut-independence.

**Honesty note (worth repeating in your report):** during the verification
for this demo, an inconsistency was discovered in the original
unsymmetric-bending shear-flow coefficients: the Ixy cross-terms had the
wrong sign (`+Vx·Ixy`, `+Vy·Ixy` instead of `−Vx·Ixy`, `−Vy·Ixy`). The old
version passed its own tests only because those tests applied a pure-Vy load
at θ = 0°, where Ixy = 0 hides the error. The coefficients were corrected in
all three implementations, the test suite was extended to check equilibrium
for pure-Vx loads too, and an independent "truth-solve" (solve for the
coefficients whose flow resultant equals the applied force) confirmed the
fix to machine precision. MODEL_MATH.md §5.1 documents the correct formula.
Moral: equilibrium checks must exercise the terms they claim to verify.

---

## 7. FAQ

**Q: Why does the twist chart show such small angles?**
A: K_t for a 4-cell closed box is large (~2.5×10⁵ N·m/rad), and the demo
load is 1 kN. The *coupling* is the point, not the size of the twist. Raise
V_y or move x_load to see it grow; also note K_t itself *drops* as the box
deploys (its cells open up), which amplifies the effect at high θ.

**Q: The shear center jumps when I change N. Why?**
A: N changes the cell widths, the seam lift Δ = (c/N)·sin θ and the number
of walls — the whole section changes. That's the parametric study, not a bug.

**Q: Where's the warping?**
A: Deliberately excluded (uniform Saint-Venant torsion only) and stated as a
limitation in MODEL_MATH.md §2.4. The open top skin makes warping
significant at large θ — a natural "future work" slide.

**Q: Is the browser demo really the same math?**
A: Yes — a line-by-line port of the five MATLAB solver functions; the
Self-test tab re-proves it every time you open the page.

**Q: What do I hand in?**
A: `main` produces the console calculation dump, Figures 1–11, the results
table (`kirigami_results.csv`) and the animation; `self_test` output is your
validation section; MODEL_MATH.md is the theory chapter skeleton; the demo
link is your live-showcase slide.
