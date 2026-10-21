# Euclidean gradient flow versus KM: three Bloch-plane examples

This standalone package produces one row of three state-space panels, showing a comparison between the Euclidean (Hilbert–Schmidt) gradient flow and the physical KM iteration on the same likelihoods and from the same initial states. Both methods are shown in all three panels. The compact layout targets one column at 85 × 25.5 mm, with 8 pt state labels and 7 pt panel letters `(a)`, `(b)`, and `(c)` in the upper left. Since the manuscript caption supplies the scenario names and visual key, there are no descriptive headings, legend, axes, ticks, or coordinate labels. The horizontal coordinate is $x$ and the vertical coordinate is $z$ in the Bloch-plane definition below. The blue circle marks $x^2+z^2=1$. All three disks have the same physical diameter and scale. Panel-specific horizontal offsets control the spacing independently of disk sizes and point-label placement.

## View or run the figure

The supplied `outputs` folder contains `euclidean_km_bloch.png`, `euclidean_km_bloch.pdf`, and `euclidean_km_bloch.svg`. Open the PNG in VS Code for a preview, use the PDF for a manuscript, or insert the SVG into PowerPoint. PNG export defaults to **600 dpi**; PDF and SVG use vector graphics and remain sharp when enlarged. Bundled STIX fonts avoid dependence on locally installed fonts or TeX.

For a single-column REVTeX figure, include the vector PDF with `\includegraphics[width=\columnwidth]{euclidean_km_bloch.pdf}`. Its native width is 85 mm (3.34646 in). The compact spacing removes unnecessary exterior margins without changing the geometry. Review the PDF at its final printed size: extra PNG resolution does not make physically small lettering larger.

Use **Python 3.11 or newer**; the project was tested with Python 3.12. On macOS:

1. Extract the package. In VS Code, choose **File → Open Folder…** and select
   `euclidean_km_bloch`. Use **File → New Window** first to keep another project open.
2. Choose **Terminal → New Terminal** and run:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python make_figure.py
```

If necessary, replace `python3` in the first command with a newer installed interpreter such as `python3.12`. Select `.venv/bin/python` with VS Code's **Python: Select Interpreter** command. In subsequent terminals, activate the existing environment with `source .venv/bin/activate`.

`python make_figure.py` recomputes the trajectories, updates the cache and exports all three formats. After changing only styling, reuse the cached data:

```bash
python make_figure.py --plot-only
python make_figure.py --plot-only --dpi 900
```

To export a variant into a separate folder:

```bash
python make_figure.py --plot-only --output-dir alternate_outputs
```

The cached `data/scenario_a.npz`, `data/scenario_b.npz`, `data/scenario_c.npz`, and `data/diagnostics.json` record the computed curves and numerical information. Cache fingerprints guard the model and numerical settings: `--plot-only` refuses an
incompatible cache. After changing physics or numerical parameters, rerun `python make_figure.py` without `--plot-only`.

## What the figure shows

| Element | Meaning |
|---|---|
| Blue circle | Physical-state boundary in the displayed Bloch plane |
| Solid black curve | Continuous Euclidean gradient flow inside the physical disk |
| Dashed black curve | The same unconstrained flow outside the physical disk |
| Violet curve | All computed discrete KM iterates, joined by straight visual-guide segments |
| Dark red dashed curve | Initial likelihood level, ℓ(σ) = ℓ(σ₀), clipped to the physical disk |
| Faint solid red curves | Other likelihood levels; exterior contours in (b) are fainter |
| σ₀ | Initial state shared by both methods within each panel |
| $\hat{\sigma}_{\mathrm{MLE}}$ | Physical maximum-likelihood state |
| σ_hit | First boundary crossing of the Euclidean flow |
| σ_out | External unconstrained maximum in (b) |
| σ_lim^KM | Last computed KM iterate used to estimate an asymptotic boundary limit |

The black flow in (c) becomes **solid again after reentering** the disk. Dashing
therefore indicates an exterior location, not every point after the first
boundary encounter. No outside-maximum → physical-MLE projection arrow is drawn.
The MLE is calculated by constrained optimization, not by assuming a projection
rule.

The dark red contour passes through the initial state and separates higher from
lower likelihoods within the disk. If the corresponding upper level set touches
the physical circle, parts of that circle also contribute to its boundary;
the clipped red contour need not form a closed loop by itself.

The violet KM line connects **every computed discrete iterate** in sequence.
Individual iteration circles and iteration numbers are omitted in all panels.
The solid violet limit dots in (b) and (c) identify
the named limiting states. The plot shows the geometry of the trajectories;
it does not encode the number of iterations through visible markers.
Arrowheads are also omitted from all panels. Each trajectory starts at
\(\sigma_0\) and proceeds toward its labeled endpoint; the black line is a
continuous Euclidean flow, while the purple line connects discrete KM states.
The included runs stop after 915 updates in (a), 123 in (b), and 552 in (c),
using the stated small-step tolerance. In (c), the KM limit is not the MLE.
The joining segments do not describe a continuous gradient flow. A KM-limit dot
is not a claim that a finite iteration hits the boundary: tiny positive
distances from the boundary can be visually indistinguishable from zero.

For these examples:

- **(a)** Both methods approach the same interior MLE along distinct trajectories.
- **(b)** Euclidean flow crosses the physical boundary and approaches the external
  maximum; KM instead approaches a nonoptimal boundary state.
- **(c)** Euclidean flow exits and reenters the disk before approaching the interior
  MLE; KM approaches a nonoptimal boundary state.

These observations do not establish a general convergence or performance claim.
A small KM update alone does not establish likelihood optimality.

The graphic has no legend. In the manuscript caption, describe the physical
boundary, Euclidean flow, KM iteration, exterior continuation, MLE stars, and
contour through \(\sigma_0\), and introduce the three scenarios by panel letter.

## Adjust the figure

Edit `figure_config.py`:

- `STYLE`: colors, fonts, line widths, marker sizes, physical figure dimensions
  and other display settings. The compact layout uses millimetres for physical
  sizing and points for typography; adjust these controls rather than resizing
  the exported figure to fix label spacing.
- `STYLE["panel_letter_size"]` and `panel_letter_position`: size and upper-left
  placement of the three panel letters. Position uses axes fractions, with
  `(0, 1)` at the upper-left corner. Descriptive headings are not drawn.
- `STYLE["panel_x_offsets_mm"]`: optional horizontal shifts of entire panels,
  applied after computing the common scale. The default `{"b": 0.2, "c": -0.45}`
  moves (b) 0.2 mm right and (c) 0.45 mm left relative to the uniformly spaced
  layout. Use `{}` for zero offsets.
- `POINT_LABELS`: annotation text, alignment and placement.
- `STYLE["show_km_markers"]` and `STYLE["show_iteration_labels"]`: both are
  `False` by default. Optional marker and label controls
  are available for diagnostic views; enable both to display numbered
  circles. `ITERATION_LABELS` sets panel and per-index number positions.
- `STYLE["arrow_head_length_pt"]`, `arrow_head_width_pt`, and
  `arrow_head_notch_pt`: dimensions of optional concave dart arrowheads, in points.
  `STYLE["show_flow_arrows"]` is `False` by default; set it to `True`
  to display the Euclidean-flow arrowheads.
  All arrowheads, including the first one in (c), use the same dimensions.
- `STYLE["arrow_head_halo_pt_by_panel"]`: optional white outline width in
  points, used only when arrowheads are enabled. Panel (a)'s optional dart
  has a 0.45 pt outline to distinguish it from the nearby purple curve.
- `STYLE["marker_size"]`: diameter of optional KM circles (2.5 pt), used only
  when markers are enabled. Those circles select actual iterates by visual
  separation, not consecutive update indices. `marker_avoid_arrows_in` and
  `marker_arrow_clearance_points` control optional arrow clearance.
- `STYLE["arrow_fraction_by_segment"]`: optional positions along individual
  Euclidean-flow segments, separated at physical-boundary crossings. Panel
  (a)'s interior arrow sits before the bend, where it is more distinct from
  the KM curve; (b)'s exterior arrow sits at 69% of its segment's arc length;
  (c)'s exterior arrow is placed after the sharpest part of the bend.
- `NUMERICS`: ODE and KM precision, tolerance and iteration settings.

`scenarios.json` defines the measurement models and initial states. Appearance
changes leave the numerical trajectories unchanged; changes to `NUMERICS` or
the scenarios require recomputation. The figure uses equal aspect so the disk
and trajectories are not stretched.

## Mathematical and numerical construction

The invariant plane is
\(\rho(x,z)=(I+x\sigma_x+z\sigma_z)/2\), physical when \(x^2+z^2\le1\).
Two rotated Pauli measurement axes lie in this plane; a third perpendicular
setting has balanced outcomes. For positive setting weights summing to one,

\[
\Pi_{j,\pm}=\frac{w_j}{2}(I\pm\mathbf n_j\cdot\boldsymbol\sigma),
\qquad f_{j,\pm}=\frac{w_j}{2}(1\pm m_j).
\]

Writing \(u\) in the rotated coordinates, Euclidean gradient ascent obeys
\(\dot u_j=w_j(m_j-u_j)/(1-u_j^2)\). Its mathematical continuation outside the
physical disk remains in the domain of positive measurement probabilities;
it does not represent a physical density matrix there.

With \(R(\rho)=\sum_i[f_i/\operatorname{Tr}(\Pi_i\rho)]\Pi_i\), the KM update is
\(\rho_+=\sqrt\rho R(\rho)\sqrt\rho\). It is normalized because
\(\operatorname{Tr}(\rho R)=1\), and it coincides with the Petz update **for these
measurement channels**. That equivalence is not asserted for general quantum
channels.

The numerical defaults are:

- Euclidean flow: SciPy Radau, relative tolerance `2e-11`, absolute tolerance
  `2e-13`, convergence distance `1e-10`, with boundary-event detection.
- KM: 80 Decimal digits plus 10 guard digits, step-norm tolerance `1e-26`,
  maximum 10,000 updates. A determinant recurrence tracks tiny boundary gaps;
  no projection, eigenvalue clipping or artificial ridge is applied.

`likelihood.py` implements the model and Euclidean flow; `km_iteration.py`
implements the physical iteration. `make_figure.py` is the entry point.
`layout_helpers.py` places the three panel letters.
