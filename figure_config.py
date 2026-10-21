"""Appearance and numerical controls; model parameters live in scenarios.json."""

STYLE = {
    # These are FINAL printed dimensions, not a large image to shrink in TeX.
    "figure_width_mm": 85.0, "figure_height_mm": 25.5, "dpi": 600,
    "font_size": 8.0, "label_size": 8.0,
    "panel_letter_size": 7.0,
    # Position in axes fractions: x from the left, y from the bottom.
    "panel_letter_position": (0.02, 0.98),
    "left_margin_mm": 0.8, "right_margin_mm": 0.9,
    "panel_gap_mm": 1.7, "plot_bottom_mm": 0.5,
    # Optical spacing correction after computing the common disk scale.
    # Negative values move a whole panel left without changing its geometry.
    "panel_x_offsets_mm": {"b": 0.20, "c": -0.45},
    "black": "#000000", "km": "#7A5195", "boundary": "#3B5D9A",
    "initial_contour": "#A32A31", "contour": "#DDB8BB",
    "exterior_contour": "#E6C9CB",
    "flow_width": .65, "km_width": .65, "boundary_width": .60,
    "contour_width": .50, "initial_contour_width": .55,
    "marker_size": 2.5, "marker_edge_width": .45,
    "show_iteration_labels": False, "iteration_label_size": 4.5,
    "iteration_label_halo_width": .5,
    "state_dot_size": 9.0, "mle_star_size": 25.0,
    "marker_spacing_points": 4.0, "marker_endpoint_clearance_points": 3.0,
    "marker_arrow_clearance_points": .25,
    "marker_avoid_arrows_in": ("b",),
    "label_halo_width": .8,
    "arrow_fraction": .64,
    # In (a), keep the arrow before the bend where the two curves converge.
    "first_arrow_fraction_by_panel": {"a": .40, "b": .46, "c": .83},
    # Segment indices follow boundary crossings (0 = initial interior arc).
    # Panel (c)'s optional exterior dart lies beyond the sharpest bend.
    "arrow_fraction_by_segment": {"b": {1: .69}, "c": {1: .80}},
    # Optional concave dart: length, full width and notch
    # distance behind the tip, all in typographical points.
    "arrow_head_length_pt": 3.3,
    "arrow_head_width_pt": 3.1,
    "arrow_head_notch_pt": 2.1,
    # All heads use the same dimensions; individual overrides remain available.
    "first_arrow_head_scale_by_panel": {},
    # A narrow white outline separates (a)'s black dart from the adjacent
    # purple curve without changing the head geometry or either trajectory.
    "arrow_head_halo_pt_by_panel": {"a": .45},
    "z_range": (-1.12, 1.10),
    "show_initial_contour": True, "show_contours": True,
    "show_exterior_contours_in": ("b",), "show_km_markers": False,
    "show_flow_arrows": False, "contour_resolution": 850,
}

# Panel-specific limits remove unused horizontal space, while the renderer
# keeps an identical millimetres-per-Bloch-unit scale in all three panels.
PANEL_X_RANGES = {
    "a": (-1.08, 1.04),
    "b": (-1.08, 1.55),
    "c": (-1.08, 1.45),
}

# Optional small numbers give the actual update index k (initial state k=0).
# Iteration circles and numbers are disabled in the manuscript figure.
# Placement only: indices are always taken from the computed KM trace.
# Each entry is (offset in points, horizontal anchor, vertical anchor).
ITERATION_LABELS = {
    "a": {"default": ((2, 0), "left", "center"),
          16: ((2, 1), "left", "center")},
    "b": {"default": ((0, -2), "center", "top")},
    "c": {"default": ((2.3, -2), "left", "center"),
          1: ((-2.2, 0), "right", "center")},
}

# (text, offset in typographical points, horizontal anchor, vertical anchor).
POINT_LABELS = {
    "a": {
        "initial": (r"$\sigma_0$", (0, -2), "center", "top"),
        "mle": (r"$\hat{\sigma}_{\mathrm{MLE}}$", (0, 3 - .25*72/25.4), "center", "bottom"),
    },
    "b": {
        "initial": (r"$\sigma_0$", (0, -2), "center", "top"),
        "mle": (r"$\hat{\sigma}_{\mathrm{MLE}}$", (-4.7, 3 - .25*72/25.4), "center", "bottom"),
        "hit": (r"$\sigma_{\mathrm{hit}}$", (2, 0), "left", "top"),
        "outside": (r"$\sigma_{\mathrm{out}}$", (0, 1.3), "center", "bottom"),
        "km": (r"$\sigma^{\mathrm{KM}}_{\mathrm{lim}}$", (-4, -.5), "right", "center"),
    },
    "c": {
        "initial": (r"$\sigma_0$", (2, 0), "left", "center"),
        "mle": (r"$\hat{\sigma}_{\mathrm{MLE}}$", (0, 3 - .25*72/25.4), "center", "bottom"),
        "hit": (r"$\sigma_{\mathrm{hit}}$", (2, 0), "left", "center"),
        "km": (r"$\sigma^{\mathrm{KM}}_{\mathrm{lim}}$", (-1.8, -2.0), "right", "bottom"),
    },
}

NUMERICS = {
    "flow_samples": 2400, "flow_rtol": 2e-11, "flow_atol": 2e-13,
    "flow_endpoint_tolerance": 1e-10,
    "km_decimal_digits": 80, "km_step_tolerance": "1e-26",
    "km_max_iterations": 10000,
}
