"""A single Bloch row; drawing controls never change the computed states."""
from pathlib import Path
import os
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, patheffects
from matplotlib.patches import Circle, Polygon
import numpy as np

from figure_config import STYLE, PANEL_X_RANGES, POINT_LABELS, ITERATION_LABELS
from layout_helpers import add_panel_letter

BASE = Path(__file__).resolve().parent


def fonts():
    for path in (BASE / "assets/fonts").glob("*.ttf"):
        font_manager.fontManager.addfont(str(path))
    plt.rcParams.update({
        "font.family": "STIXGeneral", "mathtext.fontset": "stix",
        "font.size": STYLE["font_size"], "axes.spines.top": False,
        "axes.spines.right": False, "axes.linewidth": .50,
        "pdf.fonttype": 42, "svg.fonttype": "path", "savefig.facecolor": "white",
    })


def label(ax, point, spec, color="black"):
    text, offset, ha, va = spec
    obj = ax.annotate(text, point, xytext=offset, textcoords="offset points",
                      ha=ha, va=va, fontsize=STYLE["label_size"], color=color,
                      zorder=12, annotation_clip=False)
    obj.set_gid("state-label")
    obj.set_path_effects([patheffects.withStroke(linewidth=STYLE["label_halo_width"], foreground="white")])


def contours(ax, case, model):
    x = np.linspace(*PANEL_X_RANGES[case["id"]], STYLE["contour_resolution"])
    z = np.linspace(*STYLE["z_range"], STYLE["contour_resolution"])
    xx, zz = np.meshgrid(x, z)
    values = model.value(np.stack((xx, zz), axis=-1))
    levels = sorted(float(model.value(case["center"])) - np.array(case["contour_drops"]))
    disk = Circle((0, 0), 1, transform=ax.transData)
    if STYLE["show_contours"]:
        if case["id"] in STYLE["show_exterior_contours_in"]:
            ax.contour(xx, zz, values, levels=levels, colors=STYLE["exterior_contour"],
                       linewidths=STYLE["contour_width"], linestyles="solid", zorder=0)
        light = ax.contour(xx, zz, values, levels=levels, colors=STYLE["contour"],
                           linewidths=STYLE["contour_width"], linestyles="solid", zorder=1)
        light.set_clip_path(disk)
    if STYLE["show_initial_contour"]:
        initial = ax.contour(xx, zz, values, levels=[float(model.value(case["initial"]))],
                             colors=STYLE["initial_contour"],
                             linewidths=STYLE["initial_contour_width"],
                             linestyles=[(0, (3, 2.4))], zorder=2)
        initial.set_clip_path(disk)


def flow_arrow(ax, points, fraction, head_scale=1.0, halo_width=0.0):
    """A sharp concave dart, tangent to the actual curve and sized in points.

    The numerical curve supplies the shaft. The four-vertex filled head has
    swept-back wings and a central notch, rather than a flat triangular base.
    """
    lengths = np.r_[0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    if lengths[-1] < .10:
        return
    tip_s = fraction*lengths[-1]
    tip = [np.interp(tip_s, lengths, points[:, j]) for j in range(2)]
    near_s = max(0, tip_s - .01)
    near = [np.interp(near_s, lengths, points[:, j]) for j in range(2)]
    tip_screen = ax.transData.transform(tip)
    direction = tip_screen - ax.transData.transform(near)
    direction /= np.linalg.norm(direction)
    normal = np.array([-direction[1], direction[0]])
    pixels_per_point = ax.figure.dpi / 72 * head_scale
    length = STYLE["arrow_head_length_pt"] * pixels_per_point
    half_width = STYLE["arrow_head_width_pt"] / 2 * pixels_per_point
    notch = STYLE["arrow_head_notch_pt"] * pixels_per_point
    vertices = np.array([
        tip_screen,
        tip_screen - length*direction + half_width*normal,
        tip_screen - notch*direction,
        tip_screen - length*direction - half_width*normal,
    ])
    head = Polygon(ax.transData.inverted().transform(vertices), closed=True,
                   facecolor=STYLE["black"], edgecolor="none", zorder=4)
    head.set_gid("euclidean-arrowhead")
    if halo_width > 0:
        head.set_path_effects([patheffects.withStroke(
            linewidth=halo_width, foreground="white", joinstyle="miter")])
    ax.add_patch(head)


def euclidean_flow(ax, data, case_id):
    times, points = data["flow_times"], data["flow_points"]
    # Each located event is included exactly in the saved time grid. Split at
    # every event so c's post-reentry section is physical/solid again.
    cuts = np.r_[times[0], data["crossing_times"], times[-1]]
    for segment_index, (start, end) in enumerate(zip(cuts[:-1], cuts[1:])):
        segment = points[(times >= start) & (times <= end)]
        if len(segment) < 2:
            continue
        outside = np.linalg.norm(segment[len(segment)//2]) > 1
        ax.plot(segment[:, 0], segment[:, 1], color=STYLE["black"],
                lw=STYLE["flow_width"], ls=(0, (2.6, 2)) if outside else "solid", zorder=3)
        if STYLE["show_flow_arrows"]:
            fraction = STYLE["arrow_fraction"]
            head_scale = 1.0
            if segment_index == 0:
                fraction = STYLE["first_arrow_fraction_by_panel"].get(case_id, fraction)
                head_scale = STYLE["first_arrow_head_scale_by_panel"].get(case_id, head_scale)
            fraction = STYLE["arrow_fraction_by_segment"].get(case_id, {}).get(segment_index, fraction)
            flow_arrow(ax, segment, fraction, head_scale,
                       STYLE["arrow_head_halo_pt_by_panel"].get(case_id, 0.0))


def km_markers(ax, points, state_points, case_id):
    """Thin markers in screen space; every retained marker is an actual iterate."""
    coords = ax.transData.transform(points)
    endpoints = ax.transData.transform(state_points)
    pixels_per_point = ax.figure.dpi / 72
    separation = STYLE["marker_spacing_points"]*pixels_per_point
    clearance = STYLE["marker_endpoint_clearance_points"]*pixels_per_point
    arrow_clearance = (STYLE["marker_size"]/2 + STYLE["marker_edge_width"]/2 +
                       STYLE["marker_arrow_clearance_points"])*pixels_per_point
    arrow_paths = [patch.get_path().transformed(patch.get_transform())
                   for patch in ax.patches if patch.get_gid() == "euclidean-arrowhead"
                   and case_id in STYLE["marker_avoid_arrows_in"]]

    def clears_arrowheads(point):
        for path in arrow_paths:
            if path.contains_point(point):
                return False
            vertices = path.vertices
            starts, ends = vertices[:-1], vertices[1:]
            edges = ends - starts
            lengths_squared = np.sum(edges*edges, axis=1)
            fractions = np.clip(np.sum((point-starts)*edges, axis=1) /
                                np.maximum(lengths_squared, 1e-30), 0, 1)
            nearest = starts + fractions[:, None]*edges
            if np.min(np.linalg.norm(nearest-point, axis=1)) < arrow_clearance:
                return False
        return True

    selected, previous = [], coords[0]
    for k in range(1, len(points)-1):
        if (np.linalg.norm(coords[k]-previous) >= separation and
                np.min(np.linalg.norm(endpoints-coords[k], axis=1)) >= clearance and
                clears_arrowheads(coords[k])):
            selected.append(k)
            previous = coords[k]
    ax.plot(points[selected, 0], points[selected, 1], ls="", marker="o",
            ms=STYLE["marker_size"], mfc="white", mec=STYLE["km"],
            mew=STYLE["marker_edge_width"], zorder=6)
    if STYLE["show_iteration_labels"]:
        placements = ITERATION_LABELS[case_id]
        for k in selected:
            offset, ha, va = placements.get(k, placements["default"])
            number = ax.annotate(str(k), points[k], xytext=offset,
                                 textcoords="offset points", ha=ha, va=va,
                                 fontsize=STYLE["iteration_label_size"],
                                 color=STYLE["km"], zorder=13,
                                 annotation_clip=False)
            number.set_gid(f"km-iteration-{case_id}-{k}")
            number.set_path_effects([patheffects.withStroke(
                linewidth=STYLE["iteration_label_halo_width"], foreground="white")])
    return selected


def panel(ax, case, model, data):
    contours(ax, case, model)
    ax.add_patch(Circle((0, 0), 1, fill=False, color=STYLE["boundary"],
                        lw=STYLE["boundary_width"], zorder=2.5))
    euclidean_flow(ax, data, case["id"])
    km = data["km_points"]
    ax.plot(km[:, 0], km[:, 1], color=STYLE["km"], lw=STYLE["km_width"], zorder=3.5)
    state_points = [case["initial"], data["mle"], km[-1]]
    state_points.extend(data["crossing_points"])
    if np.linalg.norm(model.center) > 1:
        state_points.append(model.center)
    if STYLE["show_km_markers"]:
        km_markers(ax, km, np.asarray(state_points), case["id"])
    labels = POINT_LABELS[case["id"]]
    ax.scatter(*case["initial"], s=STYLE["state_dot_size"], color="black", linewidths=0, zorder=10)
    label(ax, case["initial"], labels["initial"])
    mle_handle = ax.scatter(*data["mle"], s=STYLE["mle_star_size"], marker="*", color="black",
                            edgecolors="white", linewidths=.2, zorder=11, label="MLE")
    label(ax, data["mle"], labels["mle"])
    if len(data["crossing_points"]):
        hit = data["crossing_points"][0]
        ax.scatter(*hit, s=STYLE["state_dot_size"], color="black", linewidths=0, zorder=10)
        label(ax, hit, labels["hit"])
    if np.linalg.norm(model.center) > 1:
        ax.scatter(*model.center, s=STYLE["state_dot_size"], color="black", linewidths=0, zorder=10)
        label(ax, model.center, labels["outside"])
    if np.linalg.norm(km[-1]-data["mle"]) > .025:
        ax.scatter(*km[-1], s=STYLE["state_dot_size"]*1.1, color=STYLE["km"], linewidths=0, zorder=10)
        label(ax, km[-1], labels["km"], STYLE["km"])
    return mle_handle


def render(cases, outputs, dpi):
    fonts()
    width, height = STYLE["figure_width_mm"], STYLE["figure_height_mm"]
    fig = plt.figure(figsize=(width/25.4, height/25.4), dpi=150)
    spans = [np.ptp(PANEL_X_RANGES[case["id"]]) for case, _, _ in cases]
    usable = width-STYLE["left_margin_mm"]-STYLE["right_margin_mm"]-(len(cases)-1)*STYLE["panel_gap_mm"]
    scale = usable/sum(spans)  # identical physical scale, not stretched circles
    left = STYLE["left_margin_mm"]
    axes = []
    for (case, model, data), span in zip(cases, spans):
        panel_left = left + STYLE.get("panel_x_offsets_mm", {}).get(case["id"], 0.0)
        ax = fig.add_axes([panel_left/width, STYLE["plot_bottom_mm"]/height,
                           scale*span/width, scale*np.ptp(STYLE["z_range"])/height])
        axes.append(ax)
        left += scale*span+STYLE["panel_gap_mm"]
        ax.set(xlim=PANEL_X_RANGES[case["id"]], ylim=STYLE["z_range"], aspect="equal")
        # Keep data transforms and equal aspect, but draw no axes or ticks.
        ax.set_axis_off()
        add_panel_letter(ax, case["id"], STYLE["panel_letter_size"],
                         position=STYLE["panel_letter_position"])
    fig.canvas.draw()  # Fix transforms before sizing arrowheads in points.
    for ax, (case, model, data) in zip(axes, cases):
        panel(ax, case, model, data)
    outputs.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "svg", "png"):
        fd, temp = tempfile.mkstemp(prefix="euclidean_km_", suffix=".tmp", dir=outputs)
        os.close(fd)
        destination = outputs / f"euclidean_km_bloch.{extension}"
        try:
            # Preserve the exact column width: tight cropping/padding would alter it.
            fig.savefig(temp, format=extension, dpi=dpi)
            os.replace(temp, destination)
        finally:
            Path(temp).unlink(missing_ok=True)
        print(destination, flush=True)
    plt.close(fig)
