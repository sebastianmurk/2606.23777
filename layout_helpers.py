"""Panel-letter placement; appearance only, independent of numerical data."""


def add_panel_letter(ax, case_id, fontsize, *, position=(0.02, 0.98)):
    """Place only '(a)', '(b)' or '(c)' at the upper left of a panel.

    ``position`` uses axes fractions, so placement follows the panel rather
    than a particular Bloch vector. The text is anchored at its upper left.
    """
    letter = ax.text(*position, f"({case_id})", transform=ax.transAxes,
                     ha="left", va="top", fontsize=fontsize, color="black",
                     clip_on=False, zorder=20)
    letter.set_gid(f"panel-letter-{case_id}")
    return letter
