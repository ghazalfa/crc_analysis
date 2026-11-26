"""Shared figure style so every stage colours clusters and text the same way (light mode)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
EMPTY = "#ecebe7"
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
RED_RAMP = ["#fbd5d4", "#f3a3a2", "#e34948", "#a92e2d"]
NEUTRAL_MID = "#f0efec"

# clusters avoid the blue / orange / aqua used for mutation classes and up / down genes
CLUSTER_COLORS = [SERIES[6], SERIES[4], SERIES[5], SERIES[7], SERIES[3], SERIES[2], SERIES[0], SERIES[1]]
CLUSTER_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
NO_DRIVER = "C0_no_driver"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
    "xtick.color": MUTED, "ytick.color": INK2, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "axes.spines.top": False, "axes.spines.right": False,
})


def cluster_colors(clusters):
    """C1, C2, ... always get the same colour; the no-driver group is grey."""
    named = sorted(c for c in clusters if c != NO_DRIVER)
    colors = {c: CLUSTER_COLORS[i % len(CLUSTER_COLORS)] for i, c in enumerate(named)}
    colors[NO_DRIVER] = MUTED
    return colors


def cluster_markers(clusters):
    named = sorted(c for c in clusters if c != NO_DRIVER)
    return {c: CLUSTER_MARKERS[i % len(CLUSTER_MARKERS)] for i, c in enumerate(named)}
