"""
Fe-C cooling explorer (Streamlit app).

Shows an approximate iron-carbon (Fe-Fe3C) phase diagram with a dot that
tracks an alloy as it cools. Three panels update together:

    1. The phase diagram with the alloy's tie line at the current temperature.
    2. A schematic micrograph. The minority phase grows from grain boundaries
       as its lever-rule fraction increases, with a color legend beside it.
    3. A lever-rule seesaw. Each phase is a circle whose AREA is proportional
       to its mass fraction, sitting on a beam balanced at the alloy C0.

The micrograph and seesaw share one figure, so their combined height equals
the height of the phase diagram.

Run with:  streamlit run fe_c_cooling.py

Note: the boundary curves are hand-placed approximations, traced by eye from
a reference diagram. They are not digitized data.
"""

import time
from io import BytesIO
from pathlib import Path

import numpy as np
import matplotlib

# "Agg" is a non-interactive backend. Streamlit renders figures to images,
# so no GUI window is needed. This must be set before importing pyplot.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Rectangle  # noqa: E402
import streamlit as st  # noqa: E402

# Must be the first Streamlit command in the script.
st.set_page_config(page_title="Fe-C cooling explorer", layout="wide",
                  initial_sidebar_state="expanded")


# ---------------------------------------------------------------------------
# Physical constants (all temperatures in °C, compositions in wt% C)
# ---------------------------------------------------------------------------
T_ROOM = 25.0           # room temperature, end point of the cooling path
T_EUT = 727.0           # eutectoid temperature
T_EUTECTIC = 1147.0     # eutectic temperature (liquid -> austenite + Fe3C)
C_EUT = 0.76            # eutectoid composition (austenite -> ferrite + Fe3C)
C_ECT = 4.30            # eutectic composition
C_FERRITE_MAX = 0.022   # max carbon in ferrite, reached at 727 °C
C_CEM = 6.67            # carbon in cementite (Fe3C), a fixed compound

# ---------------------------------------------------------------------------
# Figure layout, in inches. Images are shown at DPI pixels per inch.
# ---------------------------------------------------------------------------
DPI = 100
DIAG_W, DIAG_H = 6.2, 5.0            # phase diagram panel (larger)
RIGHT_W = 3.4                        # width of the micrograph + seesaw panel
MICRO_D = 1.07                       # micrograph diameter
TITLE_H = 0.30                       # room above the micrograph for its title
SEE_H = DIAG_H - MICRO_D - TITLE_H   # seesaw height, so the stack matches
                                     # the diagram's height exactly


# ---------------------------------------------------------------------------
# Approximate boundary data, traced from the reference Fe-C diagram.
#
# Each pair of arrays describes one boundary curve as (composition, temperature)
# points. np.interp uses them to look up a boundary value between points.
#
# Replace these with digitized data before using numbers for grading.
# ---------------------------------------------------------------------------
# Liquidus: upper boundary of the liquid + austenite region.
LIQ_C = np.array([0.0, 0.17, 1.0, 2.0, 3.0, 4.30])
LIQ_T = np.array([1538, 1493, 1430, 1300, 1220, 1147.0])

# Solidus: lower boundary of the liquid + austenite region (steel side).
SOL_C = np.array([0.0, 0.17, 2.14])
SOL_T = np.array([1538, 1493, 1147.0])

# Delta-ferrite boundaries (A4 line and peritectic). Drawn for reference only.
# The phase logic does not model delta-ferrite and treats it as austenite.
DELTA_C = np.array([0.0, 0.09])
DELTA_T = np.array([1394, 1493.0])
PERI_C = np.array([0.09, 0.17])
PERI_T = np.array([1493, 1493.0])

# A3 (GS line): austenite / ferrite + austenite boundary on the hypoeutectoid side.
A3_T = np.array([727, 755, 790, 820, 855, 912.0])
A3_C = np.array([0.76, 0.60, 0.40, 0.20, 0.10, 0.0])

# Acm (SE line): austenite / austenite + cementite boundary on the hypereutectoid side.
ACM_C = np.array([0.76, 1.0, 1.2, 1.5, 1.8, 2.14])
ACM_T = np.array([727, 780, 850, 950, 1050, 1147.0])

# Ferrite solvus: carbon solubility limit in ferrite, falls as T drops.
ALPHA_T = np.array([25, 400, 600, 727, 912.0])
ALPHA_C = np.array([0.0008, 0.003, 0.010, 0.022, 0.0])

# Cast-iron side: liquidus extension from the eutectic toward Fe3C.
CEM_LIQ_C = np.array([C_ECT, C_CEM])
CEM_LIQ_T = np.array([T_EUTECTIC, 1250.0])

# Human-readable phase names used in labels and the readout.
DISPLAY = {
    "liquid": "liquid",
    "austenite": "γ austenite",
    "ferrite": "α ferrite",
    "cementite": "Fe₃C",
    "pearlite": "pearlite",
}

# Standard teaching colors for the micrograph and legend.
AUSTENITE_RGB = (1.00, 0.85, 0.20)   # yellow
FERRITE_RGB = (0.35, 0.55, 0.90)     # blue
CEMENTITE_RGB = (0.05, 0.05, 0.05)   # black

MICRO_COLORS = {
    "liquid": (1.0, 1.0, 1.0),
    "austenite": AUSTENITE_RGB,
    "ferrite": FERRITE_RGB,
    "cementite": CEMENTITE_RGB,
}
# Pearlite is drawn as stripes of ferrite (light) and Fe3C (dark).
PEARL_DARK = CEMENTITE_RGB
PEARL_LIGHT = FERRITE_RGB


def disp(name):
    """Return the display label for a phase key, such as 'ferrite' -> 'α ferrite'.

    Falls back to the raw key if it is not in DISPLAY.
    """
    return DISPLAY.get(name, name)


# ---------------------------------------------------------------------------
# Lever rule and phase state
# ---------------------------------------------------------------------------
def lever(C, left, right):
    """Apply the lever rule to a two-phase region.

    The lever rule says the mass fraction of each phase is set by the arm on
    the opposite side of the alloy composition C:

        f_left  = (C_right - C) / (C_right - C_left)
        f_right = (C - C_left)  / (C_right - C_left)

    Args:
        C: Overall alloy composition (wt% C).
        left: (phase_name, composition) of the left tie-line end.
        right: (phase_name, composition) of the right tie-line end.

    Returns:
        dict with:
            "mix":  [(name, fraction), ...] for both phases.
            "tie":  ((name, comp, fraction), (name, comp, fraction)) for the
                    tie-line ends. Used for drawing and the seesaw.
    """
    (nl, cl), (nr, cr) = left, right

    # Clip to [0, 1] for safety. max(..., 1e-9) avoids division by zero when
    # both ends have the same composition.
    fl = float(np.clip((cr - C) / max(cr - cl, 1e-9), 0.0, 1.0))
    fr = 1.0 - fl  # the two fractions must sum to 1

    return {
        "mix": [(nl, fl), (nr, fr)],
        "tie": ((nl, cl, fl), (nr, cr, fr)),
    }


def two_phase(label, C, left, right):
    """Build the state dict for a two-phase region.

    Wraps lever() and adds a human-readable label. The microconstituent list
    equals the phase mix here because there is no further decomposition.
    """
    out = lever(C, left, right)
    out["label"] = label
    out["micro"] = out["mix"]
    return out


def single(label, phase):
    """Build the state dict for a single-phase region (fraction 1.0, no tie line)."""
    return {
        "label": label,
        "mix": [(phase, 1.0)],
        "tie": None,
        "micro": [(phase, 1.0)],
    }


def micro_below_727(C):
    """Return the microconstituents below 727 °C, referenced to the 727 °C split.

    Pearlite is a microconstituent (a lamellar mix of ferrite and Fe3C), not an
    equilibrium phase. It forms from austenite at the eutectoid, so its
    fraction is found with the lever rule across the 727 °C split.

    Returns:
        None if C <= 0.022 (the phase mix is already the full answer).
        Otherwise a list of (name, fraction) for the microstructure.
    """
    if C <= C_FERRITE_MAX:
        return None

    # Hypoeutectoid (0.022 < C < 0.76): pearlite plus proeutectoid ferrite.
    if C <= C_EUT:
        p = (C - C_FERRITE_MAX) / (C_EUT - C_FERRITE_MAX)
        return [("pearlite", p), ("ferrite", 1.0 - p)]

    # Hypereutectoid (0.76 < C): pearlite plus proeutectoid cementite.
    pro = (C - C_EUT) / (C_CEM - C_EUT)
    return [("cementite", pro), ("pearlite", 1.0 - pro)]


def state(C, T):
    """Determine the phases present for alloy C at temperature T.

    The function checks regions from the top of the diagram down:
        1. Liquid, and liquid + austenite between the liquidus and solidus.
        2. Solid above the eutectoid temperature (austenite, ferrite +
           austenite, or austenite + cementite).
        3. Below 727 °C (ferrite, or ferrite + cementite).

    Each branch returns a state dict (see single() and two_phase()).
    """
    # --- Region 1: liquid -------------------------------------------------
    # Above the liquidus, the alloy is fully liquid.
    if T >= np.interp(C, LIQ_C, LIQ_T):
        return single("Liquid", "liquid")

    # Between the liquidus and solidus, liquid and austenite coexist.
    if T >= np.interp(C, SOL_C, SOL_T):
        # Look up the liquid and solid compositions at this temperature.
        # The arrays are reversed because np.interp needs increasing x values.
        cl = float(np.interp(T, LIQ_T[::-1], LIQ_C[::-1]))
        cs = float(np.interp(T, SOL_T[::-1], SOL_C[::-1]))
        return two_phase("Liquid + austenite", C,
                         ("austenite", cs), ("liquid", cl))

    # --- Region 2: solid, above the eutectoid temperature -----------------
    if T >= T_EUT:
        if C < C_EUT:
            # Hypoeutectoid alloy. Fully austenite above the GS line (A3).
            if T >= np.interp(C, A3_C[::-1], A3_T[::-1]):
                return single("Austenite", "austenite")

            # Below the GS line, ferrite appears. If the alloy carbon is at or
            # below the ferrite solvus, the alloy is pure ferrite.
            ca = float(np.interp(T, ALPHA_T, ALPHA_C))
            if C <= ca:
                return single("Ferrite", "ferrite")

            # Otherwise ferrite and austenite coexist between the solvus and GS.
            cg = float(np.interp(T, A3_T, A3_C))
            return two_phase("Ferrite + austenite", C,
                             ("ferrite", ca), ("austenite", cg))

        # Hypereutectoid alloy. Fully austenite above the SE line (Acm).
        if T >= np.interp(C, ACM_C, ACM_T):
            return single("Austenite", "austenite")

        # Below the SE line, cementite precipitates alongside austenite.
        cg = float(np.interp(T, ACM_T, ACM_C))
        return two_phase("Austenite + cementite", C,
                         ("austenite", cg), ("cementite", C_CEM))

    # --- Region 3: solid, below 727 °C ------------------------------------
    # Austenite has transformed. The equilibrium pair is ferrite + Fe3C.
    ca = float(np.interp(T, ALPHA_T, ALPHA_C))
    if C <= ca:
        return single("Ferrite", "ferrite")

    out = two_phase("Ferrite + cementite", C,
                    ("ferrite", ca), ("cementite", C_CEM))

    # Use the 727 °C microconstituents for the picture. If micro_below_727
    # returns None (very low carbon), fall back to the phase mix.
    micro = micro_below_727(C)
    out["micro"] = micro if micro else out["mix"]
    return out


# ---------------------------------------------------------------------------
# Microstructure schematic
# ---------------------------------------------------------------------------

def _grain_edges(g):
    """Return a boolean mask of pixels that touch a different grain.

    A pixel is on a boundary if its grain label differs from the pixel
    below it or to its right.
    """
    edge = np.zeros(g.shape, dtype=bool)
    edge[:-1, :] |= g[:-1, :] != g[1:, :]
    edge[:, :-1] |= g[:, :-1] != g[:, 1:]
    return edge

# Phase that grows from grain boundaries when two solid phases coexist.
# The first entry present in a region is the boundary phase. This is set by
# identity, not by fraction, so the colors do not swap when fractions cross 50%.
BOUNDARY_PRIORITY = ["cementite", "ferrite", "austenite", "pearlite"]


def micro_image(micro, size=200, n=18, seed=7):
    """Draw a schematic micrograph inside a circle.

    Three cases:
        1. Pure liquid: plain white, no grain boundaries.
        2. Liquid + austenite: austenite nucleates at each grain center and
           grows outward into the liquid. Pixels nearest a center transform
           first. Boundaries are drawn only inside the solid.
        3. Solid phases: the boundary phase (set by BOUNDARY_PRIORITY) sits on
           grain boundaries first and grows inward. Its area fraction equals
           its lever-rule fraction. Pearlite gets stripes at each grain's angle.

    Args:
        micro: list of (name, fraction) pairs, fractions summing to about 1.
        size: image resolution in pixels (size x size).
        n: number of grains (Voronoi seed points).
        seed: random seed, so the same inputs always give the same picture.

    Returns:
        An RGB float array of shape (size, size, 3) with values in [0, 1].
    """
    fractions = dict(micro)

    # Grain centers and a stripe angle for each grain.
    rng = np.random.default_rng(seed)
    seeds = rng.uniform(-1, 1, (n, 2))
    angles = rng.uniform(0, np.pi, n)

    # Pixel grid covering [-1, 1] x [-1, 1].
    xs = np.linspace(-1, 1, size)
    X, Y = np.meshgrid(xs, xs)
    inside = X ** 2 + Y ** 2 <= 1.0     # pixels inside the circle

    # Squared distance from every pixel to every seed, shape (size, size, n).
    d = (X[..., None] - seeds[:, 0]) ** 2 + (Y[..., None] - seeds[:, 1]) ** 2

    # g: grain label of each pixel (index of its nearest seed).
    g = d.argmin(axis=-1)

    # Stripe pattern for pearlite, using each pixel's own grain angle.
    stripes = np.sin((X * np.cos(angles[g]) +
                      Y * np.sin(angles[g])) * 40) > 0

    img = np.ones((size, size, 3))      # start white

    def paint(mask, name):
        """Color the pixels in `mask` with the color of phase `name`."""
        if name == "pearlite":
            img[mask & stripes] = PEARL_DARK
            img[mask & ~stripes] = PEARL_LIGHT
        else:
            img[mask] = MICRO_COLORS[name]

    # --- Case 1 and 2: liquid is present -----------------------------------
    if "liquid" in fractions:
        f_aus = fractions.get("austenite", 0.0)

        # Pure liquid: nothing to draw, so the circle stays white.
        if f_aus <= 0.0:
            img[~inside] = 1.0
            return img

        # Distance from each pixel to its own grain center. It is 0 at the
        # center and largest at the grain boundary.
        d_own = np.sqrt(d.min(axis=-1))

        # Austenite fills the pixels nearest the centers. The cutoff is the
        # f_aus quantile of d_own over the circle, so the austenite area
        # equals its lever-rule fraction.
        cutoff = np.quantile(d_own[inside], min(f_aus, 1.0))
        solid = inside & (d_own <= cutoff)

        img[solid] = AUSTENITE_RGB
        img[_grain_edges(g) & solid] = (0.2, 0.2, 0.2)
        img[~inside] = 1.0
        return img

    # --- Case 3: solid phases, boundary-first growth ------------------------
    # Only one phase present: paint the whole circle and draw the grain lines.
    if len(micro) == 1:
        paint(inside, micro[0][0])
        img[_grain_edges(g) & inside] = (0.2, 0.2, 0.2)
        img[~inside] = 1.0
        return img

    # Pick the boundary phase by priority, not by fraction. The other phase
    # becomes the matrix that fills the grain interiors.
    boundary_name = next(p for p in BOUNDARY_PRIORITY if p in fractions)
    matrix_name = next(name for name, _ in micro if name != boundary_name)
    f_boundary = fractions[boundary_name]

    # Gap is near zero on grain boundaries and large inside grains.
    two_nearest = np.sqrt(np.partition(d, 1, axis=-1)[..., :2])
    gap = two_nearest.max(axis=-1) - two_nearest.min(axis=-1)

    # The boundary phase takes the pixels with the smallest gaps, so it sits
    # on the grain boundaries first and grows inward.
    is_boundary = np.zeros((size, size), dtype=bool)
    if f_boundary > 0.0:
        inside_idx = np.flatnonzero(inside)
        k = int(round(f_boundary * inside_idx.size))
        if k > 0:
            order = np.argsort(gap.ravel()[inside_idx])
            chosen = inside_idx[order[:k]]
            is_boundary.flat[chosen] = True

    paint(inside & ~is_boundary, matrix_name)
    paint(inside & is_boundary, boundary_name)

    # Grain boundaries across the whole solid circle.
    img[_grain_edges(g) & inside] = (0.2, 0.2, 0.2)
    img[~inside] = 1.0
    return img

# ---------------------------------------------------------------------------
# Drawing functions
# ---------------------------------------------------------------------------
# Drawing-only stretch near zero composition. The phase logic never uses this.
DISPLAY_STRETCH = 3.0        # how much the 0 to 0.1 wt% C band is widened
DISPLAY_STRETCH_MAX = 0.1    # compositions above this are drawn to true scale


def xdisp(c):
    """Map a true wt% C value to an x-position on the plot.

    Compositions up to 0.1 wt% C are multiplied by DISPLAY_STRETCH, so the
    ferrite solvus is visible near the y-axis. Above 0.1 the scale is linear
    again, shifted by the same offset so the curves stay continuous.

    Accepts a scalar or a numpy array.
    """
    c = np.asarray(c, dtype=float)
    offset = DISPLAY_STRETCH_MAX * (DISPLAY_STRETCH - 1.0)
    return np.where(c <= DISPLAY_STRETCH_MAX,
                    c * DISPLAY_STRETCH,
                    c + offset)

def draw_diagram(ax, C0, T0, T, s):
    """Draw the phase diagram, cooling path, and current tie line.

    Every composition is passed through xdisp() before plotting, so the
    region near 0 wt% C is widened. Temperatures are drawn as-is.

    Args:
        ax: matplotlib Axes to draw on.
        C0: alloy composition (wt% C).
        T0: starting temperature (°C), the top of the cooling path.
        T: current temperature (°C).
        s: state dict from state(C0, T).
    """
    solid = dict(color="0.25", lw=1.2)
    dashed = dict(color="0.45", lw=0.9, ls="--")

    # --- Steel side -------------------------------------------------------
    ax.plot(xdisp(LIQ_C), LIQ_T, **solid)          # liquidus
    ax.plot(xdisp(SOL_C), SOL_T, **solid)          # solidus
    ax.plot(xdisp(DELTA_C), DELTA_T, **solid)      # A4 (delta boundary)
    ax.plot(xdisp(PERI_C), PERI_T, **solid)        # peritectic isotherm
    ax.plot(xdisp(A3_C), A3_T, **solid)            # GS line (A3)
    ax.plot(xdisp(ACM_C), ACM_T, **solid)          # SE line (Acm)
    ax.plot(xdisp(ALPHA_C), ALPHA_T, **solid)      # ferrite solvus

    # Eutectoid isotherm (727 °C), from the ferrite limit to the Fe3C line.
    ax.hlines(T_EUT, xdisp(C_FERRITE_MAX), xdisp(C_CEM), **solid)

    # --- Cast-iron side ---------------------------------------------------
    ax.hlines(T_EUTECTIC, xdisp(2.14), xdisp(C_CEM), **solid)
    ax.plot([xdisp(2.14)] * 2, [T_EUT, T_EUTECTIC], **dashed)
    ax.plot([xdisp(C_ECT)] * 2, [T_EUT, T_EUTECTIC], **dashed)
    ax.plot(xdisp(CEM_LIQ_C), CEM_LIQ_T, **dashed)

    # Vertical Fe3C line at 6.67 wt% C, from room temperature up.
    ax.plot([xdisp(C_CEM)] * 2, [T_ROOM, CEM_LIQ_T[-1]],
            color="0.25", lw=1.6)
    ax.text(xdisp(C_CEM) - 0.1, 1350, "Fe₃C\n6.67 %", fontsize=7, ha="right")

    # --- Temperature and point labels --------------------------------------
    for txt, x, y, ha in [("1538 °C", 0.05, 1555, "left"),
                          ("912 °C", 0.05, 935, "left"),
                          ("727 °C", 6.6, 745, "right"),
                          ("1147 °C", 6.6, 1165, "right")]:
        ax.text(x, y, txt, fontsize=6, ha=ha)

    # --- Region labels -----------------------------------------------------
    for txt, x, y in [("L", 3.4, 1350), ("γ + L", 0.95, 1385),
                      ("δ", 0.05, 1440), ("γ", 0.5, 1000),
                      ("α + Fe₃C", 1.5, 250)]:
        ax.text(x, y, txt, fontsize=7, color="0.35")

    # Dotted vertical line: the path the alloy follows as it cools.
    ax.plot([xdisp(C0)] * 2, [T_ROOM, T0], color="tab:red", lw=1, ls=":")

    # In a two-phase region, draw the tie line and mark its two ends.
    if s["tie"]:
        (nl, cl, _), (nr, cr, _) = s["tie"]
        xl, xr = xdisp(cl), xdisp(cr)
        ax.plot([xl, xr], [T, T], color="black", lw=1.5)
        ax.plot([xl, xr], [T, T], "o", color="black", ms=5)
        ax.text(xl, T + 40, disp(nl), fontsize=7, ha="center")
        ax.text(xr, T + 40, disp(nr), fontsize=7, ha="center")

    # The alloy's current point, drawn on top of everything else.
    ax.plot(xdisp(C0), T, "o", color="tab:red", ms=7, zorder=5)

    # --- Axes: ticks show true composition values at stretched positions ---
    ax.set_xlim(0, float(xdisp(6.8)))
    ax.set_ylim(0, 1600)
    tick_c = np.arange(7)
    ax.set_xticks(xdisp(tick_c))
    ax.set_xticklabels([str(c) for c in tick_c])
    ax.set_yticks(range(0, 1601, 200))
    ax.set_xlabel("wt% C (stretched near 0)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title(f"Fe-C diagram: T = {T:.0f} °C", fontsize=10)
    ax.grid(alpha=0.2)

def draw_seesaw(ax, C0, s):
    """Draw the lever-rule seesaw in inch units.

    The axes must be set to (0, RIGHT_W) by (0, SEE_H), so one data unit
    equals one inch and the circles are truly circular.

    Each phase end carries a circle whose AREA is proportional to its mass
    fraction (radius = r_max * sqrt(f)). The fill uses the phase color from
    MICRO_COLORS. The beam pivots at C0.

    Args:
        ax: matplotlib Axes in inch coordinates.
        C0: alloy composition (wt% C), used to place the pivot.
        s: state dict from state(C0, T). If it has no tie line, a message is drawn.
    """
    W, H = RIGHT_W, SEE_H

    # Single-phase region: no tie line, so no seesaw.
    if not s["tie"]:
        ax.text(W / 2, H / 2, "Single phase: no lever rule needed",
                ha="center", va="center", fontsize=8)
        return

    (nl, cl, fl), (nr, cr, fr) = s["tie"]

    # Arm lengths in composition units.
    arm_l, arm_r = C0 - cl, cr - C0

    # Layout in inches. The beam sits in the upper part of the panel.
    beam_y = SEE_H * 0.6
    max_arm_in = 1.1
    r_max = 0.42

    # Scale composition arms to inches. The pivot sits at the horizontal center.
    scale = max_arm_in / max(arm_l, arm_r, 1e-9)
    xp = W / 2
    xl = xp - arm_l * scale
    xr = xp + arm_r * scale

    # Beam and fulcrum.
    ax.plot([xl, xr], [beam_y, beam_y], color="0.3", lw=4,
            solid_capstyle="round")
    ax.plot([xp], [beam_y - 0.15], marker="^", ms=10, color="0.3")

    # One circle per phase, resting on the beam. Area is proportional to
    # the fraction. The fill uses the micrograph colors.
    for x, name, f in ((xl, nl, fl), (xr, nr, fr)):
        r = r_max * np.sqrt(f)
        fill = MICRO_COLORS.get(name, (0.85, 0.85, 0.85))
        ax.add_patch(Circle((x, beam_y + r), r, fc=fill, ec="0.3", lw=1))
        ax.text(x, beam_y + 2 * r + 0.06, f"{disp(name)}\n{f:.1%}",
                ha="center", va="bottom", fontsize=7)

    # Arm labels at the midpoint of each arm, just below the beam.
    ax.text((xl + xp) / 2, beam_y - 0.4, f"arm\n{arm_l:.3f}",
            ha="center", va="center", fontsize=7)
    ax.text((xp + xr) / 2, beam_y - 0.4, f"arm\n{arm_r:.3f}",
            ha="center", va="center", fontsize=7)

    # Rule reminder at the bottom of the panel.
    ax.text(W / 2, 0.08, "mass fraction × its arm must match on both sides",
            ha="center", va="bottom", fontsize=6.5)

def draw_legend(ax):
    """Draw a color key in inch units beside the micrograph.

    The axes must be sized to (RIGHT_W, MICRO_D) inches so data units are inches.

    Args:
        ax: matplotlib Axes covering the full width of the right panel.
    """
    ax.set_xlim(0, RIGHT_W)
    ax.set_ylim(0, MICRO_D)
    ax.set_axis_off()

    items = [("γ austenite", AUSTENITE_RGB),
             ("α ferrite", FERRITE_RGB),
             ("Fe₃C", CEMENTITE_RGB)]

    # Stack the swatches from the top down. Each swatch is a small square
    # with a label to its right.
    for i, (label, color) in enumerate(items):
        y = MICRO_D - 0.25 - i * 0.30
        ax.add_patch(Rectangle((0.05, y), 0.22, 0.20,
                               fc=color, ec="0.2", lw=0.6))
        ax.text(0.35, y + 0.10, label, va="center", fontsize=7)


def draw_right_panel(C0, s):
    """Build the right-hand figure: micrograph, legend, and seesaw.

    The figure is DIAG_H inches tall, the same as the phase diagram, so the
    two images stack to the same height on screen.

    Args:
        C0: alloy composition (wt% C), passed to the seesaw.
        s: state dict from state(C0, T).

    Returns:
        A matplotlib Figure. The caller shows it and closes it.
    """
    fig = plt.figure(figsize=(RIGHT_W, DIAG_H), dpi=DPI)

    # Micrograph axes: a square MICRO_D inches wide, centered horizontally,
    # sitting just above the seesaw area and below the title strip.
    micro_bottom = SEE_H / DIAG_H
    micro_left = (RIGHT_W - MICRO_D) / 2 / RIGHT_W
    micro_ax = fig.add_axes([micro_left, micro_bottom,
                             MICRO_D / RIGHT_W, MICRO_D / DIAG_H])
    micro_ax.imshow(micro_image(s["micro"]), extent=(-1, 1, -1, 1),
                    origin="lower")
    micro_ax.add_patch(Circle((0, 0), 1, fill=False, lw=0.9, ec="0.1"))
    micro_ax.set_xlim(-1.05, 1.05)
    micro_ax.set_ylim(-1.05, 1.05)
    micro_ax.set_axis_off()
    micro_ax.set_title("Microstructure (schematic)", fontsize=7, pad=2)

    # Legend axes: same height as the micrograph, spanning the full width.
    # The swatches are drawn in the empty space beside the circle.
    legend_ax = fig.add_axes([0, SEE_H / DIAG_H, 1, MICRO_D / DIAG_H])
    draw_legend(legend_ax)

    # Seesaw axes: fills the bottom SEE_H inches. Data units are inches.
    see_ax = fig.add_axes([0, 0, 1, SEE_H / DIAG_H])
    see_ax.set_xlim(0, RIGHT_W)
    see_ax.set_ylim(0, SEE_H)
    draw_seesaw(see_ax, C0, s)

    return fig


def show(placeholder, fig):
    """Save a figure to an in-memory PNG and display it at its native size.

    Showing at a fixed pixel width keeps the displayed heights exact. Streamlit
    would otherwise rescale each image to its column width.

    Args:
        placeholder: st.empty() slot to fill.
        fig: matplotlib Figure to show. It is closed after saving.
    """
    buf = BytesIO()
    width_px = int(round(fig.get_figwidth() * DPI))
    fig.savefig(buf, format="png", dpi=DPI)
    plt.close(fig)
    buf.seek(0)
    placeholder.image(buf, width=width_px)


def readout(T, C0, s):
    """Build the Markdown text shown under the panels.

    Lists the phase fractions, the tie-line compositions, a torque check, and
    (if different) the microconstituent fractions.

    Args:
        T: current temperature (°C).
        C0: alloy composition (wt% C).
        s: state dict from state(C0, T).

    Returns:
        A Markdown string for st.markdown.
    """
    lines = [f"**T = {T:.0f} °C · {s['label']}**", "",
             "**Phase fractions (lever rule)**"]

    # Skip fractions that round to zero (below 0.05%).
    for name, f in s["mix"]:
        if f > 0.0005:
            lines.append(f"- {disp(name)}: **{f:.1%}**")

    if s["tie"]:
        lines += ["", "**Phase compositions (tie-line ends)**"]
        for name, c, _ in s["tie"]:
            lines.append(f"- {disp(name)}: {c:.4f} wt% C")

        # Torque check: fraction × arm on each side should match.
        (nl, cl, fl), (nr, cr, fr) = s["tie"]
        lines += [
            "",
            "**Balance check (torques about C₀)**",
            f"- Left: {fl:.3f} × ({C0:.3f} − {cl:.4f}) = {fl * (C0 - cl):.5f}",
            f"- Right: {fr:.3f} × ({cr:.4f} − {C0:.3f}) = {fr * (cr - C0):.5f}",
        ]

    # Show the microconstituents only when they differ from the phase mix,
    # such as pearlite below 727 °C.
    if s["micro"] != s["mix"]:
        lines += ["", "**Microconstituents (schematic)**"]
        for name, f in s["micro"]:
            if f > 0.0005:
                lines.append(f"- {disp(name)}: {f:.1%}")

    return "\n".join(lines)

# ---------------------------------------------------------------------------
# README tab: reads README.md from the same folder as this script
# ---------------------------------------------------------------------------
README_PATH = Path(__file__).with_name("README.md")


def load_readme():
    """Return the README text, or a short note if the file is missing."""
    if README_PATH.exists():
        return README_PATH.read_text(encoding="utf-8")
    return ("README.md was not found next to fe_c_cooling.py. "
            "Put it in the same folder and reload the page.")


# ---------------------------------------------------------------------------
# UI: sidebar controls, tabs, and page layout
# ---------------------------------------------------------------------------
st.title("Fe-C cooling explorer")
st.caption("Approximate Fe-Fe₃C diagram for teaching. Boundaries are "
           "hand-placed, not exact. Composition runs from pure iron up to "
           "the steel/cast-iron boundary (2.14 wt% C).")

# Animation controls: where the cooling path starts, and how it steps.
st.sidebar.header("Cooling animation")
T0 = st.sidebar.slider("Initial temperature T₀ (°C)", 800, 1390, 1000, 5)
step = st.sidebar.select_slider("Cooling step (°C per frame)",
                                options=[5, 10, 20, 40], value=20)
delay = st.sidebar.slider("Frame delay (s)", 0.0, 0.3, 0.05, 0.01)

# Composition of the alloy, shared by the animation and the static view.
st.sidebar.header("Alloy")
C0 = st.sidebar.slider("Composition C₀ (wt% C)", 0.00, 2.14, 0.40, 0.01)

# Static view: the temperature shown when Play is not pressed.
st.sidebar.header("Current view")
T_NOW = st.sidebar.slider("Current temp (°C)", int(T_ROOM), 1390, 1000, 5)

tab_app, tab_readme = st.tabs(["Explorer", "README"])

with tab_readme:
    st.markdown(load_readme())

with tab_app:
    # Button returns True only on the run where it was clicked.
    play = st.button(f"▶ Play cooling from {T0} °C to {T_ROOM:.0f} °C")

    # Two columns: the diagram on the left, the micrograph, legend, and
    # seesaw on the right. st.empty() creates placeholders that render()
    # overwrites.
    left, right = st.columns([1.3, 1])
    ph_diag = left.empty()
    ph_right = right.empty()
    ph_text = st.empty()

    def render(T):
        """Draw the diagram, the right-hand panel, and the readout for temperature T.

        Each figure is created, shown in its placeholder, and then closed.
        Closing frees memory, which matters during an animation with many frames.
        """
        s = state(C0, T)

        # Panel 1: phase diagram.
        fig, ax = plt.subplots(figsize=(DIAG_W, DIAG_H), dpi=DPI)
        draw_diagram(ax, C0, T0, T, s)
        fig.tight_layout()
        show(ph_diag, fig)

        # Panels 2 and 3: micrograph and legend above, seesaw below.
        show(ph_right, draw_right_panel(C0, s))

        # Text readout below the panels.
        ph_text.markdown(readout(T, C0, s))

    # -----------------------------------------------------------------------
    # Main flow: animate on button press, otherwise show the current temperature.
    # -----------------------------------------------------------------------
    if play:
        # Build the temperature list from T0 down to room temperature, in
        # steps of `step`. The final frame is T_ROOM, so the animation always
        # ends there.
        frames = list(np.arange(T0, T_ROOM, -step)) + [T_ROOM]
        for T in frames:
            render(float(T))
            time.sleep(delay)  # pause so the eye can follow each frame
    else:
        # No animation: show the state at the temperature chosen on the slider.
        render(float(T_NOW))
    # No animation: show the state at the temperature chosen on the slider.
    render(float(T_NOW))
#-----End-of-file----------------------
