"""
Fe-C cooling explorer (Streamlit app).

Shows an approximate iron-carbon (Fe-Fe3C) phase diagram with a dot that
tracks an alloy as it cools. Three panels update together:

    1. The phase diagram with the alloy's tie line at the current temperature.
    2. A schematic micrograph. The minority phase grows from grain boundaries
       as its lever-rule fraction increases, with a color legend beside it.
    3. A lever-rule seesaw. Each phase is a circle whose AREA is proportional
       to its mass fraction, sitting on a beam balanced at the alloy C0.

The micrograph and seesaw share one right-hand figure so their combined height
equals the height of the phase diagram. Both panels are built as lightweight
inline SVGs (with the micrograph embedded as a small PNG whose Voronoi grain
geometry is precomputed once), avoiding heavy plotting libraries so the app
loads quickly and uses minimal memory.

Run with:  streamlit run fe_c_cooling_v3.py

Note: the boundary curves are hand-placed approximations, traced by eye from
a reference diagram. They are not digitized data.
"""

import base64
import html
import struct
import time
import zlib
from pathlib import Path

import numpy as np
import streamlit as st

# Must be the first Streamlit command in the script.
st.set_page_config(
    page_title="Fe-C cooling explorer",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Adjust the vertical padding of the main Streamlit container
st.markdown(
    """
    <style>
    .block-container {
        padding-top: 1.2rem;
        padding-bottom: 0.5rem;
    }
    /* Optional: tighten vertical gaps between widgets */
    [data-testid="stVerticalBlock"] > div {
        gap: 0.4rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

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
# SVG canvas geometry and font sizes (all units in SVG pixels)
# Both panels share HEIGHT = 500 so they align to the same height on screen.
# ---------------------------------------------------------------------------
HEIGHT = 500
DIAG_W = 620            # phase diagram panel (wider)
RIGHT_W = 380           # micrograph + legend + seesaw panel

# Plot area margins inside the phase diagram SVG
PAD_L, PAD_R, PAD_T, PAD_B = 72, 20, 42, 58
PLOT_W = DIAG_W - PAD_L - PAD_R
PLOT_H = HEIGHT - PAD_T - PAD_B

# Micrograph circle and seesaw layout inside the right-hand SVG
MICRO_CX, MICRO_CY, MICRO_R = 252, 108, 66   # circle center and radius
SEESAW_TOP = 192                             # top of the seesaw region
BEAM_Y = 365                                 # vertical position of the seesaw beam
R_MAX = 46                                   # max circle radius (at 100% fraction)
MAX_ARM_PX = 112                             # max half-span of the seesaw beam

# Typography (change any size here to scale labels across the app)
FONT_FAMILY = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif"
FS_TITLE = 18       # panel titles
FS_AXIS = 15        # axis titles (wt% C, Temperature)
FS_TICK = 13        # axis tick numbers
FS_REGION = 15      # phase region labels (L, γ + L, γ, α + Fe₃C)
FS_TEMP_REF = 13    # critical temperature & Fe3C callouts
FS_TIE = 14         # tie-line endpoint labels
FS_LEGEND = 15      # micrograph color key labels
FS_SEESAW = 15      # seesaw phase name & percentage labels
FS_ARM = 14         # seesaw arm length labels
FS_NOTE = 13.5      # seesaw bottom rule reminder


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

# Standard teaching colors in 8-bit RGB (0-255) for the micrograph, plus hex
# strings derived from them for the SVG legend and seesaw circles.
AUSTENITE_RGB = (255, 217, 51)     # yellow
FERRITE_RGB = (89, 140, 230)       # blue
CEMENTITE_RGB = (13, 13, 13)       # black
LIQUID_RGB = (255, 255, 255)       # white
EDGE_RGB = (51, 51, 51)            # dark grey grain boundaries

MICRO_COLORS = {
    "liquid": LIQUID_RGB,
    "austenite": AUSTENITE_RGB,
    "ferrite": FERRITE_RGB,
    "cementite": CEMENTITE_RGB,
}
PHASE_HEX = {name: f"#{r:02x}{g:02x}{b:02x}" for name, (r, g, b) in MICRO_COLORS.items()}

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

# Phase that grows from grain boundaries when two solid phases coexist.
# The first entry present in a region is the boundary phase. This is set by
# identity, not by fraction, so the colors do not swap when fractions cross 50%.
BOUNDARY_PRIORITY = ["cementite", "ferrite", "austenite", "pearlite"]


@st.cache_resource
def _micro_grid(size=160, n=18, seed=7):
    """Precompute the static Voronoi grain layout once per app process.

    Because the grain seeds, stripe orientations, and boundary distances never
    change between temperatures or compositions, caching them in float32 avoids
    recomputing (size x size x n) distance matrices on every animation frame.
    """
    rng = np.random.default_rng(seed)
    seeds = rng.uniform(-1, 1, (n, 2)).astype(np.float32)
    angles = rng.uniform(0, np.pi, n).astype(np.float32)

    # Pixel grid covering [-1, 1] x [-1, 1].
    xs = np.linspace(-1, 1, size, dtype=np.float32)
    X, Y = np.meshgrid(xs, xs)
    inside = X ** 2 + Y ** 2 <= 1.0     # pixels inside the circle
    inside_idx = np.flatnonzero(inside)

    # Squared distance from every pixel to every seed, shape (size, size, n).
    d = (X[..., None] - seeds[:, 0]) ** 2 + (Y[..., None] - seeds[:, 1]) ** 2

    # g: grain label of each pixel (index of its nearest seed).
    g = d.argmin(axis=-1)

    # A pixel is on a boundary if its grain label differs from its neighbor
    # below or to the right.
    edges = np.zeros(g.shape, dtype=bool)
    edges[:-1, :] |= g[:-1, :] != g[1:, :]
    edges[:, :-1] |= g[:, :-1] != g[:, 1:]
    edges &= inside

    # Stripe pattern for pearlite, using each pixel's own grain angle.
    stripes = np.sin((X * np.cos(angles[g]) + Y * np.sin(angles[g])) * 40) > 0

    # Distance from each pixel to its own grain center (0 at center, max at edge).
    d_own = np.sqrt(d.min(axis=-1))
    d_own_inside = d_own[inside]

    # Gap between the two nearest seeds is near zero on grain boundaries and
    # large inside grains. Pre-sorting inside pixels lets each frame slice [:k].
    two_nearest = np.sqrt(np.partition(d, 1, axis=-1)[..., :2])
    gap = two_nearest.max(axis=-1) - two_nearest.min(axis=-1)
    gap_order = inside_idx[np.argsort(gap.ravel()[inside_idx])]

    return size, inside, inside_idx, edges, stripes, d_own, d_own_inside, gap_order


def micro_image(micro):
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

    Returns:
        A uint8 RGB array of shape (size, size, 3) with values in [0, 255].
    """
    size, inside, inside_idx, edges, stripes, d_own, d_own_inside, gap_order = _micro_grid()
    fractions = dict(micro)
    img = np.full((size, size, 3), 255, dtype=np.uint8)  # start white

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
            return img

        # Austenite fills the pixels nearest the centers. The cutoff is the
        # f_aus quantile of d_own over the circle, so the austenite area
        # equals its lever-rule fraction.
        cutoff = np.quantile(d_own_inside, min(f_aus, 1.0))
        solid = inside & (d_own <= cutoff)

        img[solid] = AUSTENITE_RGB
        img[edges & solid] = EDGE_RGB
        return img

    # --- Case 3: solid phases, boundary-first growth ------------------------
    # Only one phase present: paint the whole circle and draw the grain lines.
    if len(micro) == 1:
        paint(inside, micro[0][0])
        img[edges] = EDGE_RGB
        return img

    # Pick the boundary phase by priority, not by fraction. The other phase
    # becomes the matrix that fills the grain interiors.
    boundary_name = next(p for p in BOUNDARY_PRIORITY if p in fractions)
    matrix_name = next(name for name, _ in micro if name != boundary_name)
    f_boundary = fractions[boundary_name]

    # The boundary phase takes the pixels with the smallest boundary gaps, so
    # it sits on the grain boundaries first and grows inward.
    is_boundary = np.zeros((size, size), dtype=bool)
    if f_boundary > 0.0:
        k = int(round(f_boundary * inside_idx.size))
        if k > 0:
            is_boundary.flat[gap_order[:k]] = True

    paint(inside & ~is_boundary, matrix_name)
    paint(inside & is_boundary, boundary_name)

    # Grain boundaries across the whole solid circle.
    img[edges] = EDGE_RGB
    return img


def _png_chunk(tag, data):
    """Build one PNG chunk (4-byte length, 4-byte type, payload, CRC32)."""
    return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)


def _rgb_to_png_b64(img):
    """Encode a uint8 RGB array as a base64 PNG using only the stdlib."""
    h, w, _ = img.shape
    raw = np.empty((h, 1 + w * 3), dtype=np.uint8)
    raw[:, 0] = 0  # PNG filter type 0 (None) at the start of each scanline
    raw[:, 1:] = img.reshape(h, w * 3)
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)
    idat = zlib.compress(raw.tobytes(), level=1)
    png = b"\x89PNG\r\n\x1a\n" + _png_chunk(b"IHDR", ihdr) + _png_chunk(b"IDAT", idat) + _png_chunk(b"IEND", b"")
    return base64.b64encode(png).decode("ascii")


@st.cache_data(max_entries=64)
def _micro_png_b64(micro_key):
    """Cache encoded micrograph PNGs for repeated states (e.g. below 727 °C)."""
    return _rgb_to_png_b64(micro_image(list(micro_key)))


# ---------------------------------------------------------------------------
# Drawing helpers & SVG builders
# ---------------------------------------------------------------------------
# Drawing-only stretch near zero composition. The phase logic never uses this.
DISPLAY_STRETCH = 3.0        # how much the 0 to 0.1 wt% C band is widened
DISPLAY_STRETCH_MAX = 0.1    # compositions above this are drawn to true scale
X_MAX = float(6.8 + DISPLAY_STRETCH_MAX * (DISPLAY_STRETCH - 1.0))


def xdisp(c):
    """Map a true wt% C value to a stretched composition coordinate.

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


def _sx_disp(xd):
    """Convert an already-stretched composition value into an SVG x-coordinate."""
    return PAD_L + (float(xd) / X_MAX) * PLOT_W


def _sx(c):
    """Convert a true wt% C value into an SVG x-coordinate."""
    return _sx_disp(xdisp(c))


def _sy(t):
    """Convert a temperature (°C) into an SVG y-coordinate (0 °C at bottom)."""
    return PAD_T + ((1600.0 - float(t)) / 1600.0) * PLOT_H


def _polyline(cs, ts, stroke="#404040", width=1.6, dash=""):
    """Build an SVG <polyline> from (wt% C, °C) arrays."""
    pts = " ".join(f"{_sx(c):.1f},{_sy(t):.1f}" for c, t in zip(cs, ts))
    d_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline fill="none" stroke="{stroke}" stroke-width="{width}"{d_attr} points="{pts}"/>'


def _line(x1, y1, x2, y2, stroke="#404040", width=1.6, dash=""):
    """Build an SVG <line> between two pixel coordinates."""
    d_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{width}"{d_attr}/>'


def _text(x, y, content, size, anchor="middle", fill="#222", weight="normal", halo=False, baseline="auto"):
    """Build an escaped SVG <text> element, optionally with a white outline halo."""
    halo_attr = ' stroke="white" stroke-width="3.5" paint-order="stroke" stroke-linejoin="round"' if halo else ""
    w_attr = f' font-weight="{weight}"' if weight != "normal" else ""
    b_attr = f' dominant-baseline="{baseline}"' if baseline != "auto" else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" '
        f'text-anchor="{anchor}"{w_attr}{b_attr}{halo_attr}>{html.escape(content)}</text>'
    )


@st.cache_data(max_entries=1)
def _static_diagram_body():
    """Build the unchanging background SVG elements of the Fe-C diagram once.

    Grid lines, axes, boundary curves, and reference labels never move during
    an animation, so caching them lets each frame render only the cooling line,
    tie line, current dot, and temperature title.
    """
    el = []

    # --- Axes, grid, and ticks (show true wt% C at stretched positions) ---
    for t in range(0, 1601, 200):
        y = _sy(t)
        el.append(_line(PAD_L, y, PAD_L + PLOT_W, y, stroke="#e6e6e6", width=1.0))
        el.append(_text(PAD_L - 8, y, str(t), FS_TICK, anchor="end", fill="#444", baseline="central"))

    for c in range(7):
        x = _sx(c)
        el.append(_line(x, PAD_T, x, PAD_T + PLOT_H, stroke="#e6e6e6", width=1.0))
        el.append(_text(x, PAD_T + PLOT_H + 18, str(c), FS_TICK, anchor="middle", fill="#444"))

    el.append(
        f'<rect x="{PAD_L}" y="{PAD_T}" width="{PLOT_W}" height="{PLOT_H}" '
        f'fill="none" stroke="#333" stroke-width="1.3"/>'
    )

    el.append(_text(PAD_L + PLOT_W / 2, HEIGHT - 10, "wt% C (stretched near 0)", FS_AXIS, weight="600"))
    el.append(
        f'<text x="20" y="{PAD_T + PLOT_H / 2:.1f}" font-size="{FS_AXIS}" font-weight="600" '
        f'fill="#222" text-anchor="middle" transform="rotate(-90 20 {PAD_T + PLOT_H / 2:.1f})">Temperature (°C)</text>'
    )

    # --- Steel side -------------------------------------------------------
    el.append(_polyline(LIQ_C, LIQ_T))          # liquidus
    el.append(_polyline(SOL_C, SOL_T))          # solidus
    el.append(_polyline(DELTA_C, DELTA_T))      # A4 (delta boundary)
    el.append(_polyline(PERI_C, PERI_T))        # peritectic isotherm
    el.append(_polyline(A3_C, A3_T))            # GS line (A3)
    el.append(_polyline(ACM_C, ACM_T))          # SE line (Acm)
    el.append(_polyline(ALPHA_C, ALPHA_T))      # ferrite solvus

    # Eutectoid isotherm (727 °C), from the ferrite limit to the Fe3C line.
    el.append(_line(_sx(C_FERRITE_MAX), _sy(T_EUT), _sx(C_CEM), _sy(T_EUT)))

    # --- Cast-iron side ---------------------------------------------------
    el.append(_line(_sx(2.14), _sy(T_EUTECTIC), _sx(C_CEM), _sy(T_EUTECTIC)))
    el.append(_line(_sx(2.14), _sy(T_EUT), _sx(2.14), _sy(T_EUTECTIC), stroke="#737373", width=1.2, dash="5,4"))
    el.append(_line(_sx(C_ECT), _sy(T_EUT), _sx(C_ECT), _sy(T_EUTECTIC), stroke="#737373", width=1.2, dash="5,4"))
    el.append(_polyline(CEM_LIQ_C, CEM_LIQ_T, stroke="#737373", width=1.2, dash="5,4"))

    # Vertical Fe3C line at 6.67 wt% C, from room temperature up.
    el.append(_line(_sx(C_CEM), _sy(T_ROOM), _sx(C_CEM), _sy(CEM_LIQ_T[-1]), stroke="#333", width=2.0))
    el.append(_text(_sx(C_CEM) - 8, _sy(1380), "Fe₃C", FS_TEMP_REF, anchor="end", weight="600", halo=True))
    el.append(_text(_sx(C_CEM) - 8, _sy(1380) + 17, "6.67 %", FS_TEMP_REF, anchor="end", weight="600", halo=True))

    # --- Temperature and point labels -------------------------------------
    for txt, xd, t, anchor in [
        ("1538 °C", 0.14, 1545, "start"),
        ("912 °C", 0.12, 928, "start"),
        ("727 °C", 6.72, 745, "end"),
        ("1147 °C", 6.72, 1165, "end"),
    ]:
        el.append(_text(_sx_disp(xd), _sy(t), txt, FS_TEMP_REF, anchor=anchor, fill="#333", halo=True))

    # --- Region labels ----------------------------------------------------
    for txt, xd, t, size in [
        ("L", 3.55, 1360, FS_REGION),
        ("γ + L", 1.15, 1380, FS_REGION),
        ("δ", 0.11, 1435, FS_TEMP_REF),
        ("γ", 0.68, 1000, FS_REGION),
        ("α + Fe₃C", 2.0, 270, FS_REGION),
    ]:
        el.append(_text(_sx_disp(xd), _sy(t), txt, size, anchor="middle", fill="#444", weight="600", halo=True))

    return "\n".join(el)


def diagram_svg(C0, T0, T, s):
    """Build the phase diagram SVG, cooling path, and current tie line.

    Args:
        C0: alloy composition (wt% C).
        T0: starting temperature (°C), the top of the cooling path.
        T: current temperature (°C).
        s: state dict from state(C0, T).
    """
    dyn = [_static_diagram_body()]

    # Title showing current temperature
    dyn.append(_text(PAD_L + PLOT_W / 2, 26, f"Fe-C diagram: T = {T:.0f} °C", FS_TITLE, weight="700"))

    # Dotted vertical line: the path the alloy follows as it cools.
    x0 = _sx(C0)
    dyn.append(_line(x0, _sy(T_ROOM), x0, _sy(T0), stroke="#d62728", width=1.8, dash="3,4"))

    # In a two-phase region, draw the tie line and mark its two ends.
    if s["tie"]:
        (nl, cl, _), (nr, cr, _) = s["tie"]
        xl, xr = _sx(cl), _sx(cr)
        yt = _sy(T)
        dyn.append(_line(xl, yt, xr, yt, stroke="#111", width=2.2))
        dyn.append(f'<circle cx="{xl:.1f}" cy="{yt:.1f}" r="4.5" fill="#111"/>')
        dyn.append(f'<circle cx="{xr:.1f}" cy="{yt:.1f}" r="4.5" fill="#111"/>')

        y_lbl = max(yt - 11, PAD_T + 18)
        # Keep endpoint labels inside the plot box and clear of each other when narrow.
        a_l = "start" if xl < PAD_L + 42 else "middle"
        a_r = "end" if xr > PAD_L + PLOT_W - 42 else "middle"
        if xr - xl < 75:
            a_l, a_r = "end", "start"
        dyn.append(_text(xl, y_lbl, disp(nl), FS_TIE, anchor=a_l, fill="#111", weight="600", halo=True))
        dyn.append(_text(xr, y_lbl, disp(nr), FS_TIE, anchor=a_r, fill="#111", weight="600", halo=True))

    # The alloy's current point, drawn on top of everything else.
    dyn.append(f'<circle cx="{x0:.1f}" cy="{_sy(T):.1f}" r="6.5" fill="#d62728" stroke="white" stroke-width="1.5"/>')

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {DIAG_W} {HEIGHT}" '
        f'style="width:100%;height:auto;display:block;font-family:{FONT_FAMILY};">\n'
        + "\n".join(dyn)
        + "\n</svg>"
    )


def right_panel_svg(C0, s):
    """Build the right-hand SVG: micrograph, color legend, and lever seesaw.

    Both figures share HEIGHT = 500 so they align to the same height on screen.
    Each seesaw circle has AREA proportional to its mass fraction
    (radius = R_MAX * sqrt(f)), and the beam pivots at C0.

    Args:
        C0: alloy composition (wt% C), used to place the pivot.
        s: state dict from state(C0, T).
    """
    el = []

    # --- Microstructure + Legend ------------------------------------------
    el.append(_text(RIGHT_W / 2, 24, "Microstructure (schematic)", FS_TITLE - 1, weight="700"))

    micro_key = tuple((name, round(float(f), 4)) for name, f in s["micro"])
    b64 = _micro_png_b64(micro_key)
    x_img, y_img, d_img = MICRO_CX - MICRO_R, MICRO_CY - MICRO_R, MICRO_R * 2
    el.append(
        f'<defs><clipPath id="micro-clip">'
        f'<circle cx="{MICRO_CX}" cy="{MICRO_CY}" r="{MICRO_R}"/>'
        f'</clipPath></defs>'
    )
    el.append(
        f'<image x="{x_img}" y="{y_img}" width="{d_img}" height="{d_img}" '
        f'href="data:image/png;base64,{b64}" clip-path="url(#micro-clip)"/>'
    )
    el.append(
        f'<circle cx="{MICRO_CX}" cy="{MICRO_CY}" r="{MICRO_R}" fill="none" stroke="#222" stroke-width="1.5"/>'
    )

    # Color key stacked beside the circular micrograph.
    legend_items = [
        ("γ austenite", PHASE_HEX["austenite"]),
        ("α ferrite", PHASE_HEX["ferrite"]),
        ("Fe₃C", PHASE_HEX["cementite"]),
    ]
    for i, (label, hex_color) in enumerate(legend_items):
        y = 56 + i * 38
        el.append(f'<rect x="14" y="{y}" width="22" height="20" rx="2" fill="{hex_color}" stroke="#333" stroke-width="1"/>')
        el.append(_text(44, y + 10, label, FS_LEGEND, anchor="start", fill="#222", baseline="central"))

    # --- Lever-rule seesaw ------------------------------------------------
    # Single-phase region: no tie line, so no seesaw.
    if not s["tie"]:
        el.append(_text(RIGHT_W / 2, (SEESAW_TOP + HEIGHT) / 2, "Single phase: no lever rule needed", FS_SEESAW + 1, fill="#444"))
    else:
        (nl, cl, fl), (nr, cr, fr) = s["tie"]

        # Arm lengths in composition units, scaled so the longer arm fits MAX_ARM_PX.
        arm_l, arm_r = C0 - cl, cr - C0
        scale = MAX_ARM_PX / max(arm_l, arm_r, 1e-9)
        xp = RIGHT_W / 2
        xl = xp - arm_l * scale
        xr = xp + arm_r * scale

        # Beam and fulcrum triangle.
        el.append(
            f'<line x1="{xl:.1f}" y1="{BEAM_Y}" x2="{xr:.1f}" y2="{BEAM_Y}" '
            f'stroke="#444" stroke-width="6" stroke-linecap="round"/>'
        )
        el.append(
            f'<polygon points="{xp:.1f},{BEAM_Y + 5} {xp - 10:.1f},{BEAM_Y + 23} {xp + 10:.1f},{BEAM_Y + 23}" fill="#444"/>'
        )

        # One circle per phase resting on the beam (area proportional to fraction).
        for x, name, f in ((xl, nl, fl), (xr, nr, fr)):
            r = R_MAX * float(np.sqrt(f))
            fill = PHASE_HEX.get(name, "#d9d9d9")
            el.append(f'<circle cx="{x:.1f}" cy="{BEAM_Y - r - 3:.1f}" r="{r:.1f}" fill="{fill}" stroke="#333" stroke-width="1.2"/>')
            y_lbl = BEAM_Y - 2 * r - 12
            el.append(_text(x, y_lbl - 18, disp(name), FS_SEESAW, fill="#222"))
            el.append(_text(x, y_lbl, f"{f:.1%}", FS_SEESAW, fill="#111", weight="700"))

        # Arm labels at each arm's midpoint, kept clear of the center fulcrum.
        xl_arm = min((xl + xp) / 2, xp - 44)
        xr_arm = max((xp + xr) / 2, xp + 44)
        for xa, arm in ((xl_arm, arm_l), (xr_arm, arm_r)):
            el.append(_text(xa, BEAM_Y + 46, "arm", FS_ARM, fill="#444"))
            el.append(_text(xa, BEAM_Y + 64, f"{arm:.3f}", FS_ARM, fill="#222", weight="600"))

        # Rule reminder at the bottom of the panel.
        el.append(_text(RIGHT_W / 2, HEIGHT - 12, "mass fraction × its arm must match on both sides", FS_NOTE, fill="#444"))

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {RIGHT_W} {HEIGHT}" '
        f'style="width:100%;height:auto;display:block;font-family:{FONT_FAMILY};">\n'
        + "\n".join(el)
        + "\n</svg>"
    )


def show_svg(placeholder, svg, width_px):
    """Display an SVG inside a placeholder as a self-contained data-URI image."""
    b64 = base64.b64encode(svg.encode("utf-8")).decode("ascii")
    placeholder.markdown(
        f'<img src="data:image/svg+xml;base64,{b64}" '
        f'style="width:100%;max-width:{width_px}px;height:auto;display:block;">',
        unsafe_allow_html=True,
    )


def readout(T, C0, s):
    """Build the header and 4 side-by-side Markdown sections shown under the panels.

    Lists the phase fractions, the tie-line compositions, a torque check, and
    (if different) the microconstituent fractions.

    Args:
        T: current temperature (°C).
        C0: alloy composition (wt% C).
        s: state dict from state(C0, T).

    Returns:
        (header_md, [col1_md, col2_md, col3_md, col4_md]) for st.markdown.
    """
    header = f"**T = {T:.0f} °C · {s['label']}**"

    # 1. Phase fractions (lever rule)
    frac_lines = ["**Phase fractions (lever rule)**"]
    for name, f in s["mix"]:
        if f > 0.0005:
            frac_lines.append(f"- {disp(name)}: **{f:.1%}**")

    # 2. Phase compositions & 3. Balance check (only when a tie line exists)
    if s["tie"]:
        comp_lines = ["**Phase compositions (tie-line ends)**"]
        for name, c, _ in s["tie"]:
            comp_lines.append(f"- {disp(name)}: {c:.4f} wt% C")

        (nl, cl, fl), (nr, cr, fr) = s["tie"]
        bal_lines = [
            "**Balance check (torques about C₀)**",
            f"- Left: {fl:.3f} × ({C0:.3f} − {cl:.4f}) = {fl * (C0 - cl):.5f}",
            f"- Right: {fr:.3f} × ({cr:.4f} − {C0:.3f}) = {fr * (cr - C0):.5f}",
        ]
    else:
        comp_lines = []
        bal_lines = []

    # 4. Microconstituents (only when they differ from the phase mix, e.g. pearlite)
    if s["micro"] != s["mix"]:
        micro_lines = ["**Microconstituents (schematic)**"]
        for name, f in s["micro"]:
            if f > 0.0005:
                micro_lines.append(f"- {disp(name)}: {f:.1%}")
    else:
        micro_lines = []

    cols = [
        "\n".join(frac_lines),
        "\n".join(comp_lines),
        "\n".join(bal_lines),
        "\n".join(micro_lines),
    ]
    return header, cols


# ---------------------------------------------------------------------------
# README tab: reads README.md from the same folder as this script
# ---------------------------------------------------------------------------
README_PATH = Path(__file__).with_name("README.md")


@st.cache_data(max_entries=1)
def load_readme():
    """Return the README text, or a short note if the file is missing."""
    if README_PATH.exists():
        return README_PATH.read_text(encoding="utf-8")
    return ("README.md was not found next to fe_c_cooling_v3.py. "
            "Put it in the same folder and reload the page.")


# ---------------------------------------------------------------------------
# UI: sidebar controls, tabs, and page layout
# ---------------------------------------------------------------------------
st.title("Fe-C cooling explorer")
st.caption("Approximate Fe-Fe₃C diagram for teaching. Boundaries are "
           "hand-placed, not exact. Composition runs from pure iron up to "
           "the steel/cast-iron boundary (2.14 wt% C).")

# Animation controls: initial temperature, step/delay, and playback buttons.
st.sidebar.header("Cooling animation")
T0 = st.sidebar.slider("Initial temperature T₀ (°C)", 800, 1390, 1300, 5)
step = st.sidebar.select_slider("Cooling step (°C per frame)",
                                options=[5, 10, 20, 40], value=5)
delay = st.sidebar.slider("Frame delay (s)", 0.0, 0.3, 0.05, 0.01)

if "T_NOW" not in st.session_state:
    st.session_state.T_NOW = T0
if "playing" not in st.session_state:
    st.session_state.playing = False

# Apply any pending animation step before the T_NOW slider is instantiated.
if st.session_state.pop("advance", False) and st.session_state.playing:
    if st.session_state.T_NOW > T_ROOM:
        st.session_state.T_NOW = max(int(T_ROOM), st.session_state.T_NOW - step)
    else:
        st.session_state.playing = False


def _toggle_play(t0):
    st.session_state.playing = not st.session_state.playing
    if st.session_state.playing and st.session_state.T_NOW <= T_ROOM:
        st.session_state.T_NOW = t0


def _step_back(t0, s):
    st.session_state.playing = False
    st.session_state.T_NOW = min(t0, st.session_state.T_NOW + s)


def _step_forward(s):
    st.session_state.playing = False
    st.session_state.T_NOW = max(int(T_ROOM), st.session_state.T_NOW - s)


def _restart(t0):
    st.session_state.playing = False
    st.session_state.T_NOW = t0


b_play, b_back, b_fwd, b_rst = st.sidebar.columns(4)
b_play.button("⏸" if st.session_state.playing else "▶", help="Play / Pause",
              on_click=_toggle_play, args=(T0,))
b_back.button("⏮", help=f"Step back (+{step} °C)",
              on_click=_step_back, args=(T0, step))
b_fwd.button("⏭", help=f"Step forward (-{step} °C)",
             on_click=_step_forward, args=(step,))
b_rst.button("↺", help="Restart at T₀", on_click=_restart, args=(T0,))

# Static view: the composition and temperature shown when Play is not pressed.
st.sidebar.header("Current view")
C0 = st.sidebar.slider("Current composition C₀ (wt% C)", 0.00, 2.14, 1.50, 0.01)
T_NOW = st.sidebar.slider("Current temp (°C)", int(T_ROOM), 1390, step=5, key="T_NOW")
st.sidebar.info(f"**Current state:** {state(C0, T_NOW)['label']}")

tab_app, tab_readme = st.tabs(["Explorer", "README"])

with tab_readme:
    st.markdown(load_readme())

with tab_app:
    # Two columns: the diagram on the left, the micrograph, legend, and
    # seesaw on the right. st.empty() creates placeholders that render()
    # overwrites.
    left, right = st.columns([DIAG_W, RIGHT_W])
    ph_diag = left.empty()
    ph_right = right.empty()
    ph_header = st.empty()
    ph_cols = [c.empty() for c in st.columns(4)]

    def render(T):
        """Draw the diagram, the right-hand panel, and the readout for temperature T."""
        s = state(C0, T)
        show_svg(ph_diag, diagram_svg(C0, T0, T, s), DIAG_W)
        show_svg(ph_right, right_panel_svg(C0, s), RIGHT_W)
        header, cols = readout(T, C0, s)
        ph_header.markdown(header)
        for ph, text in zip(ph_cols, cols):
            ph.markdown(text)

    # -----------------------------------------------------------------------
    # Main flow: render current temperature, and advance if playing.
    # -----------------------------------------------------------------------
    render(float(T_NOW))

    if st.session_state.playing:
        if T_NOW > T_ROOM:
            time.sleep(delay)
            st.session_state.advance = True
            st.rerun()
        else:
            st.session_state.playing = False
            st.rerun()
