import os
import numpy
import pandas as pd
import matplotlib.pyplot as plt


def angle_diff_deg(a_deg, b_deg):
    """
    Smallest signed difference a-b in degrees, in [-180, 180).
    Works with scalars, numpy arrays, pandas Series.
    """
    return ((a_deg - b_deg + 180.0) % 360.0) - 180.0


def circ_signed_error_deg(est_deg, theory_deg):
    """
    Absolute circular error magnitude in [0, 180].
    """
    return angle_diff_deg(est_deg, theory_deg)


def signed_error_to_radius(err_deg, err_min=-190.0, err_max=190.0):
    """
    Map signed error [err_min..err_max] to radius [1..0].
    err_min (-190) -> r=1 (rim)
    0 -> r=0.5 (middle)
    err_max (+190) -> r=0 (center)
    """
    e = float(err_deg)
    e = min(max(e, err_min), err_max)
    frac = (e - err_min) / (err_max - err_min)  # 0..1
    return 1.0 - frac


def plot_event_on_ax(
    ax,
    event_name,
    event_data,
    method,
    cc_threshold,
    err_max_deg,
    bar_width_deg,
    alpha_bar,
    lw_whisk,
    use_percentile_spread=True,
    spread_percentiles=(16, 84),
    draw_theory_marker=True,
    label=None,
):
    """
    Plot one event as:
      - bar at theta=true_baz with length encoding median(abs error)
      - radial whisker encoding spread of ABS error distribution
    """

    if method == "mid":
        # Match your linear code logic:
        ed = event_data[event_data["corr_value_max"] > cc_threshold]
        est_baz = ed["baz_value_midpoint"]
        true_baz = ed["true_baz"]

    elif method == "max":
        # Match your linear code logic:
        ed = event_data[event_data["corr_value_midpoint"] > cc_threshold]
        est_baz = ed["baz_value_max"]
        true_baz = ed["true_baz"]
    else:
        return

    if len(ed) == 0:
        return

    mask = numpy.isfinite(est_baz.values) & numpy.isfinite(true_baz.values)
    if mask.sum() == 0:
        return

    est_vals = est_baz.values[mask].astype(float)
    true_vals = true_baz.values[mask].astype(float)

    # theoretical baz per event (robust)
    theoretical_baz_deg = float(numpy.nanmedian(true_vals))
    theta_theory = numpy.deg2rad(theoretical_baz_deg)

    # errors relative to ONE theoretical value (important!)
    errs = circ_signed_error_deg(est_vals, theoretical_baz_deg)

    # BOX POSITION = median error (majority)
    r0 = signed_error_to_radius(0.0, err_min=err_min_deg, err_max=err_max_deg)
    box_err_deg = float(numpy.nanpercentile(errs, 50))
    r_box = signed_error_to_radius(box_err_deg, err_min=err_min_deg, err_max=err_max_deg)

    p_lo, p_hi = numpy.nanpercentile(errs, spread_percentiles)
    r_lo = signed_error_to_radius(float(p_lo), err_min=err_min_deg, err_max=err_max_deg)
    r_hi = signed_error_to_radius(float(p_hi), err_min=err_min_deg, err_max=err_max_deg)

    r_bottom = min(r0, r_box)
    bar_height = abs(r_box - r0)

    # BAR: from median-error radius to rim
    ax.bar(
        theta_theory,
        bar_height,
        width=numpy.deg2rad(bar_width_deg),
        bottom=r_bottom,
        alpha=alpha_bar,
        edgecolor="k",
        label=label,
    )

    # Whisker: radial spread (no label so legend doesn't duplicate)
    ax.plot(
        [theta_theory, theta_theory],
        [r_lo, r_hi],
        linewidth=lw_whisk,
        label=None,
    )

    # Optional theoretical marker at rim
    if draw_theory_marker:
        ax.plot([theta_theory], [r0], marker="o", markersize=2, linestyle="None", label=None)


# ------------------- user config -------------------
root = "/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/correlation_baz_csv"
root_savefig = "/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/Backazimuth_estimations/statistics"

freq_tag = "0.01_0.1"  # choose one frequency type
mincorr = 0.8          # ONLY THIS ONE now

err_min_deg = -190.0
err_max_deg = 190.0
bar_width_deg = 3.0

use_percentile_spread = True
spread_percentiles = (16, 84)

locations = ["LIV_external", "LIV_internal", "HAN_external"]
methods = ["mid", "max"]


for location in locations:
    for method in methods:

        # ------------------- read all events for this location -------------------
        data_all_events = {}
        freqmin = None
        freqmax = None

        for csv in os.listdir(root):
            if not csv.endswith(".csv"):
                continue
            if freq_tag not in csv:
                continue
            if location not in csv:
                continue

            df = pd.read_csv(os.path.join(root, csv))

            folder = csv.split("_")[0]
            freqmin = csv.split("_")[-6]
            freqmax = csv.split("_")[-5]

            data_all_events[folder] = df

        if len(data_all_events) == 0:
            continue

        # ------------------- one circular figure with ALL events -------------------
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection="polar")

        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)

        # Radial axis labels in "error degrees" (edge=-190, center=190)
        error_ticks = [-190, -180, -90, 0, 90, 180, 190]
        rticks = [signed_error_to_radius(e, err_min=err_min_deg, err_max=err_max_deg) for e in error_ticks]
        ax.set_yticks(rticks)
        ax.set_yticklabels([f"{e}°" for e in error_ticks])
        ax.set_ylim(0, 1)

        ax.set_title(
            f"{location} | {method}points | CC>{mincorr} | {freqmin} - {freqmax} Hz\n"
            f"Plotted at theoretical BAz; radial distance = |error|",
            pad=20
        )

        # Plot all events ONCE (only CC>0.8)
        for event_name, df in data_all_events.items():
            plot_event_on_ax(
                ax=ax,
                event_name=event_name,
                event_data=df,
                method=method,
                cc_threshold=mincorr,
                err_max_deg=err_max_deg,
                bar_width_deg=bar_width_deg,
                alpha_bar=0.35,
                lw_whisk=1.2,
                use_percentile_spread=use_percentile_spread,
                spread_percentiles=spread_percentiles,
                draw_theory_marker=True,
                label=event_name,  # <-- legend shows events
            )

        # Event legend (can get large!)
        ax.legend(
            loc="lower left",
            bbox_to_anchor=(1.05, 0.05),
            frameon=False,
            ncol=2,
            fontsize=8,
            title="Events"
        )

        plt.tight_layout()

        outpath = f"{root_savefig}/circular_BAz_Corr_{location}_{method}_allevents_CC{mincorr}.png"
        plt.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)

plt.show()


'''import numpy
import matplotlib.pyplot as plt

def angle_diff_deg(a_deg, b_deg):
    """
    Smallest signed difference a-b in degrees, in [-180, 180).
    """
    return ((a_deg - b_deg + 180.0) % 360.0) - 180.0

def circ_signed_error_deg(est_deg, theory_deg):
    """
    Absolute circular error magnitude in [0, 180].
    """
    return abs(angle_diff_deg(est_deg, theory_deg))

def error_to_radius(err_deg, err_max=270.0):
    """
    Map error (0..err_max) to radius (1..0).
    0 error -> r=1 (edge)
    err_max -> r=0 (center)
    """
    err_clipped = min(max(err_deg, 0.0), err_max)
    return 1.0 - (err_clipped / err_max)



# ---------------- example input ----------------
for
theoretical_baz_deg = 325.0
estimated_mean_baz_deg = 200.0
estimated_spread_deg = 70.0          # +/- around mean estimate (in azimuth)
err_max_deg = 270.0                  # edge->center scale

# Mean error magnitude (circular)
mean_err_deg = circ_signed_error_deg(estimated_mean_baz_deg, theoretical_baz_deg)
r_mean = error_to_radius(mean_err_deg, err_max=err_max_deg)

# Convert the +/- spread in estimated BAz into a *radial* spread in ERROR magnitude
est_low_deg  = (estimated_mean_baz_deg - estimated_spread_deg) % 360.0
est_high_deg = (estimated_mean_baz_deg + estimated_spread_deg) % 360.0

err_low_deg  = circ_signed_error_deg(est_low_deg,  theoretical_baz_deg)
err_high_deg = circ_signed_error_deg(est_high_deg, theoretical_baz_deg)

# Whisker should span min..max error magnitude
err_min = min(err_low_deg, err_high_deg)
err_max = max(err_low_deg, err_high_deg)

r_err_min = error_to_radius(err_min, err_max=err_max_deg)  # closer to rim
r_err_max = error_to_radius(err_max, err_max=err_max_deg)  # closer to center

# Polar angle where we plot everything (theoretical BAz)
theta_theory = numpy.deg2rad(theoretical_baz_deg)

# ---------------- plotting ----------------
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection="polar")

# Put 0° at North and increase clockwise (common for azimuth)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

# Radial axis labels in "error degrees" (edge=0, center=270)
error_ticks = [0, 90, 180, 270]
rticks = [error_to_radius(e, err_max=err_max_deg) for e in error_ticks]
ax.set_yticks(rticks)
ax.set_yticklabels([f"{e}°" for e in error_ticks])
ax.set_ylim(0, 1)

# 1) Theoretical backazimuth marker on rim (0 error)
ax.plot([theta_theory], [1.0], marker="o", markersize=5, linestyle="None",
        label=f"Theoretical BAz = {theoretical_baz_deg:.0f}°")

# 2) Mean error as a radial "box/bar" at theoretical BAz
bar_height = 1.0 - r_mean
ax.bar(theta_theory, bar_height, width=numpy.deg2rad(4.0), bottom=r_mean,
       alpha=0.35, edgecolor="k",
       label=f"Mean est. BAz={estimated_mean_baz_deg:.0f}° → |Δ|={mean_err_deg:.0f}°")

# 3) Error spread as a RADIAL whisker at theoretical BAz
ax.plot([theta_theory, theta_theory], [r_err_max, r_err_min], linewidth=1,
        label=f"Spread from est ±{estimated_spread_deg:.0f}° → error [{err_min:.0f}°, {err_max:.0f}°]")

# Optional: marker at mean-error radius (still at theoretical angle)
#ax.plot([theta_theory], [r_mean], marker="s", markersize=7, linestyle="None")

ax.set_title("Backazimuth: plot at theoretical BAz; radial distance = error", pad=20)
ax.legend(loc="lower left", bbox_to_anchor=(1.05, 0.05))
plt.tight_layout()
plt.show()'''