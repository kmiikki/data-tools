#!/usr/bin/env python3
"""
sliding_segmented_regression.ver4_onsets.py

Unified slope-based onset detector (multi-onset version).

- Computes local slope dS/dt from the chosen signal (MGI/BW or Total Counts)
- Finds nucleation onset t1 using ignition + stabilization delay
- Finds ALL significant bulk-growth bursts after t1 by scanning slope(t)
  for strong growth waves (generalization of old "t2")
- Reports and plots the first K onsets according to --onsets

Author: Kim Miikki et al., 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

# ============================================================
# Utility: path resolution
# ============================================================

def resolve_paths(infile_name, outfig_name):
    """
    Pi logic: if running in /opt/tools, data is assumed to be located in ~/rgb/analysis/.
    If running from Desktop/analysis, use relative paths.
    """
    if os.getcwd() == "/opt/tools":
        infile = os.path.expanduser("~/rgb/analysis/" + infile_name)
        outfig = os.path.expanduser("~/rgb/analysis/" + outfig_name)
    else:
        infile = infile_name
        outfig = outfig_name
    return infile, outfig

# ============================================================
# Core helpers
# ============================================================

def compute_local_slopes(df: pd.DataFrame,
                         time_col: str,
                         signal_col: str,
                         window_size_s: float = 100.0) -> pd.DataFrame:
    """
    Compute local linear slope dS/dt by fitting np.polyfit(time, signal, 1)
    in a symmetric sliding window of +/- window_size_s around each point.
    """
    t = df[time_col].to_numpy()
    y = df[signal_col].to_numpy()

    # median dt
    dt = np.median(np.diff(t))
    halfwin = max(1, int(window_size_s / dt))

    slopes = np.full_like(y, np.nan, dtype=float)
    for i in range(len(t)):
        lo = max(0, i - halfwin)
        hi = min(len(t), i + halfwin)
        if hi - lo < 5:
            continue
        m, c = np.polyfit(t[lo:hi], y[lo:hi], 1)
        slopes[i] = m

    out = df.copy()
    out["slope"] = slopes
    return out


def detect_t_end_from_slope(df: pd.DataFrame,
                            time_col: str,
                            slope_col: str = "slope",
                            slope_threshold: float = 0.002) -> float:
    """
    Find earliest time after which |slope| < slope_threshold for the remainder.
    If never fully plateaus, return last timestamp.
    """
    t = df[time_col].to_numpy()
    s = df[slope_col].to_numpy()
    for i in range(len(t)):
        if np.all(np.abs(s[i:]) < slope_threshold):
            return float(t[i])
    return float(t[-1])


def detect_t1_first_onset(df: pd.DataFrame,
                          time_col: str,
                          signal_col: str,
                          slope_col: str,
                          baseline_frac: float = 0.10,
                          k_sigma_primary: float = 4.0,
                          confirm_time_s: float = 60.0,
                          drop_tolerance_abs: float = 0.5,
                          stability_delay_s: float = 100.0):
    """
    Two-stage nucleation onset detector.

    Stage 1 (ignition, t1_raw):
      - earliest time where slope is >= (median_baseline_slope + k_sigma_primary * baseline_slope_std)
      - and the signal DOES NOT fall back down within confirm_time_s
        by more than drop_tolerance_abs in absolute units

    Stage 2 (report, t1_report):
      - shift t1_raw forward by stability_delay_s
        This corresponds to "visually stable nucleation onset" that you
        actually want to report in figures/papers.

    Returns (t1_raw, t1_report). If no onset found, returns (None, None).
    """

    t = df[time_col].to_numpy()
    s = df[slope_col].to_numpy()
    y = df[signal_col].to_numpy()

    # baseline region = first baseline_frac of data
    n_base = max(5, int(len(t) * baseline_frac))
    slope_base = s[:n_base]
    slope_mean = np.nanmedian(slope_base)
    slope_std  = np.nanstd(slope_base)
    slope_thresh = slope_mean + k_sigma_primary * slope_std

    dt = np.median(np.diff(t))
    confirm_n = max(2, int(confirm_time_s / dt))

    for i in range(len(t)):
        if not np.isfinite(s[i]):
            continue
        if s[i] < slope_thresh:
            continue

        # candidate ignition
        t_candidate = t[i]
        y_candidate = y[i]

        j_end = min(len(t), i + confirm_n)
        if j_end - i < confirm_n:
            break  # not enough future data to confirm

        y_future = y[i:j_end]

        # condition: does NOT collapse back down
        if np.all(y_future >= (y_candidate - drop_tolerance_abs)):
            t1_raw = float(t_candidate)
            t1_report = t1_raw + stability_delay_s
            return t1_raw, t1_report

    return None, None


def detect_growth_bursts(df: pd.DataFrame,
                         time_col: str,
                         slope_col: str,
                         t_start: float,
                         t_end: float,
                         min_separation_s: float = 180.0,
                         prominence_fraction: float = 0.8,
                         recovery_fraction: float = 0.2,
                         recovery_time_s: float = 200.0,
                         backtrack_s: float = 180.0,
                         backtrack_fraction: float = 0.3):
    """
    Detect major growth bursts and report their start times.

    backtrack_fraction:
        When moving backward from the slope peak, stop when slope < 
        backtrack_fraction * max_slope. This finds the *start of acceleration*,
        not the point where slope was already large.
    """

    t = df[time_col].to_numpy()
    s = df[slope_col].to_numpy()

    mask = (t >= t_start) & (t <= t_end)
    if not np.any(mask):
        return []

    t_seg = t[mask]
    s_seg = s[mask]
    s_max = np.nanmax(s_seg)
    if not np.isfinite(s_max) or s_max <= 0:
        return []

    high_thresh = prominence_fraction * s_max
    low_thresh = recovery_fraction * s_max
    backtrack_thresh = backtrack_fraction * s_max

    candidates = []
    last_time = None

    for i in range(1, len(s_seg) - 1):
        if not np.isfinite(s_seg[i]):
            continue
        if s_seg[i] < high_thresh:
            continue
        if not (s_seg[i] >= s_seg[i - 1] and s_seg[i] >= s_seg[i + 1]):
            continue

        t_peak = t_seg[i]

        # --- Backtrack to find where this burst started ---
        t_back_min = t_peak - backtrack_s
        j = i
        while j > 0 and t_seg[j] >= t_back_min and s_seg[j] >= backtrack_thresh:
            j -= 1
        t_start_burst = float(t_seg[j + 1])

        # --- Recovery check (between previous and this burst) ---
        if last_time is not None:
            mask_low = (t_seg >= last_time) & (t_seg <= t_start_burst)
            s_low = s_seg[mask_low]
            if len(s_low) == 0:
                continue
            if not np.any(s_low < low_thresh):
                continue
            if (t_start_burst - last_time) < min_separation_s:
                continue

        candidates.append(t_start_burst)
        last_time = t_start_burst

    return candidates


def time_to_temp(df: pd.DataFrame,
                 time_col: str,
                 temp_col: str,
                 t_value: float | None) -> float | None:
    """
    Given a time t_value, return the closest temperature reading from df.
    """
    if t_value is None:
        return None
    idx = (df[time_col] - t_value).abs().idxmin()
    return float(df.loc[idx, temp_col])


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Unified slope-based onset detection with multi-burst support (ver4)."
    )
    parser.add_argument(
        "--signal",
        type=str,
        default="BW",
        choices=["BW", "Total"],
        help="Which signal to analyze: BW (=MGI) or Total (=Total Counts)"
    )
    parser.add_argument(
        "--onsets",
        type=int,
        default=0,
        help=(
            "How many onsets to report:\n"
            " 0 = automatic (t1 + all bursts)\n"
            " 1 = t1 only\n"
            " 2 = t1 + first burst\n"
            " 3 = t1 + first two bursts\n"
            " etc."
        )
    )

    args = parser.parse_args()

    # Filenames
    infile_name = "merged_rgb_counts_time.csv"
    outfig_name = "onset_analysis_v4.png"

    infile, outfig = resolve_paths(infile_name, outfig_name)

    # Load CSV
    df_raw = pd.read_csv(infile, encoding="utf-8")

    # Time column auto-detect
    time_candidates = [
        "Time (s)",
        "time_s",
        "Time_s",
        "time (s)",
        "seconds",
        "t_s",
        "t (s)"
    ]
    time_col = None
    for cand in time_candidates:
        if cand in df_raw.columns:
            time_col = cand
            break
    if time_col is None:
        raise KeyError(
            f"Could not find a time column. Tried {time_candidates}, "
            f"got {list(df_raw.columns)}"
        )

    temp_col = "Tr (°C)"

    if args.signal == "BW":
        signal_col = "BW"            # legacy name, plotted as MGI
        y_label = "MGI"
    else:
        signal_col = "Total Counts"
        y_label = "Total Counts"

    # interpolate missing points
    df_raw[signal_col] = df_raw[signal_col].interpolate(method="linear")

    # 1. compute slope dataframe
    df_slope = compute_local_slopes(
        df_raw,
        time_col=time_col,
        signal_col=signal_col,
        window_size_s=100.0
    )

    # 2. detect end-of-transformation region (plateau)
    t_end = detect_t_end_from_slope(df_slope, time_col, slope_col="slope")
    print(f"Detected end of transformation at t = {t_end:.1f} s")

    # Work only up to t_end
    df_eff = df_slope[df_slope[time_col] <= t_end].copy()

    # 3. detect t1 (nucleation onset)
    t1_raw, t1_report = detect_t1_first_onset(
        df_eff,
        time_col=time_col,
        signal_col=signal_col,
        slope_col="slope",
        baseline_frac=0.10,
        k_sigma_primary=4.0,
        confirm_time_s=60.0,
        drop_tolerance_abs=0.5,
        stability_delay_s=100.0,
    )

    # (This is what we want to *show* as t1)
    t1 = t1_report

    # 4. detect all bulk growth bursts after t1
    burst_times = []
    if t1 is not None:
        burst_times = detect_growth_bursts(
            df_eff,
            time_col=time_col,
            slope_col="slope",
            t_start=t1 + 300.0,
            t_end=t_end,
            min_separation_s=180.0,
            prominence_fraction=0.8,
            recovery_fraction=0.2,
            recovery_time_s=200.0,
            backtrack_s=180.0,
            backtrack_fraction=0.5
        )
       
    # 5. Decide how many onsets to report based on --onsets
    # Always include t1 first (if found)
    report_times = []
    report_labels = []
    if t1 is not None:
        report_times.append(t1)
        report_labels.append("t₁")

    # Add bursts in chronological order
    # onsets arg is "how many total", so subtract 1 for t1
    if args.onsets <= 0:
        # auto mode: include all bursts
        max_extra = len(burst_times)
    else:
        max_extra = max(args.onsets - 1, 0)

    for bi, bt in enumerate(burst_times):
        if bi >= max_extra:
            break
        report_times.append(bt)
        # label them t2, t3, t4...
        report_labels.append(f"t{bi+2}")

    # 6. Map each reported time to temperature
    report_temps = [time_to_temp(df_raw, time_col, temp_col, tx) for tx in report_times]

    # 7. Print summary to console
    for lbl, tx, Tx in zip(report_labels, report_times, report_temps):
        if tx is not None:
            if Tx is not None:
                print(f"{lbl}: {tx:.1f} s @ {Tx:.1f} °C")
            else:
                print(f"{lbl}: {tx:.1f} s @ (no temp)")

    # 8. Plot data and onsets
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(df_raw[time_col], df_raw[signal_col], ".", color="blue", lw=0.1, alpha=0.25, label=y_label)

    # vertical lines for each reported onset
    y_top = ax.get_ylim()[1]
    for lbl, tx, Tx in zip(report_labels, report_times, report_temps):
        if tx is None:
            continue
        ax.axvline(tx, color="red", linestyle="--", alpha=0.8)
        if Tx is not None:
            txt = f"{lbl}: {tx:.0f} s\n{Tx:.1f} °C"
        else:
            txt = f"{lbl}: {tx:.0f} s"
        ax.text(
            tx,
            y_top,
            txt,
            rotation=90,
            va="top",
            ha="left",
            fontsize=8,
            color="red"
        )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel(y_label)
    ax.set_title("Unified Slope-Based Onset Detection (ver4)")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(outfig, dpi=300)
    print(f"Saved figure: {outfig}")
    plt.show()

