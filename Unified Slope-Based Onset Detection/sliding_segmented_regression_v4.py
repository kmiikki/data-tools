#!/usr/bin/env python3
"""
sliding_segmented_regression.ver4_onsets.py

Unified slope-based onset detector (multi-onset version).

- Computes local slope dS/dt from the chosen signal (MGI/BW or Total Counts)
- Finds nucleation onset t1 using ignition + stabilization delay
- Finds significant bulk-growth bursts after t1 by scanning slope(t)
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
    Pi-specific logic: if running in /opt/tools, data is assumed to be
    located in ~/rgb/analysis/. If running from Desktop/analysis or
    elsewhere, use relative paths.
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

    if len(t) < 2:
        out = df.copy()
        out["slope"] = np.full(len(df), np.nan)
        return out

    # robust median dt
    dt_vals = np.diff(t)
    dt_vals = dt_vals[np.isfinite(dt_vals)]
    if dt_vals.size == 0:
        out = df.copy()
        out["slope"] = np.full(len(df), np.nan)
        return out

    dt = np.median(dt_vals)
    if not np.isfinite(dt) or dt <= 0:
        out = df.copy()
        out["slope"] = np.full(len(df), np.nan)
        return out

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
                            slope_threshold: float = 0.002,
                            min_plateau_fraction: float = 0.2) -> float:
    """
    Find earliest time after which |slope| < slope_threshold for the remainder.
    If no clear plateau is found, return the last timestamp.

    min_plateau_fraction:
        Ignore plateau candidates that occur too early in the experiment
        (e.g. first 20% of points).
    """
    t = df[time_col].to_numpy()
    s = df[slope_col].to_numpy()

    if len(t) == 0:
        raise ValueError("Empty dataframe in detect_t_end_from_slope")

    start_index = int(len(t) * min_plateau_fraction)

    for i in range(start_index, len(t)):
        if np.all(np.abs(s[i:]) < slope_threshold):
            return float(t[i])

    # fallback: no clear plateau → last time
    return float(t[-1])


def detect_t1_first_onset(df: pd.DataFrame,
                          time_col: str,
                          signal_col: str,
                          slope_col: str,
                          baseline_frac: float = 0.10,
                          k_sigma_primary: float = 4.0,
                          confirm_time_s: float = 60.0,
                          drop_tolerance_abs: float = 0.5,
                          stability_delay_s: float = 60.0,
                          k_sigma_amp: float = 3.0):
    """
    Two-stage nucleation onset detector.

    Stage 1:
        Find a reliable growth wave using:
            - slope >= median_baseline_slope + k_sigma_primary * baseline_slope_std
            - AND signal level >= mean_baseline_level + k_sigma_amp * baseline_level_std
        plus a confirmation window where the signal does not fall back down
        by more than drop_tolerance_abs.

        This gives an index i_slope where we are definitely in the growth wave.

    Stage 2:
        Search BETWEEN the end of the baseline region and i_slope for the
        earliest time index j where:
            - signal >= amp_threshold
            - and the signal stays stable (no drop > drop_tolerance_abs)
        within the same confirm_time_s window.

        This j defines t1_raw (earliest stable amplitude crossing
        belonging to the growth wave).

    Stage 3:
        Report t1_report = t1_raw + stability_delay_s.

    Returns
    -------
    (t1_raw, t1_report) or (None, None)
    """

    t = df[time_col].to_numpy()
    s = df[slope_col].to_numpy()
    y = df[signal_col].to_numpy()

    if len(t) < 2:
        return None, None

    # --- Baseline region (first baseline_frac fraction of data) ---
    n_base = max(5, int(len(t) * baseline_frac))
    n_base = min(n_base, len(t))

    slope_base = s[:n_base]
    slope_base_finite = slope_base[np.isfinite(slope_base)]
    if slope_base_finite.size == 0:
        return None, None

    # slope baseline
    slope_mean = np.median(slope_base_finite)
    slope_std = np.std(slope_base_finite)
    slope_thresh = slope_mean + k_sigma_primary * slope_std

    # amplitude baseline
    y_base = y[:n_base]
    y_base_finite = y_base[np.isfinite(y_base)]
    if y_base_finite.size == 0:
        return None, None

    y_mean = np.mean(y_base_finite)
    y_std = np.std(y_base_finite)
    amp_thresh = y_mean + k_sigma_amp * y_std

    # --- dt & confirm_n ---
    dt_vals = np.diff(t)
    dt_vals = dt_vals[np.isfinite(dt_vals)]
    if dt_vals.size == 0:
        return None, None

    dt = np.median(dt_vals)
    if not np.isfinite(dt) or dt <= 0:
        return None, None

    confirm_n = max(2, int(confirm_time_s / dt))

    first_index_to_check = n_base

    # ======================================================
    # 1) Find a reliable growth wave (slope + amplitude)
    # ======================================================
    i_slope = None
    for i in range(first_index_to_check, len(t)):
        if not np.isfinite(s[i]):
            continue
        if s[i] < slope_thresh:
            continue
        if y[i] < amp_thresh:
            continue

        j_end = min(len(t), i + confirm_n)
        if j_end - i < confirm_n:
            break

        y_future = y[i:j_end]
        # signal must not collapse back down
        if np.all(y_future >= (y[i] - drop_tolerance_abs)):
            i_slope = i
            break

    if i_slope is None:
        return None, None

    # ======================================================
    # 2) Find earliest amplitude onset that belongs to this wave
    # ======================================================
    i_amp = i_slope  # fallback

    for j in range(first_index_to_check, i_slope + 1):
        if y[j] < amp_thresh:
            continue

        j_end = min(len(t), j + confirm_n)
        if j_end - j < confirm_n:
            break

        y_future = y[j:j_end]
        # no collapse back down
        if np.all(y_future >= (y[j] - drop_tolerance_abs)):
            i_amp = j
            break

    t1_raw = float(t[i_amp])
    t1_report = t1_raw + stability_delay_s

    return t1_raw, t1_report


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
    parser.add_argument(
        "--t2-delay",
        type=float,
        default=300.0,
        help="Delay after t₁ before searching for bulk growth bursts (seconds)."
    )
    parser.add_argument(
        "--t2-prom",
        type=float,
        default=0.8,
        help="Relative prominence threshold for growth bursts "
             "(fraction of max slope, e.g. 0.8)."
    )
    parser.add_argument(
        "--t2-back",
        type=float,
        default=0.3,
        help="Relative backtrack slope level for growth bursts "
             "(fraction of max slope, e.g. 0.3)."
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Step 1–2: ensure merged_rgb_counts_time.csv exists
    # ------------------------------------------------------------------
    time_infile_name = "merged_rgb_counts_time.csv"
    raw_infile_name = "merged_rgb_counts.csv"
    outfig_name = "onset_analysis_v4.png"

    # Resolve paths according to Pi logic
    time_infile, outfig = resolve_paths(time_infile_name, outfig_name)
    raw_infile, _ = resolve_paths(raw_infile_name, outfig_name)

    if os.path.exists(time_infile):
        # 1. time-resolved file already exists
        infile = time_infile
    else:
        # 2. create merged_rgb_counts_time.csv from raw data
        if not os.path.exists(raw_infile):
            raise FileNotFoundError(
                f"Could not find '{time_infile_name}' or '{raw_infile_name}'.\n"
                f"Tried paths:\n  {time_infile}\n  {raw_infile}"
            )

        df_base = pd.read_csv(raw_infile, encoding="utf-8")

        if "Timestamp" not in df_base.columns:
            raise KeyError(
                "Input file does not contain required 'Timestamp' column "
                f"(columns are: {list(df_base.columns)})"
            )

        # Timestamps as integers → normalized time in seconds
        ts = pd.to_numeric(df_base["Timestamp"], errors="raise")
        t0 = ts.iloc[0]
        time_s = ts - t0  # e.g. 0, 2, 4, 6, ... if timestamps are 0, 2, 4, 6, ...

        # Insert 'Time (s)' directly after 'Timestamp'
        ts_idx = df_base.columns.get_loc("Timestamp")
        df_base.insert(ts_idx + 1, "Time (s)", time_s)

        # Save new file and continue using it
        df_base.to_csv(time_infile, index=False, encoding="utf-8")
        print(f"Created '{time_infile_name}' from '{raw_infile_name}'")
        infile = time_infile

    # ------------------------------------------------------------------
    # 3. From here on the program proceeds normally
    # ------------------------------------------------------------------

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
    try:
        t_end = detect_t_end_from_slope(
            df_slope,
            time_col,
            slope_col="slope",
            slope_threshold=0.002,
            min_plateau_fraction=0.2,
        )
    except ValueError:
        t_end = float(df_raw[time_col].max())

    # If t_end is practically at the beginning, fall back to full data range
    t_min = float(df_raw[time_col].min())
    t_max = float(df_raw[time_col].max())
    if t_end <= t_min and t_max > t_min:
        print(
            f"Warning: plateau detected at first time point (t = {t_end:.1f} s). "
            "Using full data range instead."
        )
        t_end = t_max

    print(f"Detected end of transformation at t = {t_end:.1f} s")

    # Work only up to t_end
    df_eff = df_slope[df_slope[time_col] <= t_end].copy()
    if len(df_eff) < 2:
        # safety: if filtering removed almost everything, fall back to full data
        df_eff = df_slope.copy()
        t_end = t_max

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
        stability_delay_s=10.0,
    )

    # This is the onset we want to show as "t1"
    t1 = t1_report

    # 4. detect bulk growth bursts after t1 (slope-based t₂, t₃, ...)
    burst_times: list[float] = []
    if t1 is not None:
        burst_times = detect_growth_bursts(
            df_eff,
            time_col=time_col,
            slope_col="slope",
            t_start=t1 + args.t2_delay,        # default: 300 s after t1
            t_end=t_end,
            min_separation_s=180.0,
            prominence_fraction=args.t2_prom,  # default: 0.8
            recovery_fraction=0.2,
            recovery_time_s=200.0,
            backtrack_s=180.0,
            backtrack_fraction=args.t2_back,   # default: 0.3
        )

    # 5. Decide how many onsets to report based on --onsets
    report_times: list[float] = []
    report_labels: list[str] = []

    # Always include t1 if it was found
    if t1 is not None:
        report_times.append(t1)
        report_labels.append("t₁")

    # Add bursts (t₂, t₃, ...) in chronological order
    if args.onsets <= 0:
        # auto mode: include all bursts
        max_extra = len(burst_times)
    else:
        max_extra = max(args.onsets - 1, 0)

    for bi, bt in enumerate(burst_times):
        if bi >= max_extra:
            break
        report_times.append(bt)
        report_labels.append(f"t{bi+2}")

    # 6. Map each reported time to temperature
    report_temps = [time_to_temp(df_raw, time_col, temp_col, tx)
                    for tx in report_times]

    # 7. Print summary to console
    for lbl, tx, Tx in zip(report_labels, report_times, report_temps):
        if tx is not None:
            if Tx is not None:
                print(f"{lbl}: {tx:.1f} s @ {Tx:.1f} °C")
            else:
                print(f"{lbl}: {tx:.1f} s @ (no temp)")

    # 8. Plot data and onsets
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(
        df_raw[time_col],
        df_raw[signal_col],
        ".",
        color="blue",
        lw=0.1,
        alpha=0.25,
        label=y_label,
    )

    # Vertical lines for each reported onset
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
            color="red",
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
