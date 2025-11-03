# Unified Slope-Based Onset Detection (ver4)

**Purpose**  
This Python tool automatically detects *onset points* in colorimetric time-series data (MGI / BW or Total Counts), corresponding to **nucleation** and **bulk growth** events during crystallization or similar processes.  
It uses a robust, slope-based algorithm that operates directly on CSV data exported from the RGB analysis pipeline.

---

## 🔍 Method Overview

The algorithm works in the **time domain** and analyzes the **local slope** (first derivative) of the chosen signal.

### Input data

The script expects either:

- `merged_rgb_counts_time.csv` with a time column already present, or  
- `merged_rgb_counts.csv` with at least:
  - `Timestamp` (integer, increasing)
  - `Tr (°C)` (thermocouple / reactor temperature)
  - `BW` (MGI) and/or `Total Counts`

If `merged_rgb_counts_time.csv` is missing but `merged_rgb_counts.csv` is found, the script automatically creates `merged_rgb_counts_time.csv` by:

1. Reading `Timestamp`
2. Subtracting the first timestamp → `Time (s)` (0, 2, 4, … style sequence)
3. Inserting `Time (s)` immediately after `Timestamp`
4. Saving the result as `merged_rgb_counts_time.csv`

From that point on, analysis always uses `merged_rgb_counts_time.csv`.

### Processing steps

| Step | Operation | Description |
|------|-----------|-------------|
| **1** | **Load CSV** | Uses `Time (s)` (or auto-detected equivalent), `Tr (°C)`, and chosen signal (`BW` or `Total Counts`). |
| **2** | **Preprocessing** | Missing values in the chosen signal are linearly interpolated. |
| **3** | **Local slope** | For each point, a local linear regression in a ±100 s window provides `slope = d(signal)/dt`. |
| **4** | **Plateau / end time (tₑₙd)** | Find earliest time where \|slope\| stays below a small threshold for the rest of the data; if not found, use the last time point. |
| **5** | **First onset (t₁)** | Nucleation onset from a combination of baseline **slope** and **amplitude** thresholds, plus a confirmation and stabilization delay. |
| **6** | **Later onsets (t₂, t₃, …)** | Detected as strong, well-separated growth bursts in the slope signal after t₁. |
| **7** | **Reporting & plotting** | Selected onsets are printed and drawn as vertical red dashed lines with time and temperature annotations. |

---

## 🧠 t₁: Nucleation Onset Logic

The first onset `t₁` is detected in three stages:

1. **Baseline estimation (early fraction of data)**  
   - Use the first ~10% of the time range as baseline.  
   - Compute median and standard deviation of the slope (`slope_mean`, `slope_std`).  
   - Compute mean and standard deviation of the signal level (`y_mean`, `y_std`).

2. **Growth-wave detection (slope + amplitude)**  
   - Scan forward from the end of the baseline region.  
   - Find the first index where **both** conditions hold:
     - `slope >= slope_mean + k_sigma_primary * slope_std`  
     - `signal >= y_mean + k_sigma_amp * y_std`  
   - Require that, within a confirmation window (e.g. 60 s), the signal does **not** fall back down more than `drop_tolerance_abs`.  
   - This identifies a reliable **growth wave**.

3. **Earliest stable amplitude crossing**  
   - Between the baseline end and the growth-wave index, search for the earliest time where the amplitude threshold is crossed and remains stable.  
   - This time is `t₁_raw`.  
   - The reported t₁ is `t₁ = t₁_raw + stability_delay_s` (a small delay, e.g. 10 s, to match what a human would mark as visually stable onset).

---

## 🧠 t₂, t₃, …: Bulk Growth Bursts

Additional onsets (t₂, t₃, …) are detected from the **slope signal** after t₁:

1. Restrict the analysis to times between `t₁ + t2_delay` and `t_end`.
2. Compute the maximum slope within this interval (`s_max`).
3. Any local slope peak is considered a **burst candidate** if:
   - slope at the peak ≥ `prominence_fraction * s_max`  
   - it is a local maximum (≥ neighbors)
4. For each candidate:
   - **Backtrack** in time until the slope falls below `backtrack_fraction * s_max` or a maximum backtrack window (e.g. 180 s) is reached.  
     This defines the **burst start time**, which is reported as t₂, t₃, … in chronological order.
   - Enforce **recovery** and **minimum separation** rules so that multiple bursts are distinct and not just noise.

---

## ⚙️ Command-Line Usage

The script is usually installed as an executable (e.g. in `~/bin`) and can be called directly:

```bash
sliding_segmented_regression_v4.py [OPTIONS]
```
## Options

| Option       | Type  | Description                                                                                             | Default |
| ------------ | ----- | ------------------------------------------------------------------------------------------------------- | ------- |
| `--signal`   | str   | Which signal to analyze: `BW` (MGI / brightness) or `Total` (`Total Counts`).                           | `BW`    |
| `--onsets`   | int   | How many onsets to report. `0` = automatic (t₁ + all detected bursts), `1` = t₁ only, `2` = t₁+t₂, etc. | `0`     |
| `--t2-delay` | float | Time offset after t₁ before searching for slope bursts t₂, t₃, … (seconds).                             | `300.0` |
| `--t2-prom`  | float | Relative prominence threshold for bursts (fraction of max slope, e.g. `0.8`).                           | `0.8`   |
| `--t2-back`  | float | Relative backtrack level for bursts (fraction of max slope, e.g. `0.3`).                                | `0.3`   |

The input and output filenames are fixed by convention:

Input: merged_rgb_counts_time.csv (or it is derived from merged_rgb_counts.csv)
Output figure: onset_analysis_v4.png
Typical usage:

```bash
sliding_segmented_regression_v4.py --signal BW --onsets 2
```

🖼️ Example Output

Console:

```text
Created 'merged_rgb_counts_time.csv' from 'merged_rgb_counts.csv'
Detected end of transformation at t = 13560.0 s
t₁: 2220.0 s @ 54.0 °C
t₂: 4894.0 s @ 46.6 °C
Saved figure: onset_analysis_v4.png
```

Figure:

- Blue dots: MGI (or Total Counts) vs. time
- Red dashed vertical lines: detected onsets (t₁, t₂, …)
- Each line annotated with time and corresponding temperature from Tr (°C).

🧩 Dependencies

Install required libraries, for example with pip:

```bash
pip install numpy pandas matplotlib
```

The script itself does not currently depend on SciPy or other non-standard packages.

🧪 Tested Environment

Python 3.12 (conda / micromamba)
- NumPy 2.x
- Pandas 2.x
- Matplotlib 3.x
- Tested on:
- Ubuntu 24.04 (Niobium)
- Raspberry Pi OS Bookworm (Raspberry Pi 5)

🧬 Citation

If this method or code is used in your research, please cite as:

Miikki, K. (2025). Unified Slope-Based Onset Detection for Crystallization Processes.
Aalto University, School of Chemical Engineering.

🧑‍💻 Author

Kim Miikki<br>
Research Engineer – Aalto University, CHEM<br>
📧 kim.miikki@aalto.fi
