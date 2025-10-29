# Unified Slope-Based Onset Detection (ver4)

**Purpose:**  
This Python tool automatically detects *onset points* in colorimetric (MGI–time) data, corresponding to **nucleation** and **bulk growth** events during crystallization or similar processes.  
It provides a robust, slope-based algorithm that operates directly on time-series input data.

---

## 🔍 Method Overview

The algorithm analyzes the **first derivative** of the measured intensity (MGI) signal:

| Step | Operation | Description |
|------|------------|-------------|
| **1** | **Input CSV** | Columns: `Time (s)`, `Tr (°C)`, `MGI` |
| **2** | **Preprocessing** | Linear interpolation, smoothing, slope computation `d(MGI)/dt` |
| **3** | **First Onset (t₁)** | Detected when slope exceeds baseline + σ threshold for ≥N samples |
| **4** | **Additional Onsets (t₂, t₃…)** | Local slope maxima separated by drift and prominence rules |
| **5** | **Plateau Detection** | Phase end detected when slope ≈ 0 for ≥300 s |
| **6** | **Visualization** | MGI curve plotted with red dashed onset lines and temperature annotations |

---

## ⚙️ Command-Line Usage

```bash
python3 sliding_segmented_regression_v4.py [--onsets N] [--signal BW|Total]
```
## Arguments

| Option     | Description                                                          | Default                      |
| ---------- | -------------------------------------------------------------------- | ---------------------------- |
| `--onsets` | Number of onsets to detect (1–3). If not given, automatic detection. | auto                         |
| `--signal` | Select analyzed signal: `BW` (MGI brightness) or `Total Counts`.     | `BW`                         |
| `--infile` | Input CSV file containing MGI data.                                  | `merged_rgb_counts_time.csv` |
| `--outfig` | Output figure file (PNG).                                            | `onset_detection.png`        |

```bash
python3 sliding_segmented_regression_v4.py --onsets 2 --signal BW
```

📊 Example Result

Detected end of transformation at t = 13560.0 s

| Onset  | Time (s) | Temperature (°C) | Interpretation     |
| ------ | -------- | ---------------- | ------------------ |
| **t₁** | 2220.0   | 54.0             | Nucleation begins  |
| **t₂** | 4894.0   | 46.6             | Bulk growth starts |

🧠 Algorithm Summary

Compute derivative of MGI over time.

Identify first sustained slope rise above baseline noise (t₁).

Locate additional slope bursts (t₂, t₃…) using prominence and persistence tests.

Define plateau (end) when slope remains near zero.

Report onset times and corresponding temperatures.

This method combines the interpretability of classical signal processing with the reproducibility of programmatic detection.

🧩 Dependencies

Install the required libraries using:

```bash
pip install numpy pandas matplotlib scipy
```

Optionally:

```bash
pip install piecewise-regression reportlab
```

🧪 Tested Environment

Python 3.12 (conda / micromamba)

NumPy 2.x, Pandas 2.x, SciPy 1.14

Tested under Ubuntu 24.04 (Niobium) and Raspberry Pi 5 (Bookworm)

🧬 Citation

If this method or code is used in your research, please cite as:

Miikki, K. (2025). Unified Slope-Based Onset Detection for Crystallization Processes.
Aalto University, School of Chemical Engineering.

🧑‍💻 Author

Kim Miikki<br>
Research Engineer – Aalto University, CHEM<br>
📧 kim.miikki@aalto.fi
