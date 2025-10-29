# Algorithm: Unified Slope-Based Onset Detection (ver4)

## 1. Goal

The algorithm detects physically meaningful transition points (“onsets”) in crystallization / fouling experiments, where image-derived intensity (MGI) changes over time.

It returns:

* **`t_1`** = nucleation onset
* **`t_2`, `t_3`, …** = subsequent bulk-growth bursts
* (optionally) process plateau / end time

It also annotates each onset with the corresponding process temperature.

---

## 2. Inputs and Definitions

### Required columns in the input CSV

* `Time (s)` — time in seconds from start of experiment
* `Tr (°C)` — temperature
* `BW` (MGI) — mean grayscale intensity of the region of interest

  * This is the primary signal for crystallization (more crystals → higher MGI)

Optionally, the script can also work with `Total Counts` instead of `BW` via `--signal Total`.

### Notation

Let `y(t) = MGI(t)` be the intensity signal.

Let `s(t) = dy/dt` be its local slope (rate of change).

---

## 3. Preprocessing

### 3.1 Interpolation and smoothing

* Missing values in the selected signal (`BW` / `MGI` or `Total Counts`) are linearly interpolated.
* Local slopes are computed using a sliding linear regression window (default ~100 s window).
  For each time index ( i ), a line is fitted to ((t, y)) in a symmetric time window around ( t_i ), and the slope of that fit is taken as ( s(t_i) ).

This produces:

* a smoothed version of the trend, and
* a robust estimate of slope even in noisy segments.

---

## 4. Plateau detection (end-of-process)

We define the effective end of transformation `t_end` as the earliest time after which the absolute slope stays below a small threshold for the rest of the experiment:

The end of transformation `t_end` is defined as the earliest time after which the absolute slope stays below a small threshold:

`|s(t)| < ε  for all  t' ≥ t`

Typical `ε` is approximately 0.002–0.005 (MGI units / s).

If no such point exists, `t_end` defaults to the last timestamp.

All subsequent analysis is restricted to the data range:

`[0, t_end]`


---

## 5. Detecting the first onset (`t_1`: nucleation onset)

This is the moment when the system stops being quiet and starts to change irreversibly.

### 5.1 Baseline characterization

We take the first fraction of the experiment (e.g. first 10% of samples) as “baseline.”  
From that segment we compute:

- median baseline slope `μ_b`
- standard deviation of slope `σ_b`

These describe the normal drift / noise before nucleation.

### 5.2 Threshold crossing

We scan forward in time and find the earliest point `t*` such that:

1. The slope at that time is significantly above baseline:

   `s(t*) > μ_b + k·σ_b`

   with `k = 4` by default.

2. The increase is not just a spike:  
   after `t*`, the intensity `y(t)` must **not** drop back down by more than a small absolute tolerance  
   (e.g. 0.5 MGI units) within a short confirmation window (e.g. 60 s).

This rejects false positives from momentary flicker or camera noise.

We call this raw time `t1_raw`.

### 5.3 Stabilization delay

For reporting purposes, we shift forward by a small fixed stabilization delay (e.g. +100 s):

`t_1 = t1_raw + Δt_stabilize`

Reason: visually, nucleation is judged after the turbidity is clearly established,  
not at the very first detectable blip. The stabilized time `t_1` matches expert labeling  
(what a human would call “the onset”).

We also record the temperature at `t_1` by looking up the nearest `Tr (°C)` value at that timestamp.

---

## 6. Detecting later onsets (`t_2`, `t_3`, …)

These correspond to “bulk growth bursts,” i.e. rapid crystal growth phases after nucleation.

### 6.1 Region of interest

We now analyze only times after `t_1`, with an intentional guard period:

`t ≥ t_1 + 300 s`

Why?  
The immediate aftermath of `t_1` often still belongs to the same physical event.  
We do not want to call that “t_2.”

We also restrict to:

`t ≤ t_end`

### 6.2 Find candidate bursts from slope peaks

Within that interval:

1. Find local maxima of the slope `s(t)`.  
   A point counts as a candidate if:
   - it is a local peak (greater than its neighbors), and  
   - its amplitude is large — for example:

     `s(t_peak) ≥ α·s_max`

     where `s_max` is the maximum slope in the interval and `α ≈ 0.8`.  
     This enforces *prominence*: only major growth surges are considered.

2. For each candidate peak at `t_peak`, walk **backward in time** to locate when that burst truly *started accelerating*.  
   We track back until the slope falls below a fraction of the maximum slope  
   (e.g. `0.5 × s_max`) or until a backtracking time limit (e.g. 180 s) is reached.

   The first timestamp in that high-slope ramp is taken as the **burst onset time**, which we call `t_burst`.

This “backtracking” step is crucial:

- We do not report the time of maximum growth rate.  
- We report when the growth phase actually *ignited*.

That is exactly the physical moment users want to publish  
(e.g. “bulk crystallization started at 4894 s”).

### 6.3 Separation and recovery criteria

We do not want to label tiny fluctuations inside the same burst as multiple separate onsets.
So we enforce two additional rules between successive bursts:

1. **Minimum separation in time**
   Two bursts must be separated by at least `min_separation_s` (e.g. 180 s).
   This prevents near-duplicates like 4894 s and 4910 s.

2. **Recovery / quiescence requirement**
   Between two bursts, the slope must fall back to a low level and stay low long enough.
   Concretely: between the previous burst and the new candidate, slope must spend
   ≥ `recovery_time_s` seconds below a small fraction of max slope (e.g. <20%).
   This encodes: “the system calmed down, then a *new* growth wave started.”

After applying these filters, the remaining burst onset times are labeled `t_2`, `t_3`, ... in chronological order.

Each onset is paired with the temperature `Tr (°C)` at that time.

---

## 7. Practical output

For a typical run, the algorithm prints something like:

```text
Detected end of transformation at t = 13560.0 s
t_1: 2220.0 s @ 54.0 °C
t_2: 4894.0 s @ 46.6 °C
```

