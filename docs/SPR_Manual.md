# Single Pulse Response (SPR) Tab Manual

## Overview

The **SPR Tab** allows you to analyze the "Single Pulse Response" (washout characteristics) of your laser ablation system. By analyzing the shape of individual laser pulses, the software can calculate the **Wait Time** (washout time) required for the signal to decay to a specific baseline level (10% or 1%).

This specific washout time is a critical input for the **Optimiser**, ensuring that subsequent laser shots do not overlap unintentionally or that necessary delays are added between spots to prevent signal mixing.

---

## User Interface

The interface is divided into two main sections: **Controls & Results** (Left) and **Visualisation** (Right).

### 1. Left Panel: Controls & Results

#### **Global Controls**

* **Reload Data from iolite**: Refreshes the input data from the main iolite interface. Use this if you have changed selections or processing in iolite while the plugin is open.

#### **Peak Detection**

Configure how the software identifies individual laser pulses in the signal.

* **Select Channel**: Choose the isotope channel to analyze (e.g., `U238`, `Th232`). Usually, a high-abundance isotope is best for characterizing washout.
* **Peak Cutoff (Prominence)**:
  * **Auto**: Automatically estimates a threshold to identify peaks above the baseline noise.
  * **Manual**: Uncheck "Auto" to set a specific prominence value. Higher values filter out noise but may miss smaller peaks.
* **Min. Distance**: The minimum number of data points required between two peaks. Increasing this helps prevent detecting multiple points on the same peak as separate events.
* **Auto-Rescale Y**: Defaults to Checked. When enabled, the plot always scales to show the full pulse intensity.
* **Time Unit**: Toggle between **Seconds (s)** and **Milliseconds (ms)** for display.

#### **SPR Analysis Results**

Displays aggregate statistics for all detected (and included) peaks.

* **FW 0.1M (10%)**: Width of the peak at 10% of its maximum height. This is the standard definition for washout time in many applications.
  * **Average**: Mean duration across all peaks.
  * **RSD**: Relative Standard Deviation (%) – measures consistency.
  * **Max**: The longest duration observed (conservative estimate).
  * **Area**: Mean integrated area (counts) under the peaks.
* **FW 0.01M (1%)**: Width at 1% of maximum height (stricter washout criterion).

**Actions:**

* **Apply Average to Optimiser**: Sends the *Average* 10% washout time to the Optimiser tab for the selected channel.
* **Apply Max to Optimiser**: Sends the *Maximum* 10% washout time (more conservative) to the Optimiser tab.
* **Apply Composite to Optimiser**: Sends the calculated width of the **Composite Peak** (average pulse shape) to the Optimiser. This is often the most representative and stable value for washout.

### 2. Right Panel: Visualisation

#### **SPR Plot (Top)**

* Displays the raw signal intensity over time.
* **Markers**:
  * **Orange X**: Detected peaks.
  * **Grey X**: Excluded peaks.
* **Interaction**: Click on a peak marker or its label to **Exclude/Include** it from the statistics. Excluded peaks are grayed out and removed from calculations.

#### **Composite Plot (Bottom Right - if enabled)**

* Shows the "Average Peak Shape" by stacking all valid pulses.
* Useful for visualising the tailing behavior and verifying if the 10% or 1% widths are appropriate.

#### **Mouse Interactions (Plots)**

* **Zoom**: Use the **Mouse Wheel** over any plot to zoom in/out on the X-axis. If "Pan/Zoom Y" is enabled, it will also zoom the Y-axis.
* **Pan**: **Left-Click & Drag** a plot to move the view.
* **Double-Click**: Instantly resets the plot to show the full data range.
* **Auto-Rescale Y**: When checked, the plot automatically snaps to show the full pulse intensity whenever values change. Uncheck this to manually zoom into specific features (like the peak tail).

#### **Results Table (Bottom)**

* Lists every detected peak individually.
* Columns: `Peak ID`, `Time`, `10% Width`, `10% Area`, `1% Width`, `1% Area`, `Max Intensity`.
* **Interaction**: Rows corresponding to excluded peaks are shown with strikethrough text.

---

## Workflow Guide

### Step 1: Load Data

Ensure you have a selection in iolite that contains single laser pulses (e.g., a "Test" or "Tune" line with low repetition rate, like 1 Hz).

1. Open the **Advanced Spot Optimiser**.
2. Click **Reload Data from iolite** if your data is not visible.

### Step 2: Select Channel

1. Go to the **SPR** tab.
2. In the **Peak Detection** group, select a representative channel (e.g., `U238` or a matrix element).

### Step 3: Optimize Detection

1. Observe the **SPR Plot**. Are all pulses marked with an orange `X`?
2. If noise is being detected as peaks:
    * Increase the **Peak Cutoff** (uncheck Auto).
    * Increase **Min. Distance**.
3. If peaks are missed:
    * Decrease the **Peak Cutoff**.

### Step 4: Refine Data

1. Look at the **Results Table** or **SPR Plot** for outliers (e.g., double pulses, noise spikes).
2. **Exclude Outliers**: Click on the peak marker in the plot to remove it from the calculation.
    * The statistics in the Left Panel will update instantly.
    * The peak will turn grey in the plot and strikethrough in the table.

### Step 5: Apply to Optimiser

1. Decide whether you want to use the **Average**, **Maximum**, or **Composite** washout time.
    * *Average*: Good for stable systems.
    * *Max*: Safer to ensure zero overlap.
    * *Composite*: Best for noisy signals or when you want the most statistically representative shape.
2. Click **Apply Average/Max/Composite to Optimiser**.
3. A confirmation message will appear. This value is now set as the "Custom" dwell/wait time for this channel in the **Optimiser Tab**.

---

## Troubleshooting

* **No Peaks Detected**:
  * Check if the correct channel is selected.
  * Lower the **Peak Cutoff**.
  * Ensure your data actually contains single pulses.
* **"Error" in Results**:
  * Usually means the baseline is too high or the pulse never drops below 10%/1% within the available window.
  * Try selecting a cleaner channel.
* **Optimiser vs SPR Units**:
  * The Optimiser usually works in **Milliseconds**. the SPR tab handles unit conversion automatically when you click "Apply".

---

## Methodology & References

The **SPR Tab** characterizes system washout by analyzing the temporal profile of individual laser pulses, as described in:

- **Washout Characteristics**:
  - Ulianov, A., et al. (2015). "Single-pulse responses of laser ablation cells." *Journal of Analytical Atomic Spectrometry (JAAS)*, 30, 1297.

### Peak Analysis Process:
1. **Detection**: Peaks are identified using a prominence-based algorithm.
2. **Width Calculation**: The software determines the time taken for the signal to decay from its maximum to 10% (FW 0.1M) and 1% (FW 0.01M) of its peak height.
3. **Composite Peak**: An "Average Peak Shape" is generated by time-aligning and stacking all valid pulses. This provides a high-signal-to-noise profile for precise washout determination, even on trace elements.
