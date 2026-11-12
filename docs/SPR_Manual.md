# Single Pulse Response (SPR) Tab Manual

## Overview

The **SPR Tab** allows you to analyse the "Single Pulse Response" (washout characteristics) of your laser ablation system. By analysing the shape of individual laser pulses, the software can calculate the **Washout Time** required for the signal to decay to a specific baseline level (`10%` or `1%`).

This washout time is a critical input for the **Optimiser**, ensuring that subsequent laser shots are properly accounted for to prevent unintended signal mixing or aliasing.

---

## User Interface

The interface is divided into two main sections: **Controls & Results** (Left) and **Visualisation** (Right).

### 1. Left Panel: Controls & Results

#### **Global Controls**

* **Import SPR File**: Loads a single pulse response data file into the SPR tab (and automatically unloads the temporary raw channels from your active iolite session once read).

#### **Peak Detection**

Configure how the software identifies individual laser pulses in the signal.

* **Select Channel**: Choose the isotope channel to analyse (e.g., `U238`, `Th232`). Usually, a high-abundance isotope is best for characterising washout.
* **Peak Cutoff**:
  * **Auto**: Automatically sets the minimum peak height threshold to `30%` of the maximum data point (after background subtraction if enabled).
  * **Manual**: Uncheck **Auto** to enter a custom value in the spin box, or hover over and drag the dashed **Peak Cutoff** line directly on the SPR plot.
* **Min. Peak Spacing**: The minimum time spacing required between two consecutive peaks (defaults to `900 ms` for a `1 Hz` pulse train). Increasing this prevents detecting multiple points on the same peak as separate events.
* **Auto-Exclude Wide Outliers (Percentile >)**: When checked, automatically excludes unusually wide pulses whose widths exceed the specified percentile (defaults to `90.0 %`).
* **Time Unit**: Toggle between **Seconds (s)** and **Milliseconds (ms)** for peak spacing and washout display.
* **Subtract Background**: Subtracts the baseline gas blank before measuring peak heights and widths. When **Auto Select Region** is enabled in the plot controls, the software automatically detects whether at least `5 seconds` of gas blank exists at the start or end of the file, and automatically unchecks **Subtract Background** if no gas blank is present.

#### **SPR Analysis Results**

Displays aggregate statistics for all detected (and included) peaks on the currently selected channel. Both **FW 0.1M (10%)** (width at `10%` of peak maximum) and **FW 0.01M (1%)** (width at `1%` of peak maximum) panels provide:

* **Average** & **RSD**: Mean washout duration across all included peaks and its Relative Standard Deviation (%).
* **Area** & **RSD**: Mean integrated area (`Counts`) under the peaks and its Relative Standard Deviation (%).
* **Max**: The longest washout duration observed among the included peaks on the selected channel.
* **Composite**: The width of the time-aligned **Composite Peak** (average pulse profile).
* **Apply Average / Apply Composite / Apply Max**: Sends the corresponding `10%` or `1%` washout value directly to the **Washout (ms)** input on the **Method Optimiser** tab.

#### **Maximum Observed SPR**

A dedicated summary panel that scans all loaded channels to identify the channel with the longest washout duration.

* Displays the winning channel name alongside its **Average**, **Max**, and **Composite** washout times at both **FW 0.1M (10%)** and **FW 0.01M (1%)**.
* **Apply Average / Apply Composite / Apply Max**: Sends the `1%` washout time of the winning channel (`Average`, `Composite`, or `Max`) to the **Method Optimiser** tab for conservative multi-element tuning.
* **Show All Channel Stats**: Switches the bottom table from individual peaks on the current channel to a multi-channel summary table comparing `Avg`, `Max`, `Comp`, `RSD`, and `Area` at `10%` and `1%` across every loaded channel.

### 2. Right Panel: Visualisation

#### **Plot Controls**

To preserve display height on compact screens, plot controls are arranged in two rows above the plot:

* **Row 1**:
  * **Theme**: Select `Auto` (matches system theme), `Dark`, or `Light`.
  * **Pan / Zoom Y**: Enables vertical zooming and panning with the mouse.
  * **Auto-Rescale Y**: Defaults to Checked. Automatically scales the Y-axis to fit the visible signal whenever parameters or channels change.
  * **Hide Summary Stats**: Hides or shows the floating statistics overlay box on the main SPR plot.
* **Row 2** *(visible when **Subtract Background** is checked)*:
  * **Auto Select Region**: Automatically locates a `>= 5 s` baseline region at the start or end of the recording.
  * **Background (s)**: Start and end time spin boxes for the background region (can also be adjusted on the plot via **Shift + Left-Click & Drag**).

#### **SPR Plot (Top Left)**

* Displays the raw signal intensity over time along with the shaded **Background** region (when background subtraction is active) and the horizontal dashed **Peak Cutoff** line.
* **Markers**:
  * **Orange X**: Detected and included peaks.
  * **Grey X**: Excluded peaks.
* **Interactions**:
  * **Click Peak Marker or Label**: Toggles whether an individual pulse is excluded or included in the washout statistics.
  * **Drag Peak Cutoff Line**: Hover over the dashed horizontal **Peak Cutoff** line until the vertical resize cursor appears, then **Left-Click & Drag** up or down to adjust the threshold interactively.
  * **Shift + Click & Drag Background Edge**: Adjusts the start or end boundary of the blue background region.

#### **Composite Plot (Top Right)**

* Shows the normalised "Average Peak Shape" by time-aligning and stacking all included pulses, along with horizontal markers for the `10%` (`FW 0.1M`) and `1%` (`FW 0.01M`) washout widths.
* Useful for visualising tailing behaviour and verifying whether the `10%` or `1%` widths are appropriate.

#### **Mouse Navigation (Plots)**

* **Zoom**: Use the **Mouse Wheel** over any plot to zoom in/out on the X-axis (and Y-axis if **Pan / Zoom Y** is enabled).
* **Pan**: **Left-Click & Drag** on the plot canvas to pan the view.
* **Double-Click**: Resets the plot to show the full data range.
* **Legend Click**: Single-click legend items to show/hide traces.

#### **Results Table (Bottom)**

* **Single-Channel View** *(default)*: Lists every detected peak (`Peak`, `Time (s)`, `FW0.1M (10 %)`, `FW0.1M (10 %) Area`, `FW0.01M (1 %)`, `FW0.01M (1 %) Area`, `Max Intensity`). Rows for excluded peaks are shown with strikethrough text.
* **All-Channel View** *(when **Show All Channel Stats** is checked)*: Lists summary metrics for every channel in the file.

---

## Workflow Guide

### Step 1: Load Data

Ensure you have a recording of single laser pulses (e.g., a line or spot scan run at `1 Hz` repetition rate).

1. Open the **iolite Optimiser**.
2. On the **Determine SPR** tab, click **Import SPR File** and select your data file.

### Step 2: Select Channel

1. In the **Peak Detection** group, select a representative channel (e.g., `U238` or a major matrix element), or review the **Maximum Observed SPR** panel to see which channel has the longest washout.

### Step 3: Optimise Detection

1. Observe the **SPR Plot**. Are all pulses marked with an orange `X`?
2. If noise is being detected as peaks:
   * Increase the **Peak Cutoff** (drag the dashed **Peak Cutoff** line up on the plot, or uncheck **Auto** and increase the value).
   * Increase **Min. Peak Spacing**.
3. If valid peaks are missed:
   * Drag the dashed **Peak Cutoff** line down (or decrease the **Peak Cutoff** spin box value).

### Step 4: Refine Data

1. Inspect the **SPR Plot** and **Results Table** for outliers (e.g., double pulses or noise spikes).
2. **Exclude Outliers**: Click on a peak's marker or number label on the plot (or enable **Auto-Exclude Wide Outliers**) to remove it from the calculation.
   * The statistics in the left panel and the **Composite Peak** plot update immediately.
   * Excluded peaks turn grey on the plot and display with strikethrough text in the table.

### Step 5: Apply to Optimiser

1. Decide whether to apply a `10%` (**FW 0.1M**), `1%` (**FW 0.01M**), or **Maximum Observed SPR** washout time, and whether to use **Average**, **Composite**, or **Max**:
   * **Average**: Mean width across individual pulses; good for clean, high-intensity channels.
   * **Composite**: Width of the stacked average pulse profile; best for noisy signals or trace elements.
   * **Max**: Longest individual pulse observed; most conservative choice to guarantee zero overlap.
2. Click the corresponding **Apply Average**, **Apply Composite**, or **Apply Max** button.
3. The selected washout time (converted automatically to milliseconds) is applied directly to the **Washout (ms)** parameter on the **Method Optimiser** tab.

---

## Troubleshooting & Logging

### Date/Time Stamp Format Mismatch on CSV Imports

When importing CSV files, if the timestamp format does not match the active settings in iolite, the import will fail.

If this happens, the plugin intercepts the failure and displays an alert showing the specific failure details from the iolite application log, along with suggested date/time format strings (e.g., `yyyy MM dd hh mm ss` or `dd MM yyyy hh mm ss`) based on the digits in your file.

### Other Issues

* **No Peaks Detected**:
  * Check if the correct channel is selected.
  * Lower the **Peak Cutoff** by dragging the dashed threshold line down.
  * Ensure your data contains separated pulses at or below the **Min. Peak Spacing** interval.
* **Background Subtraction Automatically Unchecked**:
  * If a recording starts and ends immediately with laser pulses (fewer than `5 seconds` of continuous gas blank at both the start and end of the file), **Subtract Background** is automatically unchecked so pulse intensities are not subtracted as baseline. You can still re-check it and set a custom range manually if needed.

---

## Methodology

The **SPR Tab** characterises system washout by analysing the temporal profile of individual laser pulses:

1. **Detection**: Peaks are identified above the **Peak Cutoff** height separated by at least **Min. Peak Spacing**.
2. **Width Calculation**: Linear interpolation between adjacent samples determines the exact time span between the rising and falling crossings at `10%` (`FW 0.1M`) and `1%` (`FW 0.01M`) of peak maximum.
3. **Composite Peak**: An average pulse profile is constructed by normalising, time-aligning, and stacking all included pulses using baseline-split dynamic padding between consecutive shots.
