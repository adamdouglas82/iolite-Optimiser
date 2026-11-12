# Method Optimiser Tab Manual

## Overview

The **Method Optimiser** tab determines the optimal combination of **Spot Size**, **Laser Repetition Rate**, **Stage Scan Speed**, and **Dwell / Acquisition Times** to achieve your analytical goals. It works backwards from your target signal quality (SNR) and hardware constraints to find the most efficient sampling parameters.

---

## User Interface

### 1. Left Panel: Settings & Configuration

#### **Top Controls**

* **Import File**: Loads an ablation test or PreScan file into the Optimiser (and automatically unloads the temporary raw channels from your active iolite session once read).
* **Settings**: Opens the **Hardware Configuration** dialog.

#### **Hardware Configuration Dialog**

Accessible via the **Settings** button, this dialog establishes the physical limits of your ICP-MS and laser ablation system.

* **ICP-MS Hardware**:
  * **Manufacturer / Model**: Presets minimum dwell times, resolution steps, and instrument architecture (`Quadrupole`, `Sector-Field`, `Multi-Collector`, or `TOF`).
  * **Dwell Times in**: Choose `Auto` (uses `ms` for Quadrupole/Sector-Field/TOF and `s` for Multi-Collector), `ms`, or `s`.
  * **Custom Models**: Select `Custom Model` to manually define:
    * **Type**: `Quadrupole`, `Sector-Field`, `Multi-Collector`, or `TOF`.
    * **Minimum Dwell Time (ms)** & **Dwell Resolution (ms)**: Hardware integration speed and rounding step limits.
    * **Allowed Dwell Times (ms)**: *(Multi-Collector only)* Comma-separated list of valid integration times.
  * **Washout Margin (%)**: *(TOF only)* A percentage signal capture buffer added to the baseline washout time so the full Single Pulse Response peak is captured within the stage speed budget.
* **Laser Hardware**:
  * **Platform / Laser Source**: Presets maximum repetition rates and frequency step resolution.
  * **Cell Type**: Defines the maximum stage translation speed.
  * **Custom Laser**: Select `Custom` (Platform) or `Custom Laser` (Source) to manually configure:
    * **Mode**: `Maximum Rep-Rate` (continuous up to a maximum frequency) or `Discrete` (comma-separated list of allowed frequencies).
    * **Max Rep-Rate (Hz)** & **Rep-Rate Resolution (Hz)**: Upper frequency limit and step size.
    * **Max Speed (µm s⁻¹)**: Physical translation speed limit of the sample stage (`Custom Stage`).
* **Plot Options**:
  * **Limit Visible Plots**: Restricts the initial number of plotted channels (default `5`) to keep the plot legible when loading files with many masses.
* **Utilities & Diagnostics**:
  * **Enable Logging**: *(Default: Off, persisted across restarts)* When checked, diagnostic messages, warnings, and error tracebacks are written to the iolite **Messages** log (including startup events).
  * **GitHub Repo** & **Report Issue**: Opens the project repository or issue tracker in your default web browser.
  * **Check for Updates**: Queries GitHub for the latest release version and offers a direct link to the release notes if a newer version is available.

#### **Input Parameters**

Basic setup from your test or PreScan recording:

* **Initial Spot Size (µm)**: The beam diameter used during the test scan.
* **Washout (ms)**: The time required for a single pulse signal to decay (typically applied directly from the **Determine SPR** tab). Supports fast sub-millisecond cells down to `0.01 ms`.
* **Initial Rep-Rate (Hz)**: The laser repetition rate used during the test scan.
* **Mode**:
  * **Spot**: Static ablation at a single location (configures **Number of Shots**).
  * **Line**: Continuous translation along a transect (configures **Line Length (µm)**).
  * **Imaging**: 2D raster mapping (configures **Image Width (µm)** and **Image Height (µm)**).
* **Number of Analytes**: Total number of isotopes/masses in your proposed analytical method (used by the PreScan calculator and saved persistently).
* **Suggest PreScan Params**: Opens a calculator dialog recommending initial PreScan conditions before acquiring a test line:
  * **Dwell Time per Analyte**: `Washout / Number of Analytes`, snapped to your instrument's minimum dwell time and resolution.
  * **Target Cycle Time**: `Dwell Time * Number of Analytes`.
  * **Suggested Rep-Rate**: Sweeps upward from the oversampling threshold (`2.0 / Washout`) to find the lowest valid laser frequency achieving a steady-state Lockwood RSD `<= 5%`.
  * **Suggested Stage Speed**, **Dosage**, and **Overlap**: Derived from your spot size and cycle time, capped at the cell's maximum stage speed.

#### **Optimisation Parameters**

The analytical targets for the optimisation algorithm:

* **Sync Strategy**:
  * **Adaptive Integer Sync (Auto)**: Automatically sweeps repetition rates and acquisition times to guarantee a perfect integer pulse count per cycle (`0%` rounding error) close to your target.
  * **Strict Pulse Target (Manual)**: Directly calculates acquisition times for your exact pulse target (hardware rounding may produce non-integer pulse counts).
  * **Oversampling (Time-Driven)**: Focuses purely on meeting the washout time and target steady-state RSD constraints without integer pulse snapping.
  * **Combined Integer Sync (Auto)**: Sweeps repetition rates above the oversampling minimum to guarantee both steady-state stability and an integer pulse count per cycle.
* **Sync Target (Pulses)**: *(Synchronised modes)* Target number of laser pulses per mass spectrometer acquisition cycle (`Acq Time`).
* **Sync Target (RSD %)**: *(Oversampling modes)* Target Relative Standard Deviation limit (default `5.0%`) for steady-state signal stability.
* **Prefer Exact Pulses per Acq**: *(Adaptive Integer Sync mode)* When checked, prioritises finding a repetition rate that delivers the exact requested pulse count rather than the closest achievable integer sync.
* **Avoid Gaps**: *(Strict Pulse Target mode)* Ensures `Speed <= Spot Size * Rep-Rate` so adjacent laser shots never leave unablated gaps along the scan line.
* **Min Duty Cycle (%)**: Sets a minimum floor for measurement efficiency (`Total Dwell Time / Total Cycle Time`). If switching overhead dominates the cycle, dwell times are increased to meet this floor.
* **Target SNR (Sigma)**: Desired Signal-to-Noise Ratio for the limiting isotope:
  * **Noise Floor**: Uses the greater of theoretical Poisson counting variance and observed background variance.
  * **Critical Level (Lc)**: Computed using Stapleton's Square Root Transform rule (`95%` confidence threshold for detection, robust even at zero background counts).
  * **Sigma-Equivalent Scaling**: Scales the detection ratio by `2.326` (`z * sqrt(2)`) so a target of `10.0` corresponds to the standard `10-sigma` Limit of Quantification (LoQ).
* **Minimum SNR (Sigma)**: Lower acceptance threshold; channels that cannot reach this SNR are flagged (and set to minimum dwell on sequential instruments to preserve cycle budget).
* **Scale Signal with Rep-Rate**:
  * **Checked (Linear Scaling)**: Assumes signal intensity (CPS) scales linearly with repetition rate (appropriate for bulk geological mounts where firing faster ablates deeper into fresh material).
  * **Unchecked (Fixed Sensitivity)**: Assumes signal intensity depends only on spot size (appropriate for thin sections or bio-imaging where a single shot already ablates the full sample thickness).

#### **Optimised Settings (Results Summary)**

Displays the recommended instrument parameters:

* **Spot Size (µm)**: Recommended beam diameter (defaults to `30 µm` before data is loaded). Click the adjacent lock button (**🔓** / **🔒**) or adjust the value manually to lock the spot size so it stays constant while adjusting other optimisation parameters; click **🔒** to unlock and automatically recalculate the minimum required spot size.
* **Rep-Rate (Hz)**: Optimised laser firing frequency.
* **Acq Time (ms / s)**: Total mass spectrometer integration cycle time (dwell budget plus hardware overhead).
* **Speed (µm s⁻¹)**: *(Line / Imaging modes)* Optimised stage scan speed, derived from `Speed = Rep-Rate * (Spot Size - Overlap)`. Click the blue **`(int)` / `(1 dp)` / `(2 dp)`** link in the label to cycle speed rounding precision; the displayed Overlap and Dosage automatically adjust to reflect the rounded stage speed.
* **Overlap (µm & %)**: *(Line / Imaging modes)* Spatial overlap between consecutive laser pulses (`(1 dp)` / `(2 dp)` toggle).
* **Dosage (Pulses)**: *(Line / Imaging modes)* Effective number of laser pulses per spot diameter (`(int)` / `(1 dp)` / `(2 dp)` toggle).
* **Est. Time**: Predicted total run duration (`HH:MM:SS`).
* **Est. Data Points** *(Spot / Line modes)* or **Pixels Per Sec** *(Imaging mode)*: Predicted total sweeps or spatial pixel rate.

### 2. Right Panel: Visualisation & Details

#### **Plot Controls**

Arranged in two compact rows above the signal plot:

* **Row 1**:
  * **Theme**: Toggle `Auto`, `Dark`, or `Light` plot styling.
  * **Normalise**: Scales all channel traces to a `0` to `1` range for visual comparison.
  * **Pan / Zoom Y**: Enables vertical zooming and panning.
  * **Auto-Rescale Y**: Defaults to Checked (ephemeral; resets to On when new data is loaded). Automatically scales the Y-axis to fit visible traces.
* **Row 2**:
  * **Auto Select Regions**: Automatically detects the **Background (Blue)** and **Signal (Red)** time windows.
  * **Background (s)** & **Signal (s)**: Spin boxes for fine-tuning the start and end times of both regions.

#### **Interactive Signal Plot (Top)**

* **Region Adjustment**: Hover over any vertical region boundary (background or signal edge; the cursor changes to a horizontal resize arrow `↔`) and **Left-Click & Drag** to move it directly on the plot. If start and end lines are crossed, they are automatically swapped on release.
* **Legend Interactions**:
  * **Single-Click** a channel in the legend to show or hide its trace.
  * **Double-Click** a channel in the legend to isolate it (hiding all other channels); double-click again to restore all channels.
* **Zoom & Pan**: Scroll the **Mouse Wheel** to zoom, **Left-Click & Drag** to pan, and **Double-Click** the plot canvas to reset the view to full extents.

#### **Synchronisation Advisor**

Displays real-time diagnostic badges in the status summary box:

* **🟢 Synchronised**: Perfect integer number of laser pulses per acquisition cycle.
* **🟢 Oversampling Mode**: Time-driven oversampling active.
* **🟢 Stable Steady State**: All active channels meet the target steady-state Lockwood RSD limit.
* **🔴 Unstable Steady State**: Laser repetition rate is below the minimum oversampling threshold (`< 2.0 / Washout`).
* **🟡 Unstable Steady State**: Steady-state Lockwood RSD exceeds the target limit (displays the required acquisition time to reach stability).
* **🟡 Unsynchronised**: Non-integer pulse count per cycle; warns of potential beat-frequency aliasing.

#### **Optimised Dwell Time Distribution (Table)**

Per-channel breakdown at the bottom of the right panel:

* **Mode**:
  * **Auto**: Optimiser distributes dwell time to meet the **Target SNR** for all `Auto` channels.
  * **Even**: Splits the remaining acquisition budget equally across all `Even` channels.
  * **Exclude**: Removes the channel from the method and plot.
  * **Set to Min**: Locks the channel to the instrument's minimum dwell time.
  * **Custom**: Allows typing a custom dwell time directly in the **Optimised Dwell Times** column.
* **Optimised Dwell Times**: Calculated integration time per channel.
* **Initial SNR** & **Resultant SNR**: Comparison of baseline test SNR against predicted SNR under the optimised settings (channels with zero background counts display `(Lc)` to indicate Square Root Transform critical-level estimation).

---

## CSV File Importing & Timestamp Preferences

When importing CSV data files, iolite requires the timestamp column format to match the format configured in iolite's preferences. If an import fails due to a timestamp mismatch, the plugin automatically inspects the iolite log and displays a dialog with the failed timestamp string and suggested format layouts (e.g., `yyyy MM dd hh mm ss` or `dd MM yyyy hh mm ss`).

---

## Scientific References

* **SNR and Detection Limits**:
  * Tanner, S. D. (2010). "Shorter signals for improved signal to noise ratio, the influence of Poisson distribution", *J. Anal. At. Spectrom.*, 25, 405–407.
  * Donard, A., et al. (2015). "Determination of relative rare earth element distributions in very small quantities of uranium ore concentrates using femtosecond UV laser ablation–SF-ICP-MS coupling", *J. Anal. At. Spectrom.*, 30, 2420–2428.
* **Optimisation Workflow & Synchronisation Logic**:
  * Van Malderen, S. J., et al. (2018). "Considerations on data acquisition in laser ablation-inductively coupled plasma-mass spectrometry with low-dispersion interfaces", *Spectrochimica Acta Part B*, 140, 29–34.
  * Van Elteren, J. T., et al. (2019). "Insights into the selection of 2D LA-ICP-MS (multi) elemental mapping conditions", *J. Anal. At. Spectrom.*, 34, 1919.
  * Lockwood, E. T., et al. (2024). "Multiplexed elemental bioimaging with quadrupole ICP-MS and high-frequency laser ablation systems", *J. Anal. At. Spectrom.*, 39, 1125.
* **Single Pulse Response & Washout Analysis**:
  * Ulianov, A., et al. (2015). "The ICPMS signal as a Poisson process: a review of basic concepts", *J. Anal. At. Spectrom.*, 30, 1297–1321.
* **Statistical Framework**:
  * Currie, L. A. (1968). "Limits for qualitative detection and quantitative determination: Application to Radiochemistry", *Anal. Chem.*, 40, 586–593.
  * Currie, L. A. (1972). "The Measurement of Environmental Levels of Rare Gas Nuclides and the Treatment of Very Low-Level Counting Data", *IEEE Trans. Nucl. Sci.*, NS19, (1), 119–126.
