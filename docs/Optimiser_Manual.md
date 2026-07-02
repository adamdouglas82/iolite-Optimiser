# Advanced Optimiser Tab Manual

## Overview

The **Optimiser Tab** determines the best combination of **Spot Size**, **Laser Repetition Rate**, and **Scan Speed** (or Dwell Times) to achieve your analytical goals. It works by "working backwards" from your desired signal quality (SNR) and hardware limits to find the most efficient sampling capability.

---

## User Interface

### 1. Left Panel: Settings & Configuration

#### **Top Controls**

- **Reload Data from iolite**: Refreshes input channels from the main iolite interface.
- **Settings**: Opens the detailed **Hardware Configuration** dialog.

#### **Hardware Configuration Dialog**

Accessible via the "Settings" button. This establishes the physical limits of your system. A version label (e.g. `Version: dev`) is displayed at the bottom-left of this dialog.

- **ICP-MS Hardware**:
  - **Manufacturer/Model**: Presets minimum dwell times and switching speeds.
  - **Custom Models**: Choose "Custom Model" (ICP) or "Custom Laser" (Platform) to manually define limits:
    - **ICP Custom Fields**:
      - **Type**: Quadrupole, Sector-Field, Multi-Collector (MC), or TOF.
      - **Allowed Dwell Times**: (MC only) Comma-separated list of valid integration times.
      - **Min Dwell / Dwell Resolution**: (Quad/TOF) The hardware speed and rounding limits.
      - **Washout Margin (%)**: (TOF only) A percentage "Signal Capture Buffer". Extends the target washout time proportionally to ensure the entire Single Pulse Response is captured natively within the stage's speed budget.
    - **Laser Custom Fields**:
      - **Mode**: Choose "Maximum Rep-Rate" (continuous) or "Discrete" (specific allowed frequencies).
      - **Max Speed**: Physical limit of the sample stage.
      - **Rep-Rate Resolution**: Custom rounding resolution for the laser repetition rate.
- **Laser Hardware**:
  - **Platform/Source**: Presets maximum repetition rates and frequency precision.
  - **Cell Type**: Defines the maximum stage speed (washout behavior is handled separately).

#### **Input Parameters**

Basic setup for your experiment.

- **Initial Spot Size**: The starting point for optimization (often adjusted automatically).
- **Washout**: The time (ms) required for the signal to decay (usually determined via the **SPR Tab**). Supports extremely fast sub-1ms definitions (e.g., `0.25 ms`).
- **Initial Rep-Rate**: A baseline frequency guess (Hz).
- **Number of Analytes**: Specifies how many isotopes/analytes are measured (saved persistently).
- **Suggest PreScan Params**: Triggers a calculation dialog presenting optimal settings for a PreScan scan:
  - **Dwell Time per Analyte**: `washout_ms / num_analytes`, bounded by minimum dwell and rounded to the instrument's dwell resolution.
  - **Target Cycle Time**: `dwell_time * num_analytes`.
  - **Suggested Rep-Rate**: Dynamically swept starting from the oversampling limit ($2.0 / \text{washout}$ s) to select the lowest rate that achieves a statistical Lockwood RSD $\le 5\%$.
  - **Suggested Stage Speed**: `spot_size / cycle_time`.
  - **Suggested Dosage**: `rep_rate * cycle_time` (pulses/pixel).
  - **Suggested Overlap**: Spatial overlap distance and percentage between adjacent pulses.
- **Mode**:
  - **Spot**: Analysis of a single location.
    - **Total Pulses**: Total laser pulses per spot.
  - **Line**: Continuous scanning along a path.
    - **Line Length**: Total distance (µm).
  - **Imaging**: 2D mapping (raster scan).
    - **Image Width / Height**: Dimensions of the area (µm).

#### **Optimisation Parameters**

The "Goals" you want the algorithm to achieve.

- **Sync Strategy**:
  - **Adaptive Integer Sync (Auto)**: Automatically sweeps rep-rates and acquisition times to guarantee a perfect integer pulse count sync (0% rounding error) close to your target.
  - **Strict Pulse Target (Manual)**: Snaps directly to the exact target pulses (may lead to non-integer pulses due to hardware rounding).
  - **Oversampling (Time-Driven)**: Purely time-driven target based on the Lockwood oversampling RSD target, ignoring snapping.
  - **Combined Integer Sync (Auto)**: Combines oversampling with snapping by sweeping frequencies above the oversampling minimum to find a perfect integer count.
- **Sync Target (Pulses)**: The target number of laser shots that must occur during a single mass spectrometer integration (`Acq Time`) to ensure signal stability. Used in integer sync modes.
- **Sync Target (RSD %)**: The target Relative Standard Deviation (RSD) limit for signal stability during steady-state oversampling modes (usually set to `5.0%`).
- **Prefer Exact Pulses per Acq**: (For Auto sync strategies only). When checked, forces the auto-sync search algorithm to prioritize finding a laser rep-rate that delivers *exactly* the requested integer number of pulses per cycle, rather than settling for the closest available integer sync.
- **Target SNR (Sigma)**: The desired Signal-to-Noise ratio for the limiting isotope.
  - **Robust Methodology**: The optimiser calculates a **Detection Ratio** based on the scientific Critical Level ($L_c$), then scales it back to a familiar "Sigma SNR" for the UI. This provides a more rigorous assessment than simple background standard deviation:
    1. **Noise Floor**: Determined by the maximum of the theoretical **Poisson Variance** (counting statistics) and the **Observed Variance** (experimental jitter) of the background.
    2. **Critical Level ($L_c$)**: Calculated using the **Square Root Transform** rule (Stapleton's Rule) for variance stabilization. This defines the 95% confidence threshold for detection.
    3. **Sigma-Equivalent Scaling**: The final Detection Ratio is multiplied by $z\sqrt{2}$ ($\approx 2.326$) so that a value of **10.0** in the UI corresponds to the standard 10-sigma **Limit of Quantification (LoQ)**.
  - This approach is conservative and robust, ensuring that targets represent accurate confidence intervals rather than just multiples of potentially "lucky" low-noise baselines.
- **Min Duty Cycle**: Sets a floor for efficiency (Time Measuring / Total Time).
  - If the calculated duty cycle is too low (mostly settling time), the optimiser will increase dwell times to improve the ratio.
- **Minimum SNR**: Lower limit for acceptance; channels below this may trigger warnings.
- **Scale Signal with Rep-Rate**:
  - **Checked (Linear Scaling)**: Assumes signal intensity (CPS) increases linearly with Frequency (Hz).
    - **Why?** Higher Rep-Rate delivers more pulses per second, ablating more material per second, which generates a higher ion signal.
    - **Use Case**: **Geological Samples** (e.g., blocks/mounts) where there is "infinite" depth. As the laser fires active shots, it drills deeper into fresh material, maintaining the signal increase.
  - **Unchecked (Fixed Sensitivity)**: Assumes signal intensity is determined by Spot Size alone, regardless of frequency.
    - **Use Case**: **Thin Sections / Bio-Imaging**. If the laser already fully ablates the entire thickness of the sample (e.g., a 30µm tissue section) with a single shot, firing faster *does not* generate more signal—it simply consumes the sample faster. Checking this box here would dangerously overestimate your SNR.
- **Avoid Gaps**:
  - **Goal**: Prevents "striping" in images by ensuring the stage doesn't move faster than the laser can fire.
  - **How**: Constraints the calculation so `Speed <= Spot Size * Rep Rate`.
  - **Pros**: Guarantees full surface coverage.
  - **Cons**: May force a slower stage speed or higher Rep-Rate (potentially increasing washout issues).

#### **Optimised Settings (Results Summary)**

The calculated "Ideal" instrument settings.

- **Spot Size**: The recommended beam diameter.
- **Rep-Rate**: The optimal laser frequency (Hz).
- **Speed**: The optimal stage translation speed (µm/s).
- **Acq Time / Budget**: The total integration cycle time calculated to match the laser synchronization.
- **Est. Time**: Total predicted duration of the analysis (HH:MM:SS).
- **Est. Data Points**: (Displayed in Spot/Line modes only) The predicted total number of mass spectrometer measurements that will be collected over the run.
- **Pixels Per Sec**: (Displayed in Imaging mode only) The scanning rate of pixels per second based on the optimized scan speed and spot size.

### 2. Right Panel: Visualisation & Details

#### **Signal Plot (Top)**

- Displays real-time signal traces for the selected channels.
- **Regions**:
  - **Background (Blue)**: Time range for baseline calculation.
  - **Signal (Red)**: Time range for signal intensity calculation.
  - **Auto Select Regions**: Automatically finds the rising/falling edges of the signal.

#### **Plot Controls**

To maximize vertical display efficiency, controls are split into two stacked rows:
- **Row 1**: Theme override dropdown (Auto/Dark/Light), Normalize checkbox, Pan/Zoom Y checkbox, Auto-Rescale Y checkbox, Show Background checkbox.
- **Row 2**: Auto Select Regions checkbox, Background time spinboxes, Signal time spinboxes.

- **Theme**: Toggle Dark/Light mode.
- **Normalize**: Scales all traces to 0-1 range for easy comparison.
- **Auto-Rescale Y**:
  - **Checked (Default/Ephemeral)**: The plot automatically snaps to the full data range. This is reset to ON whenever new data is loaded or the application starts.
  - **Unchecked**: The plot keeps your current zoom level. Switch this off to "lock" your view while adjusting parameters.
- **Pan/Zoom Y**: Enables vertical zooming interaction.
- **Double-Click**: Resets the plot to show the full data range (both X and Y).
- **Interactive Region Adjustment**:
  - **Shift + Click & Drag**: Hold the **Shift** key and click near a region edge (Blue or Red line) to drag it to a new location.
  - **Auto-Swap Validation**: If you drag a region "start" line past its "end" line, the software automatically sorts the values on release. This ensures the start is always less than the end, preventing calculation errors.
  - **Cursor Feedback**: When holding Shift, the cursor will change to a horizontal resize cursor ↔ when you are within the 1% selection "hitbox" of a region edge.

#### **Synchronization Advisor**

Displays real-time diagnostic indicators based on the active Sync Strategy:
- **🟢 Synchronised**: Green indicator representing perfect snapping to integer pulse counts.
- **🟢 Oversampling Mode**: Green indicator indicating no integer snapping is required (pure time-driven).
- **🟢 Stable Steady State**: Green indicator indicating all channels are within the target RSD limits.
- **🔴 Unstable Steady State**: Red indicator informing the user that the rep rate is too low for oversampling ($< 2.0 / \text{washout}$ s).
- **🟡 Unstable Steady State**: Yellow indicator showing that the statistical Lockwood RSD is above the target limit, and listing the required dwell time adjustment to achieve stability.
- **🟡 Unsynchronised**: Yellow indicator warning of potential beating artifacts if a non-integer pulse count is detected during manual overrides.

#### **Optimised Dwell Time Distribution (Table)**

Detailed breakdown per channel.

- **Mode**:
  - **Auto**: Software calculates the optimal dwell.
  - **Even**: Software splits the entire available dwell budget perfectly among all `Even` channels. Calculates exact hardware intervals so no time is lost to drift or rounding errors.
  - **Exclude**: Channel is ignored.
  - **Set to Min**: Forces the channel to the hardware minimum dwell.
  - **Custom**: Allows manual entry of a specific dwell time.
- **Optimised Dwell Times**: The calculated integration time for each isotope.
- **SNR**: Comparison of Initial vs. Resultant Signal-to-Noise Ratios.

---

## CSV File Importing & Preferences

When loading external standard or tuning data via a CSV file, iolite requires the timestamp column format to exactly match the format configuration selected in preferences. If a mismatch occurs, the import will fail to parse and show a dialog:

- **Error popup**: Displays detailed failure information extracted from the iolite application log.
- **Smart Suggestions**: Analyzes the failed timestamp structure and suggests the exact matching preferences layout to configure in your iolite settings (e.g. suggesting `yyyy MM dd hh mm ss` or `dd MM yyyy hh mm ss`).

---

## Workflow Guide

### Step 1: Hardware Setup

1. Click **Settings**.
2. Select your **ICP-MS** and **Laser** models.
3. If your specific model isn't listed, choose "Custom" and enter the **Min Dwell Time** (ICP) and **Max Rep Rate** (Laser) manually.

### Step 2: Load Data

1. Ensure your iolite session has representative data (e.g., a "Tuning" line).
2. Click **Reload Data from iolite**.
3. **Missing Dwell Times?**
    - If your data lacks embedded dwell time metadata (common with some CSV imports), a dialog will appear.
    - **Set to Global**: Apply a single dwell time (e.g., 10ms) to all channels.
    - **Individual**: Manually enter the dwell time for specific channels if they differ.
4. Verify the plot shows the signal traces.
5. **Check Regions**: Ensure the **Background (Blue)** and **Signal (Red)** regions are correctly identifying the baseline and signal plateau.
    - **Manual Selection**: Hold **Shift** and drag the lines in the plot for the fastest adjustment.
    - **Precise Selection**: Use the spinboxes in the left panel.
    - **Auto Select**: Click "Auto Select Regions" for a computer-estimated starting point.
    - Note: Values are automatically swapped if the selection becomes inverted.

### Step 3: Define Goals

1. Set the **Mode** (usually "Imaging" or "Line").
2. Enter your **Washout** time (can be applied from the SPR Tab).
3. Set your target quality metrics (e.g., **Sync Target (Pulses)** = 10, **Target SNR** = 10).
4. Select a **Sync Strategy** to dictate how calculations are prioritized.

### Step 4: Run Optimization

The optimization runs automatically whenever you change a parameter.

1. Watch the **Optimised Settings** panel.
2. The software will display the recommended **Spot Size**, **Rep-Rate**, and **Speed**.
3. Check the **Results Table** at the bottom. Are any important isotopes failing to reach the target SNR?
    - *Orange Highlight*: Cannot reach Target SNR (set to Min).
    - *Blue Highlight*: Hardware limited (cannot go faster).
4. If the Spot Size is too large/small, adjust the **Target SNR** or **Sync Target (Pulses)** to constrain it.

### Step 5: Advanced Logic & Overrides

1. **Min SNR Constraints**:
    - **Quadrupoles**: If a channel cannot meet the Minimum SNR, its dwell time is forced to the *Minimum Dwell* allowed by hardware to save time for other channels.
    - **MC / TOF**: Since all channels are measured simultaneously, the dwell time cannot be reduced individually. Instead, the channel is flagged with an Orange warning in the table.
2. **Min Duty Cycle**:
    - Forces the total integration time to be long enough relative to the magnet settling time (overhead). Useful for ensuring you aren't wasting 50%+ of your analysis time just switching masses.
3. **Manual Overrides**:
    - If you need to force a specific spot size:
    - Change the **Spot Size** spinbox in the "Optimised Settings" panel (it normally says "Auto").
    - This locks the spot size and re-optimizes only the Rep-Rate and Speed.
    - Set it back to 0 (or "Auto") to resume full optimization.
4. If you need to force a specific rep-rate:
    - This is usually derived from the spot size and washout. To influence it, try changing the **Avoid Gaps** setting or **Washout** value.

---

## Tips & Tricks

- **Linking SPR**: Use the SPR Tab first to determine the exact washout time for your cell/tubing, then "Apply" it to the Optimiser.
- **Duty Cycle**: For Sequential ICP-MS (Quadrupoles), a higher Duty Cycle means more efficient use of time. For TOF/MC, this is always 100%.
- **Avoid Gaps**: Essential for imaging. It ensures that `Speed <= Spot Size * Rep Rate`, preventing "striping" artifacts in your map.

---

## Status Messages

The Optimiser provides feedback through the status bar and log to explain why certain settings were chosen or constrained. Understanding these messages can help you troubleshoot why a specific result was calculated.

### General Optimization Messages

- **"Optimum Spot Size [X] µm - Based on Minimum SNR of [Isotope]"**
  - Indicates that the Spot Size was calculated specifically to satisfy the SNR target for the weakest isotope listed.
- **"Optimum Spot Size [X] µm - Overridden to [Y] µm"**
  - Appears when you have manually set the Spot Size spinbox to a specific value (Y), overriding the calculated optimum (X).

### Hardware Constraints

- **"Constraint: Acq Time increased for Rep Rate Limit"**
  - The Acquisition Time had to be extended because the required Dwell Times + Washout would require a Laser Rep-Rate higher than the laser can physically fire.
- **"Constraint: Acq Time increased for Stage Speed Limit"**
  - The Acquisition Time had to be extended because the stage would need to move faster than the cell's maximum speed to cover the spot size in the calculated time.
- **"Constraint: Acq Time increased for Duty Cycle ([X]%)"**
  - The total integration time was increased to ensure that the measurement time vs. overhead (settling) time met your "Min Duty Cycle" target.
- **"Constraint: Acq Time increased for Dwell Budget ([X] ms)"**
  - The total time was increased because the minimum required dwell times for all elected isotopes exceeded the initial budget derived from the Washout Time. Note: For Simultaneous systems (TOF / Multi-Collector), this budget represents the *maximum single* active dwell assignment required. For Sequential systems (Quadrupole), this budget is the *cumulative sum* of all dwell assignments.

### Channel-Specific Warnings

These messages list specific isotopes that are affecting the optimization logic.

- **"The following channels do not meet the minimum SNR target..."**
  - **Orange Warning**: These channels have a calculated SNR below your "Minimum SNR" threshold.
  - For **MC/TOF** systems, this is just a warning.
  - For **Quadrupole** systems, these channels have been forced to the minimum dwell time to save budget.
- **"The following channels cannot be set lower than the hardware minimum ([X] ms)..."**
  - **Blue Warning**: The algorithm wants to reduce the dwell time for these channels (usually high-signal ones) to optimize speed, but is blocked by the hard limit of the instrument (e.g., 1ms or 10ms minimum dwell).
- **"The following channels had zero counts in the background. The detection limit ($L_c$) has been calculated using the Square Root Transform rule..."**
  - **Gray Warning**: Indicates that the background signal was perfectly zero (common in simulated data or extremely short baselines). To prevent division-by-zero errors in the SNR calculation, a theoretical detection limit ($L_c$) was calculated using the Square Root Transform rule for variance stabilization. The results matrix displays these estimated targets with a grey `(Lc)` subscript.

---

## Scientific References

The algorithms and logic used in the iolite Optimiser are based on established statistical principles and laser ablation literature:

- **SNR and Detection Limits**:
  - Tanner, S. D. (2010). "Shorter signals for improved signal to noise ratio, the influence of Poisson distribution" *Journal of Analytical Atomic Spectrometry (JAAS)*, 25, 405–407.
  - Donard, A., et al. "Determination of relative rare earth element distributions in very small quantities of uranium ore concentrates using femtosecond UV laser ablation– SF-ICP-MS coupling" J. Anal. At. Spectrom., (2015), 30, 2420–2428
- **Optimization Workflow & Loop Logic**:
  - Van Malderen, S. J., et al. "Considerations on data acquisition in laser ablation-inductively coupled plasma-mass spectrometry with low-dispersion interfaces" Spectrochimica Acta Part B, (2018), 140, 29–34
  - Van Elteren, J. T., et al. "Insights into the selection of 2D LA-ICP-MS (multi) elemental mapping conditions" J. Anal. At. Spectrom., (2019), 34, 1919
  - Lockwood, E. T. (2024). "Multiplexed elemental bioimaging with quadrupole ICP-MS andhigh-frequency laser ablation systems" J. Anal. At. Spectrom., 39, 1125
- **Single Pulse Response and Washout Analysis**:
  - Ulianov, A., et al. "The ICPMS signal as a Poisson process: a review of basic concepts" J. Anal. At. Spectrom., (2015), 30, 1297–1321
- **Statistical Framework**:
  - Currie, L. A. "Limits for qualitative detection and quantitative determination: Application to Radiochemistry." Anal. Chem., (1968), 40, 586-593.
  - Currie, L. A. "The Measurement of Environmental Levels of Rare Gas Nuclides and the Treatment of Very Low-Level Counting Data", IEEE Trans. Nucl. Sci., (1972), NS19, (1), 119-126.
