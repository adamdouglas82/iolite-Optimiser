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

Accessible via the "Settings" button. This establishes the physical limits of your system.

- **ICP-MS Hardware**:
  - **Manufacturer/Model**: Presets minimum dwell times and switching speeds.
  - **Custom Models**: Choose "Custom Model" (ICP) or "Custom Laser" (Platform) to manually define limits:
    - **ICP Custom Fields**:
      - **Type**: Quadrupole, Sector-Field, Multi-Collector (MC), or TOF.
      - **Allowed Dwell Times**: (MC only) Comma-separated list of valid integration times.
      - **Min Dwell / Precision**: (Quad/TOF) The hardware speed limits.
    - **Laser Custom Fields**:
      - **Mode**: Choose "Maximum Rep-Rate" (continuous) or "Discrete" (specific allowed frequencies).
      - **Max Speed**: Physical limit of the sample stage.
- **Laser Hardware**:
  - **Platform/Source**: Presets maximum repetition rates and frequency precision.
  - **Cell Type**: Defines the maximum stage speed (washout behavior is handled separately).

#### **Input Parameters**

Basic setup for your experiment.

- **Initial Spot Size**: The starting point for optimization (often adjusted automatically).
- **Washout**: The time (ms) required for the signal to decay (usually determined via the **SPR Tab**).
- **Initial Rep-Rate**: A baseline frequency guess (Hz).
- **Mode**:
  - **Spot**: Analysis of a single location.
    - **Number of Shots**: Total laser pulses per spot.
  - **Line**: Continuous scanning along a path.
    - **Line Length**: Total distance (µm).
  - **Imaging**: 2D mapping (raster scan).
    - **Image Width / Height**: Dimensions of the area (µm).

#### **Optimisation Parameters**

The "Goals" you want the algorithm to achieve.

- **Pulses per Dwell Time**: The number of laser shots that must occur during a single mass spectrometer integration (dwell) to ensure signal stability.
- **Target SNR (Sigma)**: The desired Signal-to-Noise ratio for the limiting isotope.
  - **Robust Methodology**: The optimiser calculates a **Detection Ratio** based on the scientific Critical Level ($L_c$), then scales it back to a familiar "Sigma SNR" for the UI. This provides a more rigorous assessment than simple background standard deviation:
    1. **Noise Floor**: Determined by the maximum of the theoretical **Poisson Variance** (counting statistics) and the **Observed Variance** (experimental jitter) of the background.
    2. **Critical Level ($L_c$)**: Calculated using the **Square Root Transform** rule (Stapleton's Rule) for variance stabilization. This defines the 95% confidence threshold for detection.
    3. **Sigma-Equivalent Scaling**: The final Detection Ratio is multiplied by $z\sqrt{2}$ ($\approx 2.326$) so that a value of **10.0** in the UI corresponds to the standard 10-sigma **Limit of Quantification (LoQ)**.
  - This approach is conservative and robust, ensuring that targets represent accurate confidence intervals rather than just multiples of potentially "lucky" low-noise baselines.
- **Min Duty Cycle**: Sets a floor for efficiency (Time Measuring / Total Time).
  - If the calculated duty cycle is too low (mostly settling time), the optimiser will increase dwell times to improve the ratio.
- **Minimum SNR**: Lower limit for acceptance; channels below this may trigger warnings.
- **Image Dimensions**: (See **Mode** above).
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

### 2. Right Panel: Visualisation & Details

#### **Signal Plot (Top)**

- Displays real-time signal traces for the selected channels.
- **Regions**:
  - **Background (Blue)**: Time range for baseline calculation.
  - **Signal (Red)**: Time range for signal intensity calculation.
  - **Auto Select Regions**: Automatically finds the rising/falling edges of the signal.

#### **Plot Controls**

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

#### **Optimised Dwell Time Distribution (Table)**

Detailed breakdown per channel.

- **Mode**:
  - **Auto**: Software calculates the optimal dwell.
  - **Exclude**: Channel is ignored.
  - **Set to Min**: Forces the channel to the hardware minimum dwell.
  - **Custom**: Allows manual entry of a specific dwell time.
- **Optimised Dwell Times**: The calculated integration time for each isotope.
- **SNR**: Comparison of Initial vs. Resultant Signal-to-Noise Ratios.

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
3. Set your target quality metrics (e.g., **Pulses per Dwell Time** = 10, **Target SNR** = 10).

### Step 4: Run Optimization

The optimization runs automatically whenever you change a parameter.

1. Watch the **Optimised Settings** panel.
2. The software will display the recommended **Spot Size**, **Rep-Rate**, and **Speed**.
3. Check the **Results Table** at the bottom. Are any important isotopes failing to reach the target SNR?
    - *Orange Highlight*: Cannot reach Target SNR (set to Min).
    - *Blue Highlight*: Hardware limited (cannot go faster).
4. If the Spot Size is too large/small, adjust the **Target SNR** or **Pulses per Dwell** to constrain it.

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
- **"Calculated Pulses per Dwell Time: [X] (Error: [Y]%)"**
  - Shows the *actual* number of laser pulses that will occur during one mass spectrometer integration. The error percentage shows the deviation from your target "Pulses per Dwell" setting. Small errors (<10%) are usually acceptable.

### Hardware Constraints

- **"Constraint: Acq Time increased for Rep Rate Limit"**
  - The Acquisition Time had to be extended because the required Dwell Times + Washout would require a Laser Rep-Rate higher than the laser can physically fire.

- **"Constraint: Acq Time increased for Stage Speed Limit"**
  - The Acquisition Time had to be extended because the stage would need to move faster than the cell's maximum speed to cover the spot size in the calculated time.
- **"Constraint: Acq Time increased for Duty Cycle ([X]%)"**
  - The total integration time was increased to ensure that the measurement time vs. overhead (settling) time met your "Min Duty Cycle" target.
- **"Constraint: Acq Time increased for Dwell Budget ([X] ms)"**
  - The total time was increased because the sum of the minimum required dwell times for all elected isotopes exceeded the initial budget derived from the Washout Time.

### Channel-Specific Warnings

These messages list specific isotopes that are affecting the optimization logic.

- **"The following channels do not meet the minimum SNR target..."**
  - **Orange Warning**: These channels have a calculated SNR below your "Minimum SNR" threshold.
  - For **MC/TOF** systems, this is just a warning.
  - For **Quadrupole** systems, these channels have been forced to the minimum dwell time to save budget.
- **"The following channels cannot be set lower than the hardware minimum ([X] ms)..."**
  - **Blue Warning**: The algorithm wants to reduce the dwell time for these channels (usually high-signal ones) to optimize speed, but is blocked by the hard limit of the instrument (e.g., 1ms or 10ms minimum dwell).
- **"The following channels had zero counts in the background. The detection limit (Lc) has been calculated using the Square Root Transform rule..."**
  - **Gray Warning**: Indicates that the background signal was perfectly zero (common in simulated data or extremely short baselines). To prevent division-by-zero errors in the SNR calculation, a theoretical detection limit ($L_c$) was calculated using the Square Root Transform rule for variance stabilization.

---

## Scientific References

The algorithms and logic used in the iolite Optimiser are based on established statistical principles and laser ablation literature:

- **SNR and Detection Limits**:
  - Tanner, S. D. (2010). "Shorter signals for improved signal to noise ratio, the influence of Poisson distribution" *Journal of Analytical Atomic Spectrometry (JAAS)*, 25, 405–407.
  - Donard, A., et al. "Determination of relative rare earth element distributions in very small quantities of uranium ore concentrates using femtosecond UV laser ablation– SF-ICP-MS coupling" J. Anal. At. Spectrom., (2015), 30, 2420–2428
- **Optimization Workflow & Loop Logic**:
  - Van Malderen, S. J., et al. "Considerations on data acquisition in laser ablation-inductively coupled plasma-mass spectrometry with low-dispersion interfaces" Spectrochimica Acta Part B, (2018), 140, 29–34
  - Van Elteren, J. T., et al. "Insights into the selection of 2D LA-ICP-MS (multi) elemental mapping conditions" J. Anal. At. Spectrom., (2019), 34, 1919
- **Single Pulse Response and Washout Analysis**:
  - Ulianov, A., et al. "The ICPMS signal as a Poisson process: a review of basic concepts" J. Anal. At. Spectrom., (2015), 30, 1297–1321
- **Statistical Framework**:
  - Currie, L. A. "Limits for qualitative detection and quantitative determination: Application to Radiochemistry." Anal. Chem., (1968), 40, 586-593.
  - Currie, L. A. "The Measurement of Environmental Levels of Rare Gas Nuclides and the Treatment of Very Low-Level Counting Data", IEEE Trans. Nucl. Sci., (1972), NS19, (1), 119-126.
  - Based on the Square Root Transform rule for variance stabilization in counting statistics.


