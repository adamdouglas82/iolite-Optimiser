# Optimisation Workflow Manual

## Overview

This guide provides a comprehensive, step-by-step workflow for using the **iolite Optimiser**. It combines the **Single Pulse Response (SPR)** analysis with the **Optimiser** to ensure your laser ablation parameters are perfectly tuned for your specific hardware and analytical goals.

---

## Phase 1: Determine Washout (SPR Tab)

Before optimizing spot size or scan speed, you must characterize your system's Single Pulse Response (SPR) to determine the washout time.

### Step 1.1: Acquire and Load SPR Test Data

1. **Acquire Test Data**: Run a single-pulse test ablaton on your instrument using a low laser repetition rate of **1 Hz** and a **10 ms** integration time measuring a single mass (isotope). This ensures you record cleanly separated, high-resolution profiles of individual laser pulses.
2. Open the **iolite Optimiser** plugin from the *Tools* menu in iolite.
3. Go to the **SPR Tab**.
4. Click **Reload Data from iolite** to ingest the current selection.

### Step 1.2: Select & Optimize Peak Detection

1. **Select Channel**: Choose the high-abundance isotope channel (e.g. `U238`, `Th232`) used during the single-pulse run.
2. **Verify Detection**: Look at the **SPR Plot**. Valid pulses should be marked with an orange `X`.
    - *Too much noise?* Increase **Peak Cutoff** (uncheck Auto) or **Min. Distance**.
    - *Missed peaks?* Decrease **Peak Cutoff**.
3. **Clean Data**: If you see double-peaks or artifacts, click on them in the plot to **Exclude** them from the calculation.
4. **Apply Washout**:
    - Review the **FW 0.1M (10%)** statistics in the left panel.
    - Click **Apply Average, Max, or Composite to Optimiser**.
    - *Tip*: The **Composite** value is often the most stable representative of your system's washout.

---

## Phase 2: Determine PreScan Parameters (Setting Initial Conditions)

Once the Single Pulse Response washout is determined, you can use the advisor tool to set up the initial conditions for your instrument PreScan.

### Step 2.1: Calculate Suggestions

1. Switch to the **Optimiser Tab** of the plugin.
2. In the **Input Parameters** group:
   - Enter your target **Initial Spot Size** (e.g. 15µm) and the **Washout** decay time applied from Phase 1.
   - Enter the **Number of Analytes** representing the total number of isotopes in your proposed mass spectrometer method.
3. Click the **Suggest PreScan Params** button.
4. The suggestions dialog will display:
   - **Dwell Time per Isotope**: Bounded by hardware minimums and rounded to the instrument's resolution.
   - **Target Cycle Time**: Overall integration loop time.
   - **Suggested Rep-Rate**: Firing frequency that achieves low statistical fluctuation ($\le 5\%$ RSD) under steady-state oversampling.
   - **Suggested Stage Speed** and **Overlap**: Recommended values to prevent striping.
5. Apply these settings to your instrument's control software to execute the PreScan and generate the primary calibration/tuning data file.

---

## Phase 3: Optimiser Setup (Optimiser Tab)

### Step 3.1: Hardware Configuration

1. Switch to the **Optimiser Tab**.
2. Click the **Settings** button.
3. **ICP-MS**: Select your Manufacturer and Model.
    - *Custom?* Choose "Custom Model" and enter your **Min Dwell Time** and **Dwell Resolution** (e.g., for TOF or generic Quad).
4. **Laser**: Select your Platform and Cell Type.
    - *Custom?* Choose "Custom Laser" and manually limit the **Max Rep Rate**, **Max Speed**, and **Rep-Rate Resolution**.

### Step 3.2: Load Sensitivity Data

1. Ensure your iolite session has a representative ablation signal (e.g., a "Tuning" line or standard block).
2. Click **Reload Data from iolite**.
3. **Dwell Time Check**:
    - If your data lacks dwell time metadata, a dialog will appear.
    - Choose **Set to Global** (e.g., 10ms) or enter **Individual** times per channel to match your instrument method.
4. **Region Selection**:
    - Check the Signal Plot.
    - Ensure the **Background (Blue)** region covers the gas blank.
    - Ensure the **Signal (Red)** region covers the steady ablation signal.
    - **Quick Adjustment**: Hold **Shift** and drag the region lines directly on the plot. The software will automatically swap start/end points if you cross them.
    - Use **Auto Select Regions** if they are incorrect.

---

## Phase 4: Run Optimization

### Step 4.1: Define Analytical Goals

1. **Mode**: Select your analysis type (usually **Imaging** or **Line**).
2. **Washout**: Verify this matches the value you applied from the SPR tab (Phase 1).
3. **Parameters**:
    - **Sync Strategy**: Choose how rep-rates and integration times are optimized:
      - **Adaptive Integer Sync (Auto)**: Sweeps frequencies to guarantee integer pulses per cycle.
      - **Strict Pulse Target (Manual)**: Hard-locks to the target pulses.
      - **Oversampling (Time-Driven)**: Prioritizes steady-state oversampling without snapping.
      - **Combined Integer Sync (Auto)**: Sweeps above the oversampling minimum to find a perfect integer sync count.
    - **Sync Target (Pulses)**: Set the strict integration synchronisation requirement (e.g., 10 shots/dwell).
    - **Sync Target (RSD %)**: Set the target stability margin (e.g. 5%) under oversampling modes.
    - **Prefer Exact Pulses per Acq**: Check this box to force the algorithm to search only for the exact integer target.
    - **Dosage**: Set the spatial density requirement (e.g., 10 shots/pixel).
        - Keep **Sync Dosage** checked for square pixels where Dosage = Pulses/Dwell.
        - Uncheck **Sync Dosage** to manually decouple the scan speed from the mass spectrometer's integration rate (e.g., for rectangular pixels or TOF instruments).
    - **Target SNR**: Set the desired quality (e.g., 10 Sigma).
    - **Avoid Gaps**: Check this for Imaging to ensure full surface coverage.

- **Navigation Tip**: **Double-click** any plot at any time to instantly reset the view to its full extents. Use the **Mouse Wheel** to zoom and **Left-Click & Drag** to pan.

### Step 4.2: Interpret Results

The software calculates the optimal settings automatically.

1. **Review Settings**: Look at the **Optimised Settings** panel for your ideal **Spot Size**, **Rep-Rate**, and **Speed**.
2. **Review Synchronization Advisor**:
   - Check the advisor diagnostics in the results summary. Ensure it displays green icons indicating **Synchronised**, **Oversampling Mode**, or **Stable Steady State** (all channels meeting target RSD). If a yellow/red status warning appears, adjust dwell or rep-rates accordingly.
3. **Check Constraints**: Look at the **Results Table**.
    - **Orange Rows**: These isotopes cannot meet the Target SNR even at the maximum allowable settings.
    - **Blue Rows**: These are limited by your hardware speed (Min Dwell).
4. **Refine**:
    - *Spot Size too big?* Reduce the **Target SNR** or **Sync Target (Pulses)**.
    - *Scan too slow?* Uncheck **Avoid Gaps** (if acceptable) or increase **Spot Size** manually.

### Step 4.3: Manual Overrides

- **Force Spot Size**: If you strictly need a specific beam diameter (e.g., 25µm), enter it in the **Spot Size** spinbox (toggle off "Auto"). Valid Rep-Rate and Speed will be re-calculated for that fixed spot.
