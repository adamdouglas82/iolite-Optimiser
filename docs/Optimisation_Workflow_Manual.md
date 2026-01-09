# Optimisation Workflow Manual

## Overview

This guide provides a comprehensive, step-by-step workflow for using the **iolite Optimiser**. It combines the **Single Pulse Response (SPR)** analysis with the **Optimiser** to ensure your laser ablation parameters are perfectly tuned for your specific hardware and analytical goals.

---

## Phase 1: Determine Washout (SPR Tab)

Before optimizing spot size or speed, you must know your system's washout time to prevent signal mixing.

### Step 1.1: Load SPR Data

1. Ensure your iolite session contains a file with single laser pulses (e.g., a "Test" or "Tune" line created with a low repetition rate, typically 1 Hz).
2. Open the **iolite Optimiser** plugin from the *Tools* menu.
3. Go to the **SPR Tab**.
4. Click **Reload Data from iolite** to ingest the current selection.

### Step 1.2: Select & Optimize

1. **Select Channel**: Choose a high-abundance isotope (e.g., `U238`, `Th232`) from the dropdown.
2. **Verify Detection**: Look at the **SPR Plot**. Valid pulses should be marked with an orange `X`.
    - *Too much noise?* Increase **Peak Cutoff** (uncheck Auto) or **Min. Distance**.
    - *Missed peaks?* Decrease **Peak Cutoff**.
3. **Clean Data**: If you see double-peaks or artifacts, click on them in the plot to **Exclude** them from the calculation.

4. **Apply Washout**:
    - Review the **FW 0.1M (10%)** statistics in the left panel.
    - Click **Apply Average, Max, or Composite to Optimiser**.
    - *Tip*: The **Composite** value is often the most stable representative of your system's washout.

---

## Phase 2: Optimiser Setup (Optimiser Tab)

### Step 2.1: Hardware Configuration

1. Switch to the **Optimiser Tab**.
2. Click the **Settings** button.
3. **ICP-MS**: Select your Manufacturer and Model.
    - *Custom?* Choose "Custom Model" and enter your **Min Dwell Time** (e.g., for TOF or generic Quad).
4. **Laser**: Select your Platform and Cell Type.
    - *Custom?* Choose "Custom Laser" and manually limit the **Max Rep Rate** and **Max Speed**.

### Step 2.2: Load Sensitivity Data

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

## Phase 3: Run Optimization

### Step 3.1: Define Analytical Goals

1. **Mode**: Select your analysis type (usually **Imaging** or **Line**).
2. **Washout**: Verify this matches the value you applied from the SPR tab (Phase 1).
3. **Parameters**:
    - **Pulses per Dwell**: Set the stability requirement (e.g., 10 shots/dwell).
    - **Target SNR**: Set the desired quality (e.g., 10 Sigma).
    - **Avoid Gaps**: Check this for Imaging to ensure full surface coverage.

- **Navigation Tip**: **Double-click** any plot at any time to instantly reset the view to its full extents. Use the **Mouse Wheel** to zoom and **Left-Click & Drag** to pan.

### Step 3.2: Interpret Results

The software calculates the optimal settings automatically.

1. **Review Settings**: Look at the **Optimised Settings** panel for your ideal **Spot Size**, **Rep-Rate**, and **Speed**.
2. **Check Constraints**: Look at the **Results Table**.
    - **Orange Rows**: These isotopes cannot meet the Target SNR even at the maximum allowable settings.
    - **Blue Rows**: These are limited by your hardware speed (Min Dwell).
3. **Refine**:
    - *Spot Size too big?* Reduce the **Target SNR** or **Pulses per Dwell**.
    - *Scan too slow?* Uncheck **Avoid Gaps** (if acceptable) or increase **Spot Size** manually.

### Step 3.3: Manual Overrides & Advanced Logic

- **Force Spot Size**: If you strictly need a specific beam diameter (e.g., 25µm), enter it in the **Spot Size** spinbox (toggle off "Auto"). Valid Rep-Rate and Speed will be re-calculated for that fixed spot.
- **Min SNR (Hardware Differences)**:
  - *Quadrupoles*: Channels failing SNR are set to the minimum dwell to save time.
  - *MC/TOF*: Channels failing SNR are flagged (Orange) but dwell times remain fixed as they are simultaneous detectors.
