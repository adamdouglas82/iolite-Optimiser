# Optimisation Workflow Manual

## Overview

This guide provides a step-by-step workflow for using the **iolite Optimiser**. It combines **Single Pulse Response (SPR)** washout characterisation, **PreScan** parameter calculation, **Method Optimisation**, and **Pulse Train Simulation** to tune your laser ablation ICP-MS method for your hardware and analytical goals.

---

## Phase 1: Determine Washout (Determine SPR Tab)

Before optimising spot size or scan speed, characterise your system's Single Pulse Response (SPR) to measure its washout time.

### Step 1.1: Acquire and Load SPR Test Data

1. **Acquire Test Data**: Run a single-pulse test ablation on your instrument using a low laser repetition rate of **1 Hz** and a **10 ms** integration time measuring a single high-abundance mass (or a short list of masses). This records cleanly separated, high-resolution profiles of individual laser pulses.
2. Open the **iolite Optimiser** plugin from the left sidebar or *Tools* menu in iolite.
3. Select the **Determine SPR** tab.
4. Click **Import SPR File** and select your SPR data file.

### Step 1.2: Select & Optimise Peak Detection

1. **Select Channel**: Choose a high-abundance isotope channel (e.g., `U238` or `Th232`), or check the **Maximum Observed SPR** panel to see which channel exhibits the longest washout.
2. **Verify Detection**: Inspect the **SPR Plot**. Valid pulses are marked with an orange `X`, and the active threshold is shown by the horizontal dashed **Peak Cutoff** line.
   * *Too much noise detected as peaks?* Hover over and drag the dashed **Peak Cutoff** line up on the plot (or uncheck **Auto** and increase the value), or increase **Min. Peak Spacing** (default `900 ms`).
   * *Missed peaks?* Drag the dashed **Peak Cutoff** line down.
3. **Verify Background Subtraction**: If your file has at least `5 seconds` of gas blank at the start or end of the recording, **Subtract Background** automatically selects that baseline region. If no `>= 5 s` gas blank is present, **Subtract Background** is automatically unchecked.
4. **Clean Data**: Click on any outlier peak marker or number label on the plot (or enable **Auto-Exclude Wide Outliers**) to exclude it from the statistics.
5. **Apply Washout**:
   * Review the **FW 0.1M (10%)**, **FW 0.01M (1%)**, or **Maximum Observed SPR** statistics in the left panel.
   * Click **Apply Average**, **Apply Composite**, or **Apply Max** to send that washout time directly to the **Method Optimiser** tab.
   * *Tip*: **Apply Composite** uses the stacked average pulse shape and is typically the most representative value for system washout.

---

## Phase 2: Determine PreScan Parameters (Setting Initial Conditions)

Once the system washout time is known, you can calculate ideal conditions for a quick PreScan over your sample or reference material.

### Step 2.1: Calculate Suggestions

1. Switch to the **Method Optimiser** tab and click **Settings** to select your **ICP-MS** and **Laser** hardware models.
2. In the **Input Parameters** group:
   * Enter your desired **Initial Spot Size (µm)** (e.g., `15 µm`) and verify the **Washout (ms)** value applied from Phase 1.
   * Enter the **Number of Analytes** in your proposed mass spectrometer method.
3. Click **Suggest PreScan Params**.
4. The dialog displays:
   * **Dwell Time per Analyte**: Bounded by hardware minimums and rounded to your instrument's dwell resolution.
   * **Target Cycle Time**: Total integration sweep time.
   * **Suggested Rep-Rate**: Lowest valid laser frequency achieving `<= 5%` steady-state Lockwood RSD under oversampling.
   * **Suggested Stage Speed**, **Dosage**, and **Overlap**: Recommended scan parameters (capped at the cell's maximum stage speed).
5. Run your PreScan on the instrument using these parameters to record a representative baseline and ablation signal.

---

## Phase 3: Load Sensitivity Data (Method Optimiser Tab)

### Step 3.1: Import PreScan / Tuning File

1. On the **Method Optimiser** tab, click **Import File** and select your PreScan or tuning data file.
2. **Dwell Time Check**:
   * If the file format does not contain embedded dwell time metadata, the **Configure Channel Dwell Times** dialog appears.
   * Enter a **Global Dwell Time** (or individual channel dwell times), select the channels to include (or load a saved channel preset), and click **OK**.
3. **Verify Background & Signal Regions**:
   * Inspect the signal plot to confirm the **Background (Blue)** region covers the gas blank and the **Signal (Red)** region covers the steady ablation plateau.
   * **Quick Adjustment**: Hold **Shift** and **Left-Click & Drag** the vertical region edges directly on the plot, or edit the **Background (s)** and **Signal (s)** spin boxes.

---

## Phase 4: Run Optimisation & Verify Stability

### Step 4.1: Define Analytical Goals

1. **Input Parameters**: Confirm **Initial Spot Size (µm)** and **Initial Rep-Rate (Hz)** match the conditions used to acquire the loaded file, and select your target **Mode** (`Spot`, `Line`, or `Imaging`).
2. **Optimisation Parameters**:
   * **Sync Strategy**:
     * **Adaptive Integer Sync (Auto)**: Sweeps repetition rates and acquisition times to guarantee an exact integer number of pulses per cycle.
     * **Strict Pulse Target (Manual)**: Locks directly to the requested pulse count.
     * **Oversampling (Time-Driven)**: Optimises for steady-state oversampling (`Sync Target (RSD %)`).
     * **Combined Integer Sync (Auto)**: Sweeps above the oversampling minimum to achieve both steady-state stability and integer pulse synchronisation.
   * **Sync Target (Pulses)** or **Sync Target (RSD %)**: Set your pulse count per cycle (e.g., `5` pulses) or target steady-state RSD (e.g., `5.0%`).
   * **Target SNR (Sigma)** & **Minimum SNR (Sigma)**: Set your required detection quality (e.g., `Target SNR = 10`).
   * **Scale Signal with Rep-Rate**: Check for bulk mounts where higher repetition rates proportionally increase signal intensity; leave unchecked for thin sections that exhaust cleanly in a few shots.

### Step 4.2: Interpret Results

The optimiser recalculates automatically whenever any parameter changes:

1. **Optimised Settings Panel**:
   * Review the recommended **Spot Size**, **Rep-Rate**, **Acq Time**, **Speed**, **Overlap**, and **Dosage**.
   * Click the blue **`(int)` / `(1 dp)` / `(2 dp)`** links next to **Speed**, **Overlap**, or **Dosage** to match the decimal precision accepted by your laser control software.
2. **Synchronisation Advisor**:
   * Confirm the advisor badges show green indicators (**🟢 Synchronised**, **🟢 Oversampling Mode**, or **🟢 Stable Steady State**).
3. **Optimised Dwell Time Distribution Table**:
   * Review the individual **Optimised Dwell Times** and **Resultant SNR** for each isotope.
   * **Orange Rows**: Channels that cannot reach **Minimum SNR** (forced to minimum dwell on sequential systems).
   * **Blue Rows**: High-intensity channels constrained by the hardware minimum dwell time.
   * Change any channel's **Mode** (`Auto`, `Even`, `Exclude`, `Set to Min`, `Custom`) to tailor how the cycle budget is allocated.
4. **Locking / Overriding Spot Size**:
   * To maintain the currently selected spot size while adjusting other parameters, click the lock button (**🔓** -> **🔒**) next to **Spot Size (µm)** in the **Optimised Settings** panel (or adjust the spin box to a manual diameter, which engages the lock automatically). Click **🔒** to unlock and resume automatic minimum spot size calculation.

### Step 4.3: Verify in the Pulse Train Simulator

1. Switch to the **Pulse Train Simulator** tab.
2. Inspect the simulated detector traces and the **Pulse Simulation Stability** table at the bottom to verify the predicted steady-state **RSD (%)** and confirm there are no beat-frequency artefacts before running your final experiment.
