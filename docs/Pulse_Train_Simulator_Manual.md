# Pulse Train Simulator Tab Manual

## Overview

The **Pulse Train Simulator** tab simulates the detector response signals resulting from a series of laser pulses. By modelling individual pulse profiles, it visualizes the composite signal at the detector—helping users understand how parameters like repetition rate, dwell/acquisition times, and washout decay interact. It is particularly useful for identifying signal aliasing, carryover, and transient peak resolution.

---

## User Interface

The tab is divided into two primary sections: **Simulation Settings & Overrides** (Left) and the **Interactive Simulator Plot** (Right).

### 1. Left Panel: Settings & Adjustments

#### **Simulation Settings**

Establish the temporal boundaries and pulse shape characteristics of the simulation.

- **Background Time**: The simulated duration (seconds) preceding active laser ablation, representing the baseline gas blank.
- **Signal Duration**: The active duration (seconds) of simulated laser firing/ablation.
- **Pulse Shape**:
  - **Real Composite Peak**: Reconstructs the pulse train utilizing the actual average pulse shape extracted from the **SPR tab** data. (Only available if SPR data has been loaded and analyzed).
  - **Model Washout Peak (Lognormal)**: Models the signal using a theoretical lognormal distribution curve. (Default shape if no composite peak is loaded).
  - **Sawtooth Approximation**: A geometric sawtooth approximation representing linear rise and decay signals.

#### **Overrides / Manual Adjustment**

Manually adjust the simulation parameters. By default, these parameters are initialized to match the values calculated on the **Optimiser tab**.

- **Reset to Optimum**: Instantly reverts all overrides and aligns the simulator inputs back to the optimal settings computed by the Optimiser tab.
- **Rep Rate**: The laser repetition frequency (Hz).
- **Acq Time**: The mass spectrometer acquisition integration cycle time (ms).
- **Pulses / Acq**: The average number of laser shots per acquisition cycle.
- **Washout**: The decay time (ms) for the pulse response.
- **Interactive Recalculation**: Adjusting one of **Rep Rate**, **Acq Time**, or **Pulses / Acq** will dynamically update the other parameters to satisfy the relationship:
  $$\text{Pulses} = \text{Rep Rate} \times \text{Acq Time}$$

---

### 2. Right Panel: Visualisation

#### **Plot Controls**

- **Theme**: Select `Auto` (inherits system theme), `Dark`, or `Light` styling.
- **Normalize**: Scales the trace amplitude to a 0–1 range, making it easy to examine shapes and relative carryover.
- **Pan / Zoom Y**: Enables mouse interactions to zoom or pan along the vertical Y-axis.
- **Auto-Rescale Y**: When enabled, automatically snaps the vertical axis to display the full simulated signal range whenever parameters change.
- **Show Background**: Toggles visualization of the baseline background noise level.

#### **Interactive Simulation Plot**

- Renders the resulting simulated detector signal trace over time.
- **Legend**:
  - Displays the active channels and simulated traces.
  - **Interaction**: Click on items in the legend to toggle their visibility on the canvas.
- **Interactive Dwell Boundaries**:
  - Displays vertical dotted lines representing individual integration slice boundaries.
  - Helps visualize how "jitter" and alignment changes occur from cycle to cycle when unsynchronized.

---

## Advanced Options

The simulator implements additional advanced options to model complex physical behaviors:

- **Washout Scaling**: Scales the decay rate of the pulse profiles to simulate changes in gas flow or tubing lengths.
- **Simulator dt (Resolution) Scaling**: Adjusts the resolution of the integration grid to resolve fine details of extremely fast transient signals and prevent aliasing artifacts in modeling.
- **Hide Pulse Train Channels**: Toggle to show/hide individual pulse channels to keep the layout clean when analyzing many isotopes.
- **SNR/CPS-based Intensity Scaling**: Adjusts simulated amplitudes to reflect actual counts-per-second (CPS) and signal-to-noise ratios calculated from your loaded iolite standard block.
- **Sequential Sampling**: Models sequential mass spectrometer sweeps (rotating through quadrupole mass selections) to visualize how skew/drift manifests in multi-element sequential analysis.
