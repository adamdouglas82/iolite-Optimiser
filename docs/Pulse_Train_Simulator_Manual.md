# Pulse Train Simulator Tab Manual

## Overview

The **Pulse Train Simulator** tab simulates the detector response resulting from a series of laser pulses. By modelling individual pulse profiles and sequential mass spectrometer sweeps, it visualises the composite signal at the detector — helping you evaluate how repetition rate, dwell/acquisition times, and washout decay interact to cause or prevent signal aliasing, beat patterns, and carryover.

---

## User Interface

The tab is divided into two primary sections: **Simulation Settings & Overrides** (Left) and **Visualisation & Stability Statistics** (Right).

### 1. Left Panel: Settings & Adjustments

#### **Simulation Settings**

Establish the temporal boundaries and pulse shape characteristics of the simulation:

* **BG Time**: The simulated baseline duration (default `5.0 s`) before and after active laser ablation, representing the gas blank.
* **Signal Duration**: The active duration (default `30.0 s`, minimum `15.0 s`) of simulated laser firing.
* **Pulse Shape**:
  * **Real Composite Peak**: Reconstructs the pulse train using the actual normalised composite pulse shape extracted from the **Determine SPR** tab (automatically selected when SPR data is available).
  * **Model Washout Peak (Lognormal)**: Models the signal using a theoretical lognormal washout curve (default when no SPR composite peak is loaded).
  * **Sawtooth Approximation**: A geometric triangular/sawtooth approximation representing linear rise and decay profiles.

#### **Overrides / Manual Adjustment**

By default, simulation parameters are synchronised with the optimal values computed on the **Method Optimiser** tab. You can manually adjust them here to explore "what-if" scenarios:

* **Reset to Optimum**: Reverts all manual overrides back to the current settings from the **Method Optimiser** tab.
* **Rep Rate**: Laser repetition frequency (`Hz`).
* **Acq Time**: Mass spectrometer acquisition cycle time (`ms`).
* **Pulses / Acq**: Average number of laser shots per acquisition cycle (`Pulses = Rep Rate * Acq Time`).
* **Washout**: Target decay time (`ms`) used to scale the individual pulse profile width.
* **Interactive Recalculation**: Adjusting **Rep Rate**, **Acq Time**, or **Pulses / Acq** dynamically updates the linked parameters according to `Pulses = Rep Rate * (Acq Time / 1000)`.

---

### 2. Right Panel: Visualisation & Stability Statistics

#### **Plot Controls**

* **Theme**: Select `Auto` (inherits system theme), `Dark`, or `Light` styling.
* **Normalise**: Scales each trace to a `0` to `1` range so pulse shapes and relative oscillation amplitudes can be compared across isotopes with different abundances.
* **Pan / Zoom Y**: Enables vertical zooming and panning with the mouse.
* **Auto-Rescale Y**: Automatically snaps the vertical axis to fit the visible simulated traces whenever parameters change.
* **Show Background**: Toggles display of the baseline background noise floor.

#### **Interactive Simulation Plot (Top)**

* Renders the simulated detector signal traces over time for all active channels.
* **Legend Interactions**:
  * **Single-Click** a legend item to show or hide that channel's trace.
  * **Double-Click** a legend item to isolate that single channel; double-click again to restore all channels.
* **Zoom & Pan**: Use the **Mouse Wheel** to zoom, **Left-Click & Drag** to pan, and **Double-Click** the plot canvas to reset the view to full extents.

#### **Pulse Simulation Stability (Signal Area ± 1s crop) Table (Bottom)**

Evaluates the steady-state stability of the simulated signal across the central ablation plateau (excluding a `1 s` buffer at the rising and falling edges):

* **Channel**: Isotope name.
* **Mean Intensity**: Average simulated intensity (`CPS`) across the steady-state plateau.
* **RSD (%)**: Relative Standard Deviation of the simulated signal across the plateau (reflecting any beat-frequency aliasing or incomplete washout mixing).
* **Ideal Match (%)**: Agreement between the discrete mass spectrometer sampling and the continuous theoretical pulse train integral.
* **Max Error (%)**: Maximum single-cycle deviation from the steady-state mean.
* **Duty Cycle (%)**: Percentage of the total acquisition cycle spent actively measuring that channel.

---

## Modelling Details

* **Washout Scaling**: When using a **Real Composite Peak** from the SPR tab, the pulse profile is dynamically scaled so its `10%` or `1%` width matches the active **Washout (ms)** parameter in the simulator.
* **Sequential vs. Simultaneous Sampling**: For sequential instruments (`Quadrupole`, `Sector-Field`), the simulator models the exact time-staggered dwell windows and settling overheads as the mass analyser steps from isotope to isotope across each cycle. For simultaneous instruments (`TOF`, `Multi-Collector`), all channels integrate across the same synchronised window.
* **Adaptive Simulation Resolution (`dt`)**: Automatically selects a sub-millisecond time-grid step (`dt`) fine enough to resolve narrow transient pulses and short dwell windows without numerical aliasing.
