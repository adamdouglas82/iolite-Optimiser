# iolite Optimiser - Session History & Handoff Notes

This file tracks technical notes, design decisions, and environment quirks from all previous development sessions to ensure future assistants maintain continuity.

---

## Session: Pulse Train Simulator Legend and Scaling Alignment (July 2026)

* **Conversation ID**: `9691bbdc-1788-4c84-8ebf-e5a1a1ae9d99`
* **Key Tasks & Technical Quirks**:
  - **Feature**: Aligned the **Pulse Train Simulator** (Pulse Logger) plot legend click interactions and dynamic row-aware margin scaling with the main **Optimiser** plot.
  - **Unified picker (`on_pick`)**: Connected `pulse_canvas` picker events to the unified `on_pick` handler. Supports single-click toggles, double-click isolations, double-click background restores, and auto-rescale sync.
  - **Dynamic Margin Scaling**: Updated `_on_plot_resize` and `_update_smart_margins` to query `last_pulse_rows` when target canvas is `pulse_canvas`. Automatically resizes (scales down) the plot to fit the legend as items/rows grow.
  - **Simulator Legend Reconstruction**: Rebuilt interactive legend mapping in `update_pulse_plot`. Sets pickers on both legend lines and text labels and populates `pulse_map_legend_to_line`. Tracks hidden channels in `pulse_hidden_channels` to preserve selections.
  - **ZeroDivisionError Fix**: Fixed ZeroDivisionError by setting custom min dwell time spinner to non-zero range `[0.0001, 100]` and protecting division by zero with `max(1e-9, ...)` and `max(1e-6, ...)` guards inside `calculate_constrained_at`.
  - **Dynamic Settings Dialog Sizing**: Replaced `self.setFixedSize(480, 500)` with `self.setFixedWidth(480)` in SettingsDialog, and triggered `self.settings_dlg.adjustSize()` on visibility changes. This expands/shrinks the window height dynamically when custom controls appear.
  - **Cleanup**: Removed the redundant `on_pulse_pick` method.

---

## Session: Adding GitHub Repo and Report Issue Buttons (July 2026)

* **Conversation ID**: `9691bbdc-1788-4c84-8ebf-e5a1a1ae9d99`
* **Key Tasks & Technical Quirks**:
  - **Feature**: Added **GitHub Repo** and **Report Issue** buttons to the bottom of the Hardware Configuration (Settings) dialog.
  - **Layout & Sizing**: Restored the `SettingsDialog` width to its original `480` (height `500`). Implemented a two-row stacked layout at the bottom: Row 1 containing the GitHub utility buttons, and Row 2 containing the version label, updater, and close buttons. Set all four buttons to a matching fixed width of `120px` for symmetric alignment.
  - **Qt URL Handling**: Used native Qt `QDesktopServices.openUrl(QUrl(url))` for browser redirection to align with iolite's embedded Python 3.13 event loop safety.
  - **Focus & Aesthetics**: Set `setFocusPolicy(Qt.NoFocus)` on all three utility buttons to prevent distracting automatic focus shifts/highlights (e.g., when the check updates button disables/enables or when buttons are clicked).
  - **Instrument & Laser Specs**:
    - Updated Thermo iCAP series (Q, RQ, TQ, MSX, MTX) dwell time precision from `0.001` to `0.1` to match actual hardware capability.
    - Increased `Custom Laser` maximum repetition rate from `100` to `300` Hz.
    - Added `"Custom Laser"` as an available source option for `"imageGEO"` and `"ESL193UC"` laser platform presets.
  - **PreScan Parameter Suggestions**:
    - Aligned the Suggested PreScan Parameters window (`show_prescan_suggestions`) to respect instrument constraints: it now checks and snaps recommendations to `allowed_dwells`, `allowed_rr`, `laser_rr_prec` step sizes, and caps speed at `max_speed` (displaying a `(capped at max)` label in the UI when active).
    - Replaced raw property access on spin boxes with the safe retrieval helper `self._get(...)`.
  - **Adaptive SNR formatting**:
    - Addressed rounding errors where low SNR was reported as `0.0` but grew post-optimisation. Increased internal rounding resolution of Initial SNR and Sigma Sep in output rows from `2` to `6` decimal places.
    - Implemented adaptive formatting in the table display layer: shows values `< 0.001` if positive and below threshold, uses 3 decimals below `0.1`, 2 decimals below `1.0`, and 1 decimal otherwise.
  - **Custom Spin Box Inputs**:
    - Fixed an issue where the custom minimum dwell time and rep-rate resolution input boxes automatically appended trailing zeros while typing. Removed the redundant dynamic `setDecimals` calls from the UI update function (`_update_ui_precisions`) so that values are not reformatted on active text edits.
  - **Manuals & Compilation**: Documented the additions in `docs/Optimiser_Manual.md` and recompiled all documentation PDFs.

---

## Session: Adding Check For Updates (July 2026)

* **Conversation ID**: `23c10865-3d71-495e-b153-879793992ed7`
* **Key Tasks & Technical Quirks**:
  * **Feature**: Added a "Check for Updates" button to the Settings dialog that queries the GitHub repository release tag.
  * **PythonQt Threading Quirk**: Discovered that PythonQt crashes if `QTimer.singleShot` is called from a non-Qt thread. Implemented a background python thread `threading.Thread` that writes to a thread-safe `result_holder` list, with a main-thread polling `QTimer` checking for updates every 100ms.
  * **PythonQt Base Class Delegation**: Fixed `CopyableTableWidget.keyPressEvent` crash when using Python's standard `super()`. Explicitly delegated to the parent class: `QTableWidget.keyPressEvent(self, event)`.
  * **OS Shell/Fork Access Violation**: Replaced standard library `webbrowser.open()` with Qt's native `QDesktopServices.openUrl(QUrl(url))` because standard Python webbrowser calls caused access violation crashes (`0xC0000005`) in the embedded Python 3.13 interpreter inside the Qt event loop on Windows.
  * **UI & Precision Fixes**:
    * Preserved sub-millisecond precision by setting decimals before the value inside `DataConfigDialog`.
    * Corrected value property retrieval error on spinboxes (`_get(self.spin_spr_prom, 'value')` instead of directly calling `.value()`).
    * Prevented UI hangs in the Pulse Train Simulator when switching to a Multi-Collector instrument.
    * Improved cross-platform file importing by tracking file count changes and resolving target paths with matching basenames.

---

## Session: Pulse Train Simulator Improvements (July 2026)

* **Conversation ID**: `73854a79-ab52-410e-9480-2a7925cb17d8`
* **Key Tasks & Technical Quirks**:
  * **Advisor Updates**: Bypassed steady-state oversampling stability and RSD warnings for **Combined Integer Sync** (or any strategy that yields perfect integer synchronization).
  * **Status Line Cleanup**: Cleaned up UI status warning clutter by removing the `"Optimised to match Dwell Time..."` warning from the UI status notes list.
  * **Documentation Updates**:
    * Created `docs/Pulse_Train_Simulator_Manual.md`.
    * Updated `docs/Optimiser_Manual.md` (detailed PreScan suggestions, Sync Strategies, and CSV imports).
    * Updated `docs/SPR_Manual.md` (stacked layouts, timestamp details).
    * Updated `docs/Optimisation_Workflow_Manual.md` (PreScan phase insertion).

---

## Session: SPR Tab Manual Creation (Jan 2026)

* **Conversation ID**: `c305e18f-d7d4-49fd-aeed-aa867cd91b2d`
* **Key Tasks**:
  * Researched and documented key SPR UI elements, inputs, outputs, and workflow steps in `docs/SPR_Manual.md`.

---

## Session: Localisation to British English (Jan 2026)

* **Conversation ID**: `503996db-c15f-4786-be53-103c300c3532`
* **Key Tasks**:
  * Localised internal naming conventions and user-facing labels to British English (e.g., `run_optimization` -> `run_optimisation`, `analyze_washout_peaks` -> `analyse_washout_peaks`, `normalize_rgb` -> `normalise_rgb`).

---

## Session: Minimum SNR Logic Fix (Jan 2026)

* **Conversation ID**: `b0d3aca1-e874-4d39-92ae-8b343907d3e8`
* **Key Tasks & Technical Quirks**:
  * Corrected Minimum SNR logic for Multi-Collector (MC) and Time-of-Flight (TOF) systems.
  * Ensured that instead of modifying the dwell time dynamically (since these systems use fixed integrations/simultaneous reads), the UI flags low-SNR channels and displays warning messages without changing the dwell budget itself.

---

## Session: Composite Peak Feature for SPR (Jan 2026)

* **Conversation ID**: `106beadb-39f5-45f1-9a4c-75956d9491bb`
* **Key Tasks & Technical Quirks**:
  * Implemented average peak shape visualization ("Composite Peak") in the Single Pulse Response tab.
  * **Baseline-Split Dynamic Padding**: To capture baselines cleanly before/after each peak, segment extraction windows are calculated by dynamically splitting the gap between successive laser pulses.

---

## Session: SPR Results Alignment & UI Refinements (Jan 2026)

* **Conversation ID**: `f43581b5-f4e9-443b-aeb0-e2ab92d31e08`
* **Key Tasks & Technical Quirks**:
  * Aligned SPR statistics panel values (Average, Max, RSD, Area) into a strict four-column grid with standardized 8pt font and explicit spacing to prevent layout shifting.
  * Added conditional summation for raw counts data in `Logic.analyse_washout_peaks`.
  * Restored undefined variables to fix a `NameError` crash in `run_spr_analysis`.

---

## Session: Intelligent Unit Awareness & Theme Controls (Dec 2025)

* **Conversation ID**: `2de3e3de-040c-4346-978c-83a9d93e651a`
* **Key Tasks**:
  * Added plot control options (Theme selector, Pan/Zoom Y axis, Auto-Rescale Y) directly above the Matplotlib canvas in the SPR tab.
  * Added intelligent unit awareness: when raw CPS data is loaded, Peak Area is integrated and formatted as "Counts" instead of static "cts*s".

---

## Session: Refinement & Code Cleanups (Dec 2025)

* **Conversation ID**: `9468776b-9a36-4194-aa98-34fc712447eb`
* **Key Tasks**:
  * Added multi-mode selection ("Spot", "Line", "Imaging") to the Method Optimiser.
  * Set persistent settings to write to JSON configuration only on window close rather than on every spinbox change to prevent I/O disk throttling.
  * Ensured custom color palettes apply instantly on startup.
