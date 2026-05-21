* 2026-05-21: **v0.9.13** - GPLv3 Release & Pulse Train Simulator Integration
  - **Licensing & Citations**: Codebase licensed under GNU GPL v3.0, introducing official license files, copyright headers, and Section 7 additional terms for scientific citations (citing Tanner 2010, Lockwood 2024, etc.).
  - **Pulse Train Simulator**: Implemented sequential sampling, SNR/CPS-based intensity scaling, automatic debounce triggering, and a signal-only comparison mode (overlay normalization, cropping, and padding).
  - **RSD & Advanced Metrics**: Added a simulator results table calculating RSD, Ideal Match, Max Error, and Duty Cycle statistics, alongside a "Reset to Optimum" shortcut.
  - **Axis Interactions**: Upgraded Matplotlib charts to support targeted, relative zooming and panning by scrolling/dragging directly on axis labels.
  - **Release Automation**: Created GitHub Actions release workflow (Node 24 native) and a pre-commit hook to automate packaging. Modernized Playwright PDF compiler with MathJax support for LaTeX math equations.
* 2026-02-26: **v0.9.11** - Fixed Repetition Rate scaling and Avoid Gaps tolerances.
* 2026-02-23: **v0.9.10** - Added Dosage separation, TOF Washout Margin percentage buffer, sub-1ms precision, and fixed Simultaneous vs Sequential Dwell Budget constraints for TOF. Documentation fully updated.
* 2026-02-18: **v0.9.9** - Added 'Maximum Observed SPR' panel to determine absolute worst-case washout times and fixed QComboBox count property access.
* 2026-02-17: **v0.9.8** - Added 'Even' dwell distribution mode for perfect division without rounding drift, fixed dosage error reporting, and improved rep-rate rounding stability.
* 2026-02-16: **v0.9.7** - Consolidate calculation logic, centralize "Net Measurement Time" (Dwell Budget) directly in the synchronization helper for Multi-Collector/TOF systems. Restored Drift Redistribution for sequential systems. Enhanced UI layout with scroll areas and enforced minimum heights. Updated notation to $L_c$.
* 2026-01-26: **v0.9.6** - Fixed empty database error by passing parameter `has_data=False`, improved empty data rendering on charts, and disabled chart interactions on empty DataFrames.
* 2026-01-22: **v0.9.5** - Fixed Refresh Error startup crash by handling tuple returns from Logic.analyse_washout_peaks(); Improved SPR chart rendering to show raw signal when no peaks detected; Fixed TypeError regression ("str object is not callable").
* 2026-01-12: v0.9.4 - New Data Import Dialog, Channel Presets, 10-Sigma Autodetection, and SPR Scaling Refinements
* 2026-01-08: Enhance Optimiser/SPR interactions, plot layout, and manuals
* 2026-01-08: refactor: Rename to iolite Optimiser and enhance calculation logic
* 2026-01-06: Initial Release: iolite Optimiser
