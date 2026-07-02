#/ Name: iolite Optimiser
#/ Description: Characterises your system's single-pulse washout to calculate the optimal balance of spot size, scan speed, and repetition rate that achieves your target data quality.
#/ Type: UI
#/ Authors: Adam Douglas
#/ Version: dev
#/ Contact: Adam.Douglas@icpms.com

# ==============================================================================
# COPYRIGHT AND GPLv3 LICENSE NOTICE
# ==============================================================================
# Copyright (C) 2026 [Adam Douglas <adamdouglas82@gmail.com>]
#
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <https://gnu.org>.
# ==============================================================================

# ==============================================================================
# SCIENTIFIC CITATION REQUIREMENT (GPLv3 Section 7 Additional Terms)
# ==============================================================================
# Under Section 7 of the GNU General Public License v3.0, you must preserve the 
# author attributions and literature citation requirements stipulated below. 
#
# If you use this script, or elements of the mathematical code/logic herein, 
# to process data for any peer-reviewed publication, abstract, or presentation, 
# you are required to cite the following primary foundational literature:
# 
#     - Tanner, S. D. (2010). "Shorter signals for improved signal to noise ratio, the influence of Poisson distribution." Journal of Analytical Atomic Spectrometry (JAAS), 25, 405–407.
#     - Donard, A., et al. (2015). "Determination of relative rare earth element distributions in very small quantities of uranium ore concentrates using femtosecond UV laser ablation– SF-ICP-MS coupling." J. Anal. At. Spectrom., 30, 2420–2428.
#     - Van Malderen, S. J., et al. (2018). "Considerations on data acquisition in laser ablation-inductively coupled plasma-mass spectrometry with low-dispersion interfaces." Spectrochimica Acta Part B, 140, 29–34.
#     - Van Elteren, J. T., et al. (2019). "Insights into the selection of 2D LA-ICP-MS (multi) elemental mapping conditions." J. Anal. At. Spectrom., 34, 1919.
#     - Ulianov, A., et al. (2015). "The ICPMS signal as a Poisson process: a review of basic concepts." J. Anal. At. Spectrom., 30, 1297–1321.
#     - Currie, L. A. (1968). "Limits for qualitative detection and quantitative determination: Application to Radiochemistry." Anal. Chem., 40, 586-593.
#     - Currie, L. A. (1972). "The Measurement of Environmental Levels of Rare Gas Nuclides and the Treatment of Very Low-Level Counting Data." IEEE Trans. Nucl. Sci., NS19, (1), 119-126.
#     - Lockwood, E. T. (2024). "Multiplexed elemental bioimaging with quadrupole ICP-MS and high-frequency laser ablation systems" J. Anal. At. Spectrom., 39, 1125
# ==============================================================================

# pyright: reportMissingImports=false, reportMissingModuleSource=false

import math
import os
import json
import traceback
import pandas as pd
import numpy as np
import random
import colorsys

try:
    from scipy.signal import find_peaks, peak_widths
except ImportError:
    find_peaks = None
    peak_widths = None

try:
    import matplotlib
    import matplotlib.pyplot as plt
    
    # Try unified QtAgg backend (supports both Qt 5 & Qt 6)
    try:
        matplotlib.use('qtagg')
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    except Exception:
        matplotlib.use('Qt5Agg')
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator, EngFormatter, ScalarFormatter
    from cycler import cycler
    from matplotlib.lines import Line2D
    import matplotlib.colors as mcolors
except Exception as e:
    print(f"Import Error: {e}")
    print(traceback.format_exc())


# iolite-specific imports (PythonQt)
from iolite.QtGui import (QAction, QSizePolicy, QWidget, QLabel, QDoubleSpinBox, QSpinBox, QCheckBox, 
                          QComboBox, QGridLayout, QHBoxLayout, QVBoxLayout, QGroupBox, QFormLayout, 
                          QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QMenu, 
                          QColorDialog, QDialog, QPushButton, QScrollArea, QSplitter, QFrame, 
                          QLineEdit, QTabWidget, QStackedWidget, QApplication, QPalette, QColor,
                          QListWidget, QListWidgetItem, QTextEdit, QInputDialog, QMessageBox,
                          QFileDialog)
from iolite.QtCore import Qt, QTimer, QSize, QEvent, QUrl

try:
    data
except NameError:
    from iolite import data

try:
    IoLog
except NameError:
    from iolite import IoLog

# Backward-compatible alias for trapezoidal integration (NumPy 1.x vs 2.x)
trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))

# --- FEATURE ENABLE ---
SHOW_PULSE_TRAIN_SIMULATOR = True  # Currently Beta

# --- EMBEDDED CONSTANTS ---
VERSION = "dev"

ICP_SPECS = {
    "Agilent": {
        "7900": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.1},
        "7700": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.1},
        "8800": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.1},
        "8900": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.1},
        "7850": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.1},
    },
    "Thermo": {
        "iCAP Q":   {"type": "Quadrupole", "min_dwell": 1.0, "prec": 0.001},
        "iCAP RQ":  {"type": "Quadrupole", "min_dwell": 1.0, "prec": 0.001},
        "iCAP TQ":  {"type": "Quadrupole", "min_dwell": 1.0, "prec": 0.001},
        "iCAP MSX": {"type": "Quadrupole", "min_dwell": 0.5, "prec": 0.001},
        "iCAP MTX": {"type": "Quadrupole", "min_dwell": 0.5, "prec": 0.001},
        "Element 2":    {"type": "Sector-Field", "min_dwell": 2.0, "prec": 1.0},
        "Element 2 XR": {"type": "Sector-Field", "min_dwell": 0.1, "prec": 0.1},
        "Neoma":   {"type": "Multi-Collector", "min_dwell": 0.0, "prec": 0.001, "allowed_dwells": [50, 100, 250, 500, 1000]},
        "Neptune": {"type": "Multi-Collector", "min_dwell": 0.0, "prec": 0.001, "allowed_dwells": [66, 131, 262, 524, 1049]},
    },
    "Perkin Elmer": {
        "5000": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.001},
        "2000": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.001},
        "1000": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.001},
    },
    "TOFWERK": {
        "ICPTOF": {"type": "TOF", "min_dwell": 0.1, "prec": 0.001},
    },
    "Nu Instruments": {
        "Vitesse": {"type": "TOF", "min_dwell": 0.1, "prec": 0.001},
        "Sapphire": {"type": "Multi-Collector", "min_dwell": 0.0, "prec": 0.001, "allowed_dwells": [100, 200, 500, 1000, 2000, 3000]},
    },

    "Custom": {
        "Custom Model": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.1}
    }
}


MFR_DEFAULT_UNITS = {
    "Agilent": "s",
    "Thermo": "s",
    "Perkin Elmer": "ms",
    "TOFWERK": "ms",
    "Nu Instruments": "ms",
    "Custom": "ms"
}

TYPE_TO_TECH_MAP = {
    "Quadrupole": "Quadrupole",
    "Sector-Field": "Sector-Field",
    "Multi-Collector": "Multi-Collector",
    "TOF": "TOF"
}

STAGE_SPECS = {"TV2": 2000, "TV3": 20000, "Custom Stage": 2000}

LASER_SOURCES = {
    "MLase 1000":   {"max_rr": 1000, "allowed_rr": None, "rr_prec": 1},
    "Coherent 500": {"max_rr": 500,  "allowed_rr": None, "rr_prec": 1},
    "Polaris":      {"max_rr": 100,  "allowed_rr": [1, 2, 4, 5, 10, 20, 50, 100], "rr_prec": 1},
    "Wedge":        {"max_rr": 1000, "allowed_rr": None, "rr_prec": 1},
    "Coherent 200": {"max_rr": 200,  "allowed_rr": None, "rr_prec": 1},
    "Tempest 213":  {"max_rr": 20,   "allowed_rr": [1, 2, 4, 5, 10, 20], "rr_prec": 1},
    "Pharos":       {"max_rr": 1000, "allowed_rr": None, "rr_prec": 1},
    "Tempest 266":  {"max_rr": 10,   "allowed_rr": [1, 2, 4, 5, 10], "rr_prec": 1},
    "Custom Laser": {"max_rr": 100,  "allowed_rr": None, "rr_prec": 1}
}

LASER_PLATFORMS = {
    "imageGEO": {
        "stages": ["TV2", "TV3"], 
        "lasers": ["MLase 1000", "Coherent 500"]
    },
    "imageBIO": {
        "stages": ["TV2", "TV3"], 
        "lasers": ["Polaris", "Wedge"]
    },
    "ESL193UC": {
        "stages": ["TV2"], 
        "lasers": ["Coherent 200"]
    },
    "ESL213": {
        "stages": ["TV2", "TV3"], 
        "lasers": ["Tempest 213"]
    },
    "ESLfemto": {
        "stages": ["TV2", "TV3"], 
        "lasers": ["Pharos"]
    },
    "ESL266": {
        "stages": ["TV2"], 
        "lasers": ["Tempest 266"]
    },
    "Custom": {
        "stages": ["TV2", "TV3", "Custom Stage"],
        "lasers": sorted(list(LASER_SOURCES.keys()))
    }
}

PLOT_COLORS = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
    "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
    "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
    "#17becf", "#9edae5"
]

def get_settings_dir():
    try:
        home = os.path.expanduser("~")
        
        # Determine the official iolite 4.11 AppData path
        if os.name == 'nt' or 'LOCALAPPDATA' in os.environ:
            # Windows: AppData/Local/iolite-software/iolite4/iolite Optimiser
            app_data = os.environ.get('LOCALAPPDATA', os.path.join(home, 'AppData', 'Local'))
            new_dir = os.path.join(app_data, 'iolite-software', 'iolite4', 'iolite Optimiser')
        else:
            # macOS: Library/Application Support/iolite-software/iolite4/iolite Optimiser
            new_dir = os.path.join(home, 'Library', 'Application Support', 'iolite-software', 'iolite4', 'iolite Optimiser')

        # Setup the old paths for migration checks
        old_dir = os.path.join(home, "Documents", "iolite", "iolite Optimiser")
        old_alt_dir = os.path.join(home, "OneDrive", "Documents", "iolite", "iolite Optimiser")

        # Create the new directory if it doesn't exist
        os.makedirs(new_dir, exist_ok=True)

        # Migrate files if needed
        for filename in ["iolite Optimiser Settings.json", "iolite Optimiser Channel Presets.json"]:
            new_file = os.path.join(new_dir, filename)
            if not os.path.exists(new_file):
                for src_dir in [old_dir, old_alt_dir]:
                    src_file = os.path.join(src_dir, filename)
                    if os.path.exists(src_file):
                        try:
                            import shutil
                            shutil.copy2(src_file, new_file)
                            break
                        except Exception:
                            pass
        return new_dir
    except Exception:
        try:
            home = os.path.expanduser("~")
            fallback = os.path.join(home, "Documents", "iolite", "iolite Optimiser")
            os.makedirs(fallback, exist_ok=True)
            return fallback
        except Exception:
            return ""

# --- EMBEDDED LOGIC ---
class Logic:
    @staticmethod
    def calculate_lockwood_rsd(f, W_ms, dt_ms):
        """
        Calculates the predicted %RSD of aliasing ripple based on Thomas Lockwood (2024) Eq 1-3.
        f: repetition rate in Hz
        W_ms: full-width washout time in milliseconds
        dt_ms: analyte dwell time in milliseconds
        """
        if f <= 0 or W_ms <= 0 or dt_ms <= 0:
            return 0.0
            
        W_s = W_ms / 1000.0
        dt_s = dt_ms / 1000.0
        T = f * dt_s
        T_floor = math.floor(T)
        T_ceil = math.ceil(T)
        
        # Guard against exact integers to avoid float round-off issues
        if abs(T - round(T)) < 1e-9:
            return 0.0
            
        A_min = ((T - T_floor) ** 2 - T_floor) / 2.0
        A_max = (T_ceil - (T_ceil - T) ** 2) / 2.0
        
        if A_max <= 0:
            return 0.0
            
        # Peak-to-peak ripple percentage
        ripple_pct = (100.0 / (f * W_s)) * (abs(A_max - A_min) / A_max)
        
        # Convert peak-to-peak ripple of a triangular wave to statistical RSD (standard deviation)
        # For a triangular wave, std_dev = peak_to_peak / sqrt(12)
        statistical_rsd = ripple_pct / math.sqrt(12.0)
        return statistical_rsd

    @staticmethod
    def get_lockwood_min_dwell(f, W_ms, target_rsd=5.0):
        """
        Sweeps cycle times starting from the washout time (W_ms) upwards
        to find the minimum time that achieves a predicted Lockwood aliasing RSD below the target_rsd.
        """
        if f <= 0 or W_ms <= 0:
            return 2.0  # Fallback to standard 2.0 ms
            
        # Sweep from W_ms upwards in steps of 0.1 ms up to W_ms + 100 ms
        start_steps = int(round(W_ms * 10))
        for steps in range(start_steps, start_steps + 1000):
            dt_ms = steps / 10.0
            rsd = Logic.calculate_lockwood_rsd(f, W_ms, dt_ms)
            if rsd <= target_rsd:
                return dt_ms
        return W_ms  # Safe fallback if never reached

    @staticmethod
    def analyse_washout_peaks(df, isotope, prominence_threshold=100.0, min_distance=1, is_cps=True, wlen=None, apply_smoothing=False, smooth_window=5, bg_sub=False, bg_start_s=2.0, bg_end_s=8.0):
        """
        Analyses washout peaks (Single Pulse Response) for a given isotope.
        Calculates FW0.1M (10% Height) and FW0.01M (1% Height).
        """
        if find_peaks is None:
            return pd.DataFrame(), {"Error": "Scipy (scipy.signal) is not installed."}

        if isotope not in df.columns or "Time" not in df.columns:
            return pd.DataFrame(), {"Error": f"Columns missing: {isotope} or Time"}

        y_proc = df[isotope].values
        time_proc = df['Time'].values

        t0 = time_proc[0] if len(time_proc) > 0 else 0
        t_zeroed = time_proc - t0 if len(time_proc) > 0 else time_proc

        # --- BACKGROUND SUBTRACTION ---
        if bg_sub:
            bg_mask = (t_zeroed >= bg_start_s) & (t_zeroed <= bg_end_s)
            if bg_mask.any():
                bg_mean = np.mean(y_proc[bg_mask])
                y_proc = np.clip(y_proc - bg_mean, 0, None)

        # Calculate typical sample rate (dt)
        dt = np.median(np.diff(time_proc)) if len(time_proc) > 1 else 1.0

        if apply_smoothing and smooth_window > 0:
            smooth_window_pts = max(3, int(round(smooth_window / dt)))
            y_proc_series = pd.Series(y_proc)
            y_proc = y_proc_series.rolling(window=smooth_window_pts, center=True, min_periods=1).median().values

        # 1. Find Peaks
        try:
            min_distance_pts = max(1, int(round(min_distance / dt)))
            wlen_pts = max(1, int(round(wlen / dt))) if wlen is not None and wlen > 0 else None
            peaks, properties = find_peaks(y_proc, prominence=prominence_threshold, distance=min_distance_pts, wlen=wlen_pts)
        except Exception as e:
            return pd.DataFrame(), {"Error": str(e)}
        
        if len(peaks) == 0:
            return pd.DataFrame(), {"Count": 0}
            
        # 2. Calculate Widths
        widths_01, width_heights_01, left_ips_01, right_ips_01 = peak_widths(y_proc, peaks, rel_height=0.9)
        widths_001, width_heights_001, left_ips_001, right_ips_001 = peak_widths(y_proc, peaks, rel_height=0.99)
        
        peak_data = []
        
        for i, p_idx in enumerate(peaks):
            # Simple interpolation helper
            def get_time(idx):
                idx_int = int(idx)
                idx_frac = idx - idx_int
                if idx_int + 1 < len(time_proc):
                    return time_proc[idx_int] + (time_proc[idx_int+1] - time_proc[idx_int]) * idx_frac
                return time_proc[-1]

            t_01_left = get_time(left_ips_01[i])
            t_01_right = get_time(right_ips_01[i])
            w_01 = t_01_right - t_01_left
            
            t_001_left = get_time(left_ips_001[i])
            t_001_right = get_time(right_ips_001[i])
            w_001 = t_001_right - t_001_left
            
            # Calculate Area 1%
            idx_start_001 = max(0, int(np.floor(left_ips_001[i])))
            idx_end_001 = min(len(time_proc), int(np.ceil(right_ips_001[i])))
            if idx_end_001 > idx_start_001:
                if is_cps:
                    area_001 = trapz(y_proc[idx_start_001:idx_end_001], x=time_proc[idx_start_001:idx_end_001])
                    area_01 = trapz(y_proc[max(0, int(np.floor(left_ips_01[i]))):min(len(time_proc), int(np.ceil(right_ips_01[i])))], 
                                    x=time_proc[max(0, int(np.floor(left_ips_01[i]))):min(len(time_proc), int(np.ceil(right_ips_01[i])))])
                else:
                    area_001 = np.sum(y_proc[idx_start_001:idx_end_001])
                    area_01 = np.sum(y_proc[max(0, int(np.floor(left_ips_01[i]))):min(len(time_proc), int(np.ceil(right_ips_01[i])))])
            else:
                area_001 = 0.0
                area_01 = 0.0

            peak_data.append({
                'Peak Index': i + 1,
                'Peak Time (s)': time_proc[p_idx],
                'Max Intensity': y_proc[p_idx],
                'FW0.1M (s)': w_01,
                'FW0.01M (s)': w_001,
                'Area 10%': area_01,
                'Area 1%': area_001,
                'p_idx': p_idx,
                'l01': t_01_left, 'r01': t_01_right, 'h01': width_heights_01[i],
                'l001': t_001_left, 'r001': t_001_right, 'h001': width_heights_001[i]
            })
            
        results_df = pd.DataFrame(peak_data)
        return results_df

    @staticmethod
    def summarize_peaks(results_df):
        """
        Calculates summary statistics (Mean, RSD, Max) for a DataFrame of detected SPR peaks.
        """
        if results_df.empty:
            return {"Count": 0}

        stats = {
            'Count': len(results_df),
            'FW0.1M Mean': results_df['FW0.1M (s)'].mean(),
            'FW0.1M RSD': (results_df['FW0.1M (s)'].std() / results_df['FW0.1M (s)'].mean() * 100) if len(results_df) > 1 else 0.0,
            'FW0.01M Mean': results_df['FW0.01M (s)'].mean(),
            'FW0.01M RSD': (results_df['FW0.01M (s)'].std() / results_df['FW0.01M (s)'].mean() * 100) if len(results_df) > 1 else 0.0,
            'FW0.1M Max': results_df['FW0.1M (s)'].max(),
            'FW0.01M Max': results_df['FW0.01M (s)'].max(),
            'Area 10% Mean': results_df['Area 10%'].mean(),
            'Area 10% RSD': (results_df['Area 10%'].std() / results_df['Area 10%'].mean() * 100) if len(results_df) > 1 else 0.0,
            'Area 1% Mean': results_df['Area 1%'].mean(),
            'Area 1% RSD': (results_df['Area 1%'].std() / results_df['Area 1%'].mean() * 100) if len(results_df) > 1 else 0.0
        }
        return stats

    @staticmethod
    def generate_composite_peak(df, isotope, peaks_df, bg_sub=False, bg_start_s=2.0, bg_end_s=8.0):
        """
        Generates a composite (averaged) peak shape by aligning and normalizing detected peaks.
        Uses asymmetric windows that split the baseline gap between adjacent pulses.
        Aligned such that Relative Time 0 is at the Maximum Intensity Time (Peak Time (s)).
        """
        if peaks_df.empty or isotope not in df.columns or "Time" not in df.columns:
            return None

        y_raw = df[isotope].values.copy()
        t_raw = df['Time'].values
        
        if bg_sub:
            t0 = t_raw[0] if len(t_raw) > 0 else 0
            t_zeroed = t_raw - t0 if len(t_raw) > 0 else t_raw
            bg_mask = (t_zeroed >= bg_start_s) & (t_zeroed <= bg_end_s)
            if bg_mask.any():
                bg_mean = np.mean(y_raw[bg_mask])
                y_raw = np.clip(y_raw - bg_mean, 0, None)
        
        # Calculate typical sample rate
        dt = np.median(np.diff(t_raw)) if len(t_raw) > 1 else 1.0
        
        # 1. Calculate the core pulse dimensions
        # Alignment is at Peak Time (s). Use robust percentiles instead of max to avoid outliers pulling in other peaks.
        dist_left = peaks_df['Peak Time (s)'] - peaks_df['l001']
        dist_right = peaks_df['r001'] - peaks_df['Peak Time (s)']
        
        core_left = np.percentile(dist_left, 90)
        core_right = np.percentile(dist_right, 90)
        core_fw = core_left + core_right
        
        # 2. Asymmetric baseline buffers
        if len(peaks_df) > 1:
            pdf = peaks_df.sort_values('Peak Time (s)')
            gaps = pdf['l001'].values[1:] - pdf['r001'].values[:-1]
            median_gap = max(0, np.median(gaps))
            p2p = np.median(np.diff(pdf['Peak Time (s)'].values))
            
            # 10 % Margin Left, 20 % Margin 
            b_left = core_left + (0.1 * core_fw)
            b_right = core_right + (0.20 * core_fw)
            
            # Safety 1: Ensure we don't clip adjacent peaks via gap allocation
            b_left = min(b_left, core_left + (0.4 * median_gap))
            b_right = min(b_right, core_right + (0.4 * median_gap))
            
            # Safety 2: Hard limit to never exceed halfway to the next peak
            b_left = min(b_left, 0.45 * p2p)
            b_right = min(b_right, 0.45 * p2p)
        else:
            # Single peak fallbacks
            b_left = core_left + (0.1 * core_fw)
            b_right = core_right + (0.20 * core_fw)
            
        # Global minimums (0.1ms / 0.1ms) - Fixes "needle" peaks
        b_left = max(b_left, 0.0001)
        b_right = max(b_right, 0.0001)
        
        peak_segments = []
        
        for i, row in peaks_df.iterrows():
            # Align at Peak Time (s) absolute time
            t_start_abs = row['Peak Time (s)'] - b_left
            t_end_abs = row['Peak Time (s)'] + b_right
            
            # Convert to indices
            idx_start = np.searchsorted(t_raw, t_start_abs)
            idx_end = np.searchsorted(t_raw, t_end_abs, side='right')
            
            if idx_end <= idx_start: continue
            
            seg_y = y_raw[idx_start:idx_end]
            seg_t = t_raw[idx_start:idx_end] - row['Peak Time (s)'] # Center t=0 at Peak Time (s)
            
            # Normalise intensity
            mx = row['Max Intensity']
            if mx > 0:
                seg_y_norm = seg_y / mx
                peak_segments.append((seg_t, seg_y_norm))

        if not peak_segments:
            return None

        # Interpolation to a common grid
        # Grid covers from -b_left to b_right
        grid_t = np.arange(-b_left, b_right + dt, dt)
        sum_y = np.zeros_like(grid_t)
        count_y = np.zeros_like(grid_t)
        
        for seg_t, seg_y in peak_segments:
            interp_y = np.interp(grid_t, seg_t, seg_y, left=0, right=0)
            sum_y += interp_y
            count_y += 1
            
        avg_y = sum_y / count_y if len(peak_segments) > 0 else sum_y
        
        # Shift X-axis so that 0 represents the 1% rise (FW0.01M left edge)
        # This preserves the highest-point alignment between all pulses, but normalizes the visual start time.
        lvl = 0.01
        above_idx = np.where(avg_y >= lvl)[0]
        if len(above_idx) > 0:
            idx_l = above_idx[0]
            if idx_l > 0:
                t_l = np.interp(lvl, [avg_y[idx_l-1], avg_y[idx_l]], [grid_t[idx_l-1], grid_t[idx_l]])
            else:
                t_l = grid_t[idx_l]
            grid_t = grid_t - t_l
        
        return pd.DataFrame({'Relative Time (s)': grid_t, 'Normalised Intensity': avg_y})

    @staticmethod
    def get_isotope_stats(name, df_sig, df_blk, dwell_s, is_raw_counts):
        """Standardized conversion between counts and CPS for isotope statistics."""
        m_sig, m_blk = df_sig[name].mean(), df_blk[name].mean()
        s_blk = df_blk[name].std()
        
        sig_cps = m_sig if not is_raw_counts else m_sig / dwell_s
        blk_cps = m_blk if not is_raw_counts else m_blk / dwell_s
        std_counts = s_blk if is_raw_counts else s_blk * dwell_s
        
        return {
            "m_sig": m_sig, "m_blk": m_blk, "s_blk": s_blk,
            "sig_cps": sig_cps, "blk_cps": blk_cps, "std_counts": std_counts
        }

    @staticmethod
    def calculate_robust_snr(sig_cps, blk_cps, bl_dt, std_counts, snr_threshold=0.0):
        """
        Standardized calculation for 'Robust SNR' (Sigma Separation).
        Centralizes Poisson vs. Observed variance and Critical Level (Lc).
        """
        # 1. Determine the Noise Floor (Count Domain)
        total_bg_counts = blk_cps * bl_dt
        # Theoretical Poisson Variance vs Observed Variance
        effective_noise_var = max(total_bg_counts, std_counts**2)
        
        # 2. Calculate Critical Level (Lc) for Detection (95% Confidence)
        z = 1.645
        if total_bg_counts == 0:
            L_c_counts = (z**2)/2 + z * math.sqrt(2 * 0.4) 
        else:
            L_c_counts = z * math.sqrt(2 * effective_noise_var)

        # 3. Calculate "Robust SNR" (Detection Metric)
        # Factor = z * sqrt(2) ~ 2.326 to restore familiar Sigma scaling
        sigma_factor = z * math.sqrt(2)
        net_counts = (sig_cps - blk_cps) * bl_dt
        
        base_snr = 0
        if L_c_counts > 0:
            detection_ratio = net_counts / L_c_counts
            base_snr = detection_ratio * sigma_factor 

        return {
            'base_snr': base_snr,
            'effective_noise_var': effective_noise_var,
            'is_zero_bg': (total_bg_counts == 0),
            'L_c_counts': L_c_counts
        }

    @staticmethod
    def get_best_integration_time(target_at, overhead_s, allowed_dwells_s, precision_s, tech):
        """
        Centralizes integration time selection for MC/TOF and hardware snapping.
        """
        if tech == "Multi-Collector" and allowed_dwells_s:
            for t in sorted(allowed_dwells_s):
                if (t + overhead_s) >= (target_at - 1e-7): 
                    return t + overhead_s
            return max(allowed_dwells_s) + overhead_s
        
        # Standard snap to hardware precision
        return round(target_at / precision_s) * precision_s
    @staticmethod
    def generate_pulse_train_v3(composite_df, rep_rate_hz, pulses_per_cycle, cycle_time_s, channel_specs, background_s=5.0, signal_s=10.0, dt_s=1e-5):
        """
        Generates simulated pulse train with sequential sampling support.
        
        Args:
            channel_specs: List of dicts [{'name': 'U238', 'dwell': 0.1, 'offset': 0.0}, ...]
        """
        try:
            if composite_df is None or composite_df.empty:
                return None, None
                
            IoLog.information("Logic.generate_pulse_train_v3: Starting...")
                
            # 1. Prepare Composite Pulse (Interpolant)
            comp_t = composite_df['Relative Time (s)'].values
            comp_y = composite_df['Normalised Intensity'].values
            
            # 2. Setup Time Axis
            total_duration = background_s + signal_s + background_s
            t_axis = np.arange(0, total_duration, dt_s)
            y_theoretical = np.zeros_like(t_axis)
            
            # 3. Add Pulses
            # Pulses centered in the signal region
            pulse_period_s = 1.0 / rep_rate_hz if rep_rate_hz > 0 else 1.0
            
            # Number of pulses to simulate
            # Usually Total pulses = RepRate * Signal Duration
            # But the 'pulses' arg in v2 was 'pulses per pixel/cycle'.
            # Here let's just fill the signal_s window with pulses at rep_rate
            n_pulses = int(signal_s * rep_rate_hz)
            
            start_time = background_s
            
            n_total_samples = len(y_theoretical)
            
            for i in range(n_pulses):
                t_pulse = start_time + i * pulse_period_s
                
                # Find the grid indices in t_axis that span the pulse duration
                t_start_abs = t_pulse + comp_t.min()
                t_end_abs = t_pulse + comp_t.max()
                
                idx_start = np.searchsorted(t_axis, t_start_abs)
                idx_end = np.searchsorted(t_axis, t_end_abs, side='right')
                
                if idx_start >= n_total_samples: break
                if idx_end > n_total_samples:
                    idx_end = n_total_samples
                    
                if idx_end > idx_start:
                    t_segment = t_axis[idx_start:idx_end]
                    # Interpolate the pulse shape directly onto the exact grid points of this segment
                    y_segment = np.interp(t_segment, comp_t + t_pulse, comp_y, left=0, right=0)
                    y_theoretical[idx_start:idx_end] += y_segment

            # 4. Integrate (Measured Signal) for EACH Channel
            channel_results = {}
            
            # Pre-compute cumulative sum for fast integration
            y_cumsum = np.cumsum(y_theoretical) * dt_s # Integral
            
            # Ensure cycle time is valid
            if cycle_time_s <= 0: cycle_time_s = 0.1
            
            for spec in channel_specs:
                name = spec['name']
                dwell = spec['dwell']
                offset = spec['offset']
                
                # Generate sample windows
                # Start from t=0 (or -offset?)
                # We want samples to cover total_duration
                
                # N cycles
                n_cycles = int(total_duration / cycle_time_s) + 1
                
                # Measurements
                m_times = []
                m_y = []
                
                for n in range(n_cycles):
                    cycle_start = n * cycle_time_s
                    t_start = cycle_start + offset
                    t_end = t_start + dwell
                    
                    if t_start >= total_duration: break
                    
                    # Integration limits indices
                    idx_s = int(t_start / dt_s)
                    idx_e = int(t_end / dt_s)
                    
                    # Bounds check
                    if idx_s < 0: idx_s = 0
                    if idx_e > len(y_cumsum) - 1: idx_e = len(y_cumsum) - 1
                    
                    if idx_e > idx_s:
                        integral = y_cumsum[idx_e] - y_cumsum[idx_s]
                        # Normalized intensity = Integral / Dwell
                        # Apply Scaling (SNR)
                        scale = spec.get('scale', 1.0)
                        val = (integral / (dwell if dwell > 0 else 1.0)) * scale
                        
                        # Store point at END of Dwell Time (as requested)
                        # User Request: "place channels x data point at the end of it's measured dwell time"
                        m_times.append(t_end) 
                        m_y.append(val)
                
                # Store Dwell in DF for plotting
                dwells = [dwell] * len(m_times)
                channel_results[name] = pd.DataFrame({'Time': m_times, 'Intensity': m_y, 'Dwell': dwells})

            IoLog.information("Logic.generate_pulse_train_v3: Complete")
            
            theo_df = pd.DataFrame({'Time': t_axis, 'Intensity': y_theoretical})
            return theo_df, channel_results
            
        except Exception as e:
            IoLog.error(f"Logic.generate_pulse_train_v3 CRASH: {e}")
            IoLog.error(traceback.format_exc())
            return None, None

    @staticmethod
    def generate_pulse_train_v2(composite_df, rep_rate_hz, pulses, acq_time_s, background_s=5.0, signal_s=10.0, dt_s=1e-5):
        """
        Generates a simulated pulse train based on the composite peak shape.
        Returns:
            theoretical_df: High-res instantaneous signal
            measured_df: Signal integrated over Acquisition Time buckets
        """
        try:
            if composite_df is None or composite_df.empty:
                return None, None
                
            IoLog.information("Logic.generate_pulse_train: Starting...")
                
            # 1. Prepare Composite Pulse (Interpolant)
            comp_t = composite_df['Relative Time (s)'].values
            comp_i = composite_df['Normalised Intensity'].values
            
            # Shift comp_t so it starts at 0 for easier handling
            t_min = comp_t.min()
            comp_t_shifted = comp_t - t_min 
            pulse_duration = comp_t_shifted.max()
            
            # 2. Setup Time Axis
            IoLog.information("Logic.generate_pulse_train: Step 2 - Setup Time Axis")
            # User requested equivalent background at rear
            total_duration = background_s + signal_s + background_s
            t_axis = np.arange(0, total_duration, dt_s)
            y_theoretical = np.zeros_like(t_axis)
            
            # 3. Generate Pulse Train
            IoLog.information("Logic.generate_pulse_train: Step 3 - Generate Pulses")
            # Start firing after background_s
            start_time = background_s
            
            # Calculate number of pulses to fire based on signal duration and rep rate
            if rep_rate_hz > 0:
                period = 1.0 / rep_rate_hz
                num_pulses = int(signal_s * rep_rate_hz)
                
                # Vectorized Pulse Addition
                # We add the pulse shape to the y_theoretical array at calculated start indices
                
                # Interpolate composite to the simulation grid resolution (dt_s) once
                # This works if dt_s is finer than composite resolution, which it should be (10us)
                pulse_grid_t = np.arange(0, pulse_duration + dt_s, dt_s)
                pulse_grid_y = np.interp(pulse_grid_t, comp_t_shifted, comp_i, left=0, right=0)
                
                n_pulse_samples = len(pulse_grid_y)
                n_total_samples = len(y_theoretical)
                
                for i in range(num_pulses):
                    t_fire = start_time + i * period
                    idx_start = int(t_fire / dt_s)
                    idx_end = idx_start + n_pulse_samples
                    
                    if idx_start >= n_total_samples: break
                    
                    # Clip if end extends beyond buffer
                    if idx_end > n_total_samples:
                        valid_len = n_total_samples - idx_start
                        y_theoretical[idx_start:n_total_samples] += pulse_grid_y[:valid_len]
                    else:
                        y_theoretical[idx_start:idx_end] += pulse_grid_y
                        
            # 4. Generate Measured (Integrated) Signal
            IoLog.information("Logic.generate_pulse_train: Step 4 - Generate Measured")
            # Integration buckets of size acq_time_s
            n_buckets = int(total_duration / acq_time_s)
            if n_buckets == 0: n_buckets = 1
            measured_t = np.arange(n_buckets) * acq_time_s
            measured_y = np.zeros(n_buckets)
            
            # Reshape theoretical into chunks of acq_time (approx)
            # Note: This is an approximation. For exact sync logic, we'd need to be careful with sample boundaries.
            # But for visualization, this is sufficient.
            samples_per_bucket = int(acq_time_s / dt_s)
            
            if samples_per_bucket > 0:
                # fast sum using reshape if length matches perfectly, else loop
                # Loop is safer for edge cases
                for i in range(n_buckets):
                    s_idx = i * samples_per_bucket
                    e_idx = s_idx + samples_per_bucket
                    if s_idx >= len(y_theoretical): break
                    # Sum * dt gives Area (Integrated Intensity)
                    # But mass spec usually reports Counts or CPS.
                    # If we assume y_theoretical is "Intensity/sec (CPS equivalent of pulse)", then
                    # this integration gives Counts. 
                    # Let's normalize it so max is somewhat relatable or just return raw integration.
                    chunk = y_theoretical[s_idx:min(e_idx, len(y_theoretical))]
                    measured_y[i] = np.sum(chunk)
            
            IoLog.information("Logic.generate_pulse_train: Step 5 - Return DataFrame")
            return pd.DataFrame({'Time': t_axis, 'Intensity': y_theoretical}), pd.DataFrame({'Time': measured_t, 'Intensity': measured_y})

        except Exception as e:
            IoLog.error(f"Logic.generate_pulse_train CRASH: {e}")
            IoLog.error(traceback.format_exc())
            return None, None

    @staticmethod
    def calculate_constrained_at(inputs):
        """Helper to determine constrained Acquisition Time and synchronized rep rate."""
        fw_s = inputs['washout_ms'] / 1000.0
        bs = inputs['spot_size_um']
        n_val = inputs['pulses_per_pixel']
        n_dose = inputs.get('dosage', n_val)
        rr_max = inputs['max_rr_hz']
        ss_max = inputs['max_speed_um_s']
        is_mc = inputs['icp_technology'] == "Multi-Collector"
        valid_times_s = [d/1000.0 for d in (inputs.get('allowed_dwells') or [])]
        p_s = inputs['precision_ms'] / 1000.0
        strategy = inputs.get('sync_strategy', 'Adaptive Integer Sync (Auto)')
        is_oversampling = "Oversampling" in strategy

        at_from_washout = fw_s
        if is_oversampling:
            at_from_rr = at_from_washout
            at_from_ss = at_from_washout
        else:
            rr_ideal = n_val / fw_s if fw_s > 0 else rr_max
            at_from_rr = n_val / rr_max if rr_ideal > rr_max else at_from_washout
            ss_ideal = (bs * rr_ideal) / n_dose if n_dose > 0 else 0
            at_from_ss = (bs * n_val) / (ss_max * n_dose) if (n_dose > 0 and ss_ideal > ss_max) else at_from_washout

        min_duty = inputs.get('min_duty_cycle', 0)
        overhead_s = inputs.get('overhead_ms', 0) / 1000.0
        min_dwell_req_s = inputs.get('min_dwell_needed_ms', 0) / 1000.0
        
        at_from_duty = at_from_washout
        harmonic = 1
        
        # State tracking for reasons why Acq Time was increased
        at_reasons = []
        if at_from_rr > (at_from_washout + 1e-7): at_reasons.append("Rep Rate Limit")
        if ss_max > 0 and at_from_ss > (at_from_washout + 1e-7): at_reasons.append("Stage Speed Limit")

        # --- HARMONIC SCALING ---
        rr_h1 = max(1, math.floor(n_val / fw_s if fw_s > 0 else rr_max))
        at_h1 = n_val / rr_h1
        base_washout = round(at_h1 / p_s) * p_s 
        base = max(base_washout, at_from_rr, at_from_ss)
        
        at_from_duty_direct = 0
        if min_duty > 0 and min_duty < 1.0:
            min_dwell_for_duty = (min_duty * overhead_s) / (1.0 - min_duty)
            required_cycle_for_duty = min_dwell_for_duty + overhead_s
            if is_mc:
                at_from_duty_direct = required_cycle_for_duty
                if required_cycle_for_duty > (base + 1e-7): at_reasons.append(f"Duty Cycle ({min_duty*100:.0f}%)")
            else:
                calc_harmonic = math.ceil((required_cycle_for_duty / base) - 1e-7)
                if calc_harmonic > 1: at_reasons.append(f"Duty Cycle ({calc_harmonic}x)")
                harmonic = max(harmonic, calc_harmonic)

        at_from_dwell_direct = 0
        if min_dwell_req_s > 0:
            required_cycle = min_dwell_req_s + overhead_s
            if is_mc:
                 at_from_dwell_direct = required_cycle
                 if required_cycle > (base + 1e-7): at_reasons.append(f"Dwell Budget ({min_dwell_req_s*1000:.0f} ms)")
            else:
                 harmonic_dwell = math.ceil((required_cycle / base) - 1e-7)
                 if harmonic_dwell > 1: at_reasons.append(f"Dwell Budget ({harmonic_dwell}x)")
                 harmonic = max(harmonic, harmonic_dwell)

        notes = []
        if at_reasons:
            res_str = at_reasons[0]
            if len(at_reasons) > 1: res_str = ", ".join(at_reasons[:-1]) + " and " + at_reasons[-1]
            notes.append(f"Constraint: Acq Time increased for {res_str}")
            
        at_from_duty = base * harmonic
        at_needed = max(at_from_washout, at_from_rr, at_from_ss, at_from_duty, at_from_duty_direct, at_from_dwell_direct)
        
        # Use Shared Integration Logic
        at_actual_s = Logic.get_best_integration_time(at_needed, overhead_s, valid_times_s, p_s, inputs['icp_technology'])
        
        # --- REFINEMENT: Calculate Synchronized Rep Rate ---
        avoid_gaps = inputs.get('avoid_gaps', False)
        rr_theoretical = n_val / at_actual_s
        rr_prec = inputs.get('rr_prec_hz', 1.0)
        if rr_prec <= 0: rr_prec = 1.0
        allowed_rr = inputs.get('allowed_rr', None)
        strategy = inputs.get('sync_strategy', 'Adaptive Integer Sync (Auto)')
        
        # Apply Lockwood Eqn (3) threshold for Oversampling and Combined modes via dynamic sweep
        if "Oversampling" in strategy or "Combined" in strategy:
            target_rsd = inputs.get('target_rsd', 5.0)
            washout_ms = inputs.get('washout_ms', 30.0)
            
            # Speed and hardware limits
            max_speed = inputs.get('max_speed_um_s', 0.0)
            dosage = inputs.get('dosage', 10.0)
            spot_size = inputs.get('spot_size_um', 10.0)
            max_allowed_rr = inputs.get('max_rr_hz', 10000.0)
            if max_speed > 0 and dosage > 0 and spot_size > 0:
                rr_speed_limit = (max_speed * dosage) / spot_size
                max_allowed_rr = min(max_allowed_rr, rr_speed_limit)
                
            # Determine active channel dwell times (in ms)
            active_dwells = inputs.get('active_dwells_ms', [])
            if not active_dwells:
                if is_mc:
                    active_dwells = [(at_actual_s - overhead_s) * 1000.0]
                else:
                    active_dwells = [inputs.get('min_dwell_ms', 2.0)]
            
            # Dynamic sweep over candidate repetition rates f to find the minimum f satisfying target RSD.
            # Require at least 2.0 pulses per washout (f * W_s >= 2.0) to ensure a stable steady state (overlapping pulses).
            f_min = max_allowed_rr
            min_oversampling_f = 2.0 / (washout_ms / 1000.0) if washout_ms > 0 else 1.0
            sweep_start = max(1.0, min_oversampling_f)
            sweep_end = max_allowed_rr
            step = rr_prec
            if step <= 0:
                step = 1.0
                
            curr_f = sweep_start
            found = False
            while curr_f <= sweep_end + 1e-9:
                all_satisfy = True
                for dt_ms in active_dwells:
                    if dt_ms <= 0:
                        continue
                    rsd = Logic.calculate_lockwood_rsd(curr_f, washout_ms, dt_ms)
                    if rsd > target_rsd:
                        all_satisfy = False
                        break
                if all_satisfy:
                    f_min = curr_f
                    found = True
                    break
                curr_f += step
                
            if not found:
                f_min = max_allowed_rr
                
            rr_theoretical = max(rr_theoretical, f_min)
            rr_theoretical = min(rr_theoretical, max_allowed_rr)

        if "Oversampling" in strategy:
            # Oversampling: Bypasses pulse‑matching cycle snapping.
            # Respect allowed hardware integration times for Multi-Collector / TOF systems, otherwise snap to precision.
            if is_mc or inputs['icp_technology'] == "TOF":
                at_actual_s = Logic.get_best_integration_time(at_needed, overhead_s, valid_times_s, p_s, inputs['icp_technology'])
            else:
                at_actual_s = round(at_needed / p_s) * p_s
                prec_decimals = max(0, int(round(-math.log10(p_s)))) if p_s > 0 else 3
                at_actual_s = round(at_actual_s, prec_decimals)
            # In oversampling mode, the repetition rate is driven solely by the minimum
            # steady‑state frequency calculated from washout and target RSD (f_min).
            # Ignore n_val (pulses per acquisition) so the result is independent of the UI value.
            rr_theoretical = f_min
            if allowed_rr:
                valid = [r for r in allowed_rr if r >= rr_theoretical - 1e-9]
                rr_actual = min(valid) if valid else max(allowed_rr)
            else:
                rr_actual = math.ceil(rr_theoretical / rr_prec - 1e-9) * rr_prec
                
            if rr_actual <= 1e-9: rr_actual = rr_prec
            optimised_n = rr_actual * at_actual_s
        elif is_mc:
            # MC: Time-Driven
            if "Combined" in strategy:
                # Find the smallest F >= rr_theoretical (which is f_min) such that F is a multiple of rr_prec
                # and F * at_actual_s is an integer (within 1e-9 tolerance)
                start_f = math.ceil(rr_theoretical / rr_prec - 1e-9) * rr_prec
                curr_f = start_f
                found_sync = False
                max_allowed_rr = inputs.get('max_rr_hz', 10000.0)
                while curr_f <= max_allowed_rr + 1e-9:
                    pulses = curr_f * at_actual_s
                    if abs(pulses - round(pulses)) < 1e-9:
                        rr_actual = curr_f
                        found_sync = True
                        break
                    curr_f += rr_prec
                if not found_sync:
                    rr_actual = start_f
            elif avoid_gaps: 
                rr_actual = math.ceil(rr_theoretical / rr_prec - 1e-9) * rr_prec
            else: 
                rr_actual = math.floor(rr_theoretical / rr_prec + 1e-9) * rr_prec
            
            if "Oversampling" in strategy or "Combined" in strategy:
                optimised_n = rr_actual * at_actual_s
            else:
                optimised_n = n_val
        else:
            # Standard: Pulse-Driven. Re-calculate AT to match an integer RR if not avoid_gaps.
            if allowed_rr:
                # Calculate natural floor limits first
                rr_floor = max([r for r in allowed_rr if r <= rr_theoretical + 1e-9]) if any(r <= rr_theoretical + 1e-9 for r in allowed_rr) else min(allowed_rr)
                
                # Only enforce gap avoidance if the natural floor produces a negative error
                actually_avoid = avoid_gaps and (rr_floor * at_actual_s) < n_val
                
                valid = [r for r in allowed_rr if (r >= rr_theoretical - 1e-9 if actually_avoid else r <= rr_theoretical + 1e-9)]
                rr_actual = (min(valid) if actually_avoid else max(valid)) if valid else (max(allowed_rr) if actually_avoid else min(allowed_rr))
            else:
                rr_floor = math.floor(rr_theoretical / rr_prec + 1e-9) * rr_prec
                actually_avoid = avoid_gaps and (rr_floor * at_actual_s) < n_val
                rr_actual = (math.ceil(rr_theoretical / rr_prec - 1e-9) if actually_avoid else math.floor(rr_theoretical / rr_prec + 1e-9)) * rr_prec
            
            if rr_actual <= 1e-9: rr_actual = rr_prec
            
            if avoid_gaps:
                optimised_n = rr_actual * at_actual_s
            else:
                at_actual_s = round((n_val / rr_actual) / p_s) * p_s
                prec_decimals = max(0, int(round(-math.log10(p_s)))) if p_s > 0 else 3
                at_actual_s = round(at_actual_s, prec_decimals)
                optimised_n = n_val

        # --- PERFECT EVEN SPLIT: Adjust AT to avoid drift in Even mode ---
        # Final adjustment happens after all pulse-driven rounding to ensure it sticks.
        num_even = inputs.get('num_even', 0)
        fixed_costs_s = inputs.get('fixed_costs_s', 0.0)
        
        if num_even > 0 and not is_mc:
            # The remaining POOL budget (after fixed costs) must be perfectly divisible by (num_even * p_s)
            quantum = num_even * p_s
            available_budget_s = at_actual_s - overhead_s
            pool_budget_s = available_budget_s - fixed_costs_s
            
            if pool_budget_s < (num_even * p_s):
                at_actual_s = fixed_costs_s + (num_even * p_s) + overhead_s
            else:
                steps = math.ceil(pool_budget_s / quantum - 1e-9)
                at_actual_s = fixed_costs_s + (steps * quantum) + overhead_s
            
            # Recalculate RR to maintain synchronization with the nudged AT
            rr_theoretical_even = optimised_n / at_actual_s if at_actual_s > 0 else 0
            rr_rounded = round(rr_theoretical_even / rr_prec) * rr_prec
            
            if avoid_gaps:
                if (rr_rounded * at_actual_s) >= n_val:
                    rr_actual = rr_rounded
                else:    
                    rr_actual = math.ceil(rr_theoretical_even / rr_prec - 1e-9) * rr_prec
            else:
                rr_actual = rr_rounded

        # --- AUTO-SYNC SWEEP LOGIC ---
        strategy = inputs.get('sync_strategy', 'Adaptive Integer Sync (Auto)')
        if "Auto" in strategy and not is_mc:
            F_raw = rr_actual
            AT_raw_s = at_actual_s
            
            # Hardware resolution constants (always needed)
            at_hw_min_ms = inputs['min_dwell_ms'] + (overhead_s * 1000.0)
            at_res_ms    = inputs.get('min_dwell_ms', 0.1)
            washout_ms   = inputs.get('washout_ms', 30.0)
            allowed_rr   = inputs.get('allowed_rr', None)

            rr_prec_int  = max(1, int(round(rr_prec)))

            # Fractional pulse count before any sync adjustment
            pulses_raw = F_raw * AT_raw_s

            best_pair = None
            min_diff_score = float('inf')

            # User's original target pulse count (from the UI spinner)
            n_val_sweep  = int(round(inputs.get('pulses_per_pixel', round(pulses_raw))))
            prefer_exact = inputs.get('prefer_exact_pulses', False)

            f_min_combined = f_min if 'f_min' in locals() else 10.0

            def _try_candidate(F, n):
                nonlocal best_pair, min_diff_score
                if n < 1:
                    return
                AT_exact_ms = n * 1000.0 / F
                # Snap to hardware resolution
                AT_hw_ms = round(AT_exact_ms / at_res_ms) * at_res_ms
                # Must satisfy both hardware minimum and washout time
                if AT_hw_ms < max(at_hw_min_ms, washout_ms):
                    return
                # Reject if hardware rounding breaks integer sync
                pulses_check = F * AT_hw_ms / 1000.0
                if abs(pulses_check - round(pulses_check)) > 1e-6:
                    return
                AT_s    = AT_hw_ms / 1000.0
                if "Combined" in strategy:
                    # In Combined mode, we want the lowest frequency >= f_min_combined that achieves sync
                    # and keeps the acquisition time AT close to AT_raw_s.
                    # We avoid using F_raw (which is based on the hidden default target pulse count of 10).
                    score = (F - f_min_combined) / max(f_min_combined, 1.0) + abs(AT_s - AT_raw_s) / max(AT_raw_s, 1e-9)
                else:
                    diff_F  = abs(F - F_raw) / max(F_raw, 1.0)
                    diff_AT = abs(AT_s - AT_raw_s) / max(AT_raw_s, 1e-9)
                    pulse_penalty = (0.0 if (not prefer_exact or n == n_val_sweep)
                                     else 1.0 + abs(n - n_val_sweep))
                    score = diff_F + diff_AT + pulse_penalty
                if score < min_diff_score:
                    min_diff_score = score
                    best_pair = (F, AT_s)

            if "Combined" in strategy:
                # ── COMBINED (AUTO) ──────────────────────────────────────────
                # 1. Use the dynamically calculated oversampling f_min
                f_min_combined = f_min
                at_min_ms = max(at_hw_min_ms, washout_ms)

                # 2. Candidate rep‑rates at or above the oversampling minimum
                if allowed_rr:
                    rr_candidates = sorted(r for r in allowed_rr
                                           if r >= f_min_combined)
                else:
                    rr_start = max(1, int(math.ceil(f_min_combined / rr_prec_int)) * rr_prec_int)
                    rr_candidates = list(range(rr_start, int(rr_max) + 1, rr_prec_int))
                if not rr_candidates:
                    rr_candidates = [max(1, int(math.ceil(f_min_combined)))]

                # 3. For each rep‑rate, try n_min and n_min+1 above the floor
                for F in rr_candidates:
                    n_min = int(math.ceil(F * at_min_ms / 1000.0))
                    _try_candidate(F, n_min)
                    _try_candidate(F, n_min + 1)

            else:
                # Search the entire operational range of the laser (1 to rr_max)
                if allowed_rr:
                    rr_candidates = sorted(allowed_rr)
                else:
                    rr_candidates = list(range(rr_prec_int, int(rr_max) + 1, rr_prec_int))
                if not rr_candidates:
                    rr_candidates = [1]

                for F in rr_candidates:
                    # Ideal (fractional) pulse count for this rep-rate at the original AT
                    n_ideal = F * AT_raw_s
                    n_lo = max(1, int(math.floor(n_ideal)))
                    n_hi = n_lo + 1
                    _try_candidate(F, n_lo)
                    _try_candidate(F, n_hi)
                    # When prefer_exact is on, also explicitly try the user's
                    # target pulse count
                    if prefer_exact and n_val_sweep not in (n_lo, n_hi):
                        _try_candidate(F, n_val_sweep)

            if best_pair is not None:
                rr_actual = best_pair[0]
                at_actual_s = best_pair[1]
                optimised_n = rr_actual * at_actual_s

        # Derived parameters (CALCULATED LAST)
        actual_pulses = rr_actual * at_actual_s
        speed = (bs * rr_actual) / n_dose if n_dose > 0 else 0
        overlap_um = bs - (speed / rr_actual) if rr_actual > 0 else 0
        overlap_pct = (overlap_um / bs) * 100 if bs > 0 else 0
        
        # Consistent dwell budget calculation (Net measurement time)
        # For MC/TOF: This is the individual integration time.
        # For Quad/Sector: This is the total time for all dwells.
        dwell_budget_s = max(0.0, at_actual_s - overhead_s)

        # Sync error and warning (moved from calculate_laser_sync)
        n_target = inputs['pulses_per_pixel']
        sync_error = actual_pulses - n_target
        err_pct = (sync_error / n_target) * 100 if n_target != 0 else 0.0
        
        warning = f"Optimised to match {'Integration Time' if inputs['icp_technology'] == 'Multi-Collector' else 'Dwell Time'}: {actual_pulses:.4f} pulses (Error: {err_pct:+.4f}%)."

        return {
            "at_actual_s": at_actual_s,
            "at_needed_s": at_needed,
            "rr_actual": rr_actual,
            "speed": speed,
            "actual_pulses": actual_pulses,
            "overlap_pct": overlap_pct,
            "overlap_um": overlap_um,
            "harmonic": harmonic,
            "valid_times_s": valid_times_s,
            "notes": notes,
            "optimised_n": optimised_n,
            "n_dose": n_dose,
            "warning": warning,
            "err_pct": err_pct,
            "dwell_budget_s": dwell_budget_s
        }


    @staticmethod
    def calculate_sigma_for_spot(spot_size, config, isotope_data):
        target_sigma = config['lower_sigma_limit']
        precision_ms = config['precision_ms']
        precision_s = precision_ms / 1000.0
        snr_threshold = config.get('snr_threshold', 0.0)
        
        # 1. Sync machine parameters for this spot size
        c_inputs = config.copy()
        c_inputs['spot_size_um'] = spot_size
        

            
        res = Logic.calculate_constrained_at(c_inputs)
        
        if 'override_at_s' in config and 'override_rr' in config:
            res['rr_actual'] = config['override_rr']
            res['at_actual_s'] = config['override_at_s']
            res['dwell_budget_s'] = max(0.0, config['override_at_s'] - (config['overhead_ms'] / 1000.0))
            res['optimised_n'] = res['rr_actual'] * res['at_actual_s']
            res['n_dose'] = config['dosage'] if not config.get('sync_dosage', True) else res['optimised_n']
            res['speed'] = (spot_size * res['rr_actual']) / res['n_dose'] if res['n_dose'] > 0 else 0.0
            res['overlap_um'] = spot_size - (res['speed'] / res['rr_actual']) if res['rr_actual'] > 0 else spot_size
            res['actual_pulses'] = res['rr_actual'] * res['at_actual_s']
            n_target = config['pulses_per_pixel']
            sync_error = res['actual_pulses'] - n_target
            res['err_pct'] = (sync_error / n_target) * 100 if n_target != 0 else 0.0
            res['warning'] = f"Optimised to match {'Integration Time' if config['icp_technology'] == 'Multi-Collector' else 'Dwell Time'}: {res['actual_pulses']:.4f} pulses (Error: {res['err_pct']:+.4f}%)."
            
        at_constrained = res['at_actual_s']
        rr_actual = res['rr_actual']
        dwell_budget_s = res['dwell_budget_s']
        min_dwell_s = c_inputs['min_dwell_ms'] / 1000.0
        
        decimals = 0 if precision_ms >= 1.0 else int(math.ceil(-math.log10(precision_ms)))
        scaling_factor = (spot_size / config['ref_spot_size_um'])**2
        
        # 2. Determine RR scaling using centralized sync logic
        rr_scaling = 1.0
        if config.get('scale_signal', False):
            initial_rr = config.get('initial_rr', 0)
            if initial_rr > 0:
                rr_scaling = rr_actual / initial_rr

        processed_iso = []
        for iso in isotope_data:
            # Use centralized SNR logic
            snr_res = Logic.calculate_robust_snr(iso['sig_cps'], iso['blk_cps'], iso['baseline_dt_s'], iso['stdev_blank_counts'])
            
            # Apply spot/rep-rate scaling to the SNR
            # Scaling: Net counts scale with scaling_factor * rr_scaling
            # SNR scales linearly with net counts for a fixed background
            scaled_base_snr = snr_res['base_snr'] * (scaling_factor * rr_scaling)
            
            initial_snr_display = iso.get('initial_snr', 0.0)
            is_optimisable = (scaled_base_snr > 0 and iso['status'] not in ["Exclude", "Set to Min", "Custom"])

            processed_iso.append({
                **iso, 
                'base_snr': scaled_base_snr,  # This is the scaled Detection Ratio at bl_dt
                'initial_snr_display': initial_snr_display,
                'is_optimisable': is_optimisable, 
                'final_dt': 0.0, 
                'snapped_dt': 0.0,
                'is_zero_bg': snr_res['is_zero_bg']
            })

        # Selection of integration time/dwell share
        # Master dwell budget from synchronized logic ensures MC/TOF integration times match the UI
        dwell_budget_s = res['dwell_budget_s']

        if config['icp_technology'] == "Multi-Collector" or config['icp_technology'] == "TOF":
            for iso in processed_iso:
                iso['snapped_dt'] = 0 if iso['status'] == "Exclude" else dwell_budget_s
        else:
            # Quad/Optimization Mode
            for _ in range(len(processed_iso) + 5):
                current_budget = dwell_budget_s
                inv_snr_sum = 0
                for iso in processed_iso:
                    if iso['status'] == "Exclude": iso['final_dt'] = 0.0
                    elif iso['status'] == "Set to Min": iso['final_dt'] = min_dwell_s; current_budget -= min_dwell_s
                    elif iso['status'] == "Custom": val = iso.get('custom_time_s', min_dwell_s); iso['final_dt'] = val; current_budget -= val
                    elif not iso['is_optimisable']: iso['final_dt'] = min_dwell_s; current_budget -= min_dwell_s
                    else: 
                        iso['final_dt'] = 0.0
                        inv_snr_sum += (math.sqrt(iso['baseline_dt_s']) / iso['base_snr']) if iso['base_snr'] > 0 else 0
                
                valid_indices = [i for i, x in enumerate(processed_iso) if x['is_optimisable'] and x['final_dt'] == 0]
                num_valid = len(valid_indices)
                reserved_time = num_valid * min_dwell_s
                extra_budget = current_budget - reserved_time
                
                if extra_budget < 0:
                    for idx in valid_indices: processed_iso[idx]['final_dt'] = min_dwell_s
                elif num_valid > 0:
                    # Check if we are using "Even" distribution
                    is_even = any(processed_iso[idx]['status'] == "Even" for idx in valid_indices)
                    
                    for idx in valid_indices:
                        if is_even:
                            share_ratio = 1.0 / num_valid
                        else:
                            # Auto (SNR-weighted)
                            share_ratio = (math.sqrt(processed_iso[idx]['baseline_dt_s']) / processed_iso[idx]['base_snr']) / inv_snr_sum if inv_snr_sum > 0 else (1.0 / num_valid)
                            
                        processed_iso[idx]['final_dt'] = min_dwell_s + (extra_budget * share_ratio)

                failures = []
                for idx in valid_indices:
                    res_snr = processed_iso[idx]['base_snr'] * math.sqrt(max(0.0, processed_iso[idx]['final_dt']) / processed_iso[idx]['baseline_dt_s']) if processed_iso[idx]['baseline_dt_s'] > 0 else 0
                    if res_snr < snr_threshold: failures.append((idx, res_snr))
                
                if not failures: break
                else:
                    failures.sort(key=lambda x: x[1])
                    processed_iso[failures[0][0]]['is_optimisable'] = False

            total_snapped = 0.0
            for iso in processed_iso:
                if iso['status'] == "Exclude": snapped = 0.0
                else:
                    snapped = max(min_dwell_s, math.floor(iso['final_dt'] / precision_s + 1e-9) * precision_s)
                    if not iso['is_optimisable'] and iso['status'] in ["Auto", "Even"]: iso['constraint'] = "Min SNR"
                    elif abs(snapped - min_dwell_s) < 1e-9 and iso['status'] in ["Auto", "Even"]: iso['constraint'] = "Min ICP"
                    elif iso['status'] == "Set to Min": iso['constraint'] = "Min ICP"
                    else: iso['constraint'] = ""
                iso['snapped_dt'] = snapped
                total_snapped += snapped

            # Redistribute 'drift' (leftover budget after rounding down)
            drift = dwell_budget_s - total_snapped
            if drift >= (precision_s * 0.9):
                # Only give drift to isotopes part of the 'Auto' or 'Even' optimization
                candidates = [i for i, x in enumerate(processed_iso) if x['is_optimisable'] and x['status'] not in ["Exclude", "Custom", "Set to Min"]]
                if candidates:
                    # Sort by current Sigma Separation (give to the worst performers first)
                    candidates.sort(key=lambda i: processed_iso[i]['base_snr'] * math.sqrt(max(0.0, processed_iso[i]['snapped_dt']) / processed_iso[i]['baseline_dt_s']) if processed_iso[i]['baseline_dt_s'] > 0 else 0)
                    num_steps = int(round(drift / precision_s))
                    for idx in range(num_steps):
                        processed_iso[candidates[idx % len(candidates)]]['snapped_dt'] += precision_s

        separations = []; output_rows = []
        for iso in processed_iso:
            sep = iso['base_snr'] * math.sqrt(max(0.0, iso['snapped_dt']) / iso['baseline_dt_s']) if iso['baseline_dt_s'] > 0 else 0
            if iso['status'] not in ["Exclude", "Set to Min", "Custom"]: separations.append(sep)
            val_ms = round(iso['snapped_dt'] * 1000.0, decimals)
            output_rows.append({"Isotope": iso['name'], "Final Dwell (ms)": val_ms if decimals > 0 else int(val_ms),
                                "Initial SNR": round(iso['initial_snr_display'], 2), "Sigma Sep": round(sep, 2),
                                "Status": iso['status'], "Constraint": iso.get('constraint', ""), "IsZeroBG": iso.get('is_zero_bg', False),
                                "Signal CPS": iso.get('sig_cps', 0.0)})

        return (min(separations) if separations else 999.0), output_rows, res


    @staticmethod
    def calculate_minimum_required_spot_size(config, isotope_data):
        target_sigma = config['lower_sigma_limit']
        low, high = 1, 300
        best_spot, best_res = high, None
        while low <= high:
            mid = (low + high) // 2 
            sigma, rows, _ = Logic.calculate_sigma_for_spot(mid, config, isotope_data)
            if sigma >= target_sigma: best_spot = mid; best_res = rows; high = mid - 1
            else: low = mid + 1
        if best_res is None:
            _, best_res, _ = Logic.calculate_sigma_for_spot(300, config, isotope_data)
        return best_spot, best_res

    @staticmethod
    def calculate_signal_statistics(df, isotope_cols, bg_range, sig_range, edited_dwells_df, estimated_dwell_ms, is_raw_counts, show_counts_check):
        df_blk = df[(df["Time"] >= bg_range[0]) & (df["Time"] <= bg_range[1])]
        df_sig = df[(df["Time"] >= sig_range[0]) & (df["Time"] <= sig_range[1])]
        temp_data = []
        for col in isotope_cols:
            dwell_s = edited_dwells_df[col].iloc[0] / 1000.0 if col in edited_dwells_df.columns else estimated_dwell_ms / 1000.0
            
            stats = Logic.get_isotope_stats(col, df_sig, df_blk, dwell_s, is_raw_counts)
            snr_res = Logic.calculate_robust_snr(stats['sig_cps'], stats['blk_cps'], dwell_s, stats['std_counts'])
            
            factor = 1.0 if (is_raw_counts == show_counts_check) else (dwell_s if show_counts_check else 1/dwell_s)
            temp_data.append({
                "name": col, "current_dwell_ms": dwell_s*1000, 
                "disp_sig": stats['m_sig']*factor, "disp_std_sig": df_sig[col].std()*factor, 
                "disp_blk": stats['m_blk']*factor, "disp_std_blk": stats['s_blk']*factor, 
                "initial_snr": snr_res['base_snr']
            })
        return temp_data

    @staticmethod
    def prepare_rows_for_optimization(edited_stats_df, iso_status_map, iso_custom_dt_map, df, bg_range, sig_range, is_raw_counts):
        df_sig = df[(df["Time"] >= sig_range[0]) & (df["Time"] <= sig_range[1])]
        df_blk = df[(df["Time"] >= bg_range[0]) & (df["Time"] <= bg_range[1])]
        rows = []
        for row in edited_stats_df:
            name = row['name']
            row['status'] = iso_status_map.get(name, "Auto")
            row['custom_time_s'] = iso_custom_dt_map.get(name, 0.0) / 1000.0
            dwell_s = row['current_dwell_ms'] / 1000.0
            
            stats = Logic.get_isotope_stats(name, df_sig, df_blk, dwell_s, is_raw_counts)
            row['sig_cps'] = stats['sig_cps']
            row['blk_cps'] = stats['blk_cps']
            row['stdev_blank_counts'] = stats['std_counts']
            row['baseline_dt_s'] = dwell_s
            rows.append(row)
        return rows, df_sig.copy(), df_blk.copy()

    @staticmethod
    def auto_detect_regions(df):
        try:
            numeric_cols = [c for c in df.columns if c != 'Time' and np.issubdtype(df[c].dtype, np.number)]
            if not numeric_cols: return (0.0, 1.0), (2.0, 3.0)
            
            y = df[numeric_cols].sum(axis=1).values
            time_vals = df['Time'].values
            
            # 1. Very sensitive baseline estimation
            init_window = max(10, int(len(y) * 0.05))
            bg_est = y[:init_window]
            bg_mean = np.mean(bg_est)
            bg_std = np.std(bg_est) if len(bg_est) > 1 else 0
            
            y_max = np.max(y)
            # Threshold: 10.0 sigma above noise, or 10% of max range (User Request: Stricter filtering)
            thresh = max(bg_mean + 10.0 * bg_std, bg_mean + (y_max - bg_mean) * 0.1)
            
            # 2. Narrower rolling max to handle pulsed data without excessive smear
            dt = time_vals[1] - time_vals[0] if len(time_vals) > 1 else 0.1
            window_size = max(3, int(0.5 / dt)) # 0.5s window
            y_series = pd.Series(y)
            y_smooth = y_series.rolling(window=window_size, center=True).max().fillna(y_series).values
            
            sig_indices = np.where(y_smooth > thresh)[0]
            
            if len(sig_indices) > 5:
                # 3. Island-based noise filtering (User Request: Ignore pre-pulse noise)
                # Group into contiguous clusters (islands)
                islands = []
                if len(sig_indices) > 0:
                    curr = [sig_indices[0]]
                    for i in range(1, len(sig_indices)):
                        if sig_indices[i] == sig_indices[i-1] + 1:
                            curr.append(sig_indices[i])
                        else:
                            islands.append(curr)
                            curr = [sig_indices[i]]
                    islands.append(curr)
                
                # Filter noise spikes: Real signal has "mass" (duration or intensity)
                y_max_excess = y_max - bg_mean
                valid_islands = []
                for isl in islands:
                    duration = time_vals[isl[-1]] - time_vals[isl[0]]
                    peak_ex = np.max(y[isl]) - bg_mean
                    # Laser pulses usually have duration > 0.5s OR > 10% of global intensity
                    if duration > 0.5 or peak_ex > 0.1 * y_max_excess:
                        valid_islands.append(isl)
                
                # If we have valid islands, use them. Otherwise fallback to raw clusters to be safe.
                filtered_indices = [idx for isl in (valid_islands if valid_islands else islands) for idx in isl]
                
                s_idx, e_idx = filtered_indices[0], filtered_indices[-1]
                
                # Signal Region: Symmetric 5% crop
                sig_len = e_idx - s_idx
                sig_crop = int(sig_len * 0.05)
                rec_sig = (time_vals[s_idx + sig_crop], time_vals[e_idx - sig_crop])
                
                # Background Region: Asymmetric crop to ensure clearing
                # We pull back significantly (15%) from the detected signal start
                gap_len = s_idx
                bg_start_idx = int(gap_len * 0.05)
                bg_end_idx = s_idx - int(gap_len * 0.15)
                
                if bg_end_idx > bg_start_idx + 2:
                    rec_bg = (time_vals[bg_start_idx], time_vals[bg_end_idx])
                else:
                    # Fallback if space is very tight
                    rec_bg = (time_vals[0], time_vals[max(1, s_idx // 2)])
                
                return rec_bg, rec_sig
        except Exception as e:
            # IoLog.warning(f"iolite Optimiser: Auto-detect error: {e}")
            pass
            
        t = df['Time'].values
        return (t[0], t[len(t)//5]), (t[len(t)//3], t[-1])




class CopyableTableWidget(QTableWidget):
    def _get(self, obj, attr):
        return getattr(obj, attr)() if callable(getattr(obj, attr)) else getattr(obj, attr)

    def keyPressEvent(self, event):
        # IoLog.information(f"Key Press: {event.key()} Mod: {event.modifiers()}")
        if (event.key() == Qt.Key_C and (event.modifiers() & Qt.ControlModifier)):
            self.copy_selection()
            event.accept()
            return
        super().keyPressEvent(event)

    def contextMenuEvent(self, event):
        from iolite.QtGui import QMenu, QAction
        menu = QMenu(self)
        copy_action = QAction("Copy Selection", self)
        copy_action.triggered.connect(self.copy_selection)
        menu.addAction(copy_action)
        menu.exec_(event.globalPos())

    def copy_selection(self):
        try:
            selection = self.selectedRanges()
            if selection:
                # Fix: Rename variable 'range' to 'sel_range' to avoid shadowing built-in range()
                rows = sorted(list(set(r for sel_range in selection for r in range(sel_range.topRow(), sel_range.bottomRow() + 1))))
                cols = sorted(list(set(c for sel_range in selection for c in range(sel_range.leftColumn(), sel_range.rightColumn() + 1))))
                
                if not rows or not cols:
                    return

                text_grid = []
                for r in rows:
                    row_data = []
                    for c in cols:
                        is_sel = False
                        for sel_range in selection:
                            if sel_range.topRow() <= r <= sel_range.bottomRow() and sel_range.leftColumn() <= c <= sel_range.rightColumn():
                                is_sel = True; break
                        
                        if is_sel:
                            item = self.item(r, c)
                            row_data.append(self._get(item, 'text') if item else "")
                        else:
                            row_data.append("")
                    text_grid.append("\t".join(row_data))
                
                final_text = "\n".join(text_grid)
                QApplication.clipboard().setText(final_text)
                IoLog.information(f"Copied {len(rows)} rows to clipboard.")
        except Exception as e:
            IoLog.error(f"Copy Error: {e}")


class SettingsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Hardware Configuration")
        # Fixed size: Compacted width and stable height (Reduced to 500)
        self.main_layout = QVBoxLayout(self)
        self.setFixedSize(480, 500)
        self.main_layout.setContentsMargins(10, 10, 10, 10)
        self.main_layout.setSpacing(5)

        # Added stretch to keep group boxes at the top and prevent them from
        # being stretched vertically to fill the window height.
        self.main_layout.addStretch(1)

        # We will add groupboxes to this layout from initUI
        
        # Add a Close button and version label at the bottom
        btns = QHBoxLayout()
        lbl_version = QLabel(f"Version: {VERSION}")
        lbl_version.setStyleSheet("color: gray; font-size: 8pt;")
        btns.addWidget(lbl_version)
        btns.addStretch()
        self.btn_close = QPushButton("Close")
        self.btn_close.setFixedWidth(120)
        self.btn_close.setFixedHeight(30)
        self.btn_close.clicked.connect(lambda: self.done(1))
        btns.addWidget(self.btn_close)
        self.main_layout.addLayout(btns)


class PresetManagerDialog(QDialog):
    """Dialog to list, view, and delete saved channel presets."""
    def __init__(self, presets, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Manage Presets")
        self.resize(450, 400)
        self.presets = presets.copy()
        
        layout = QVBoxLayout(self)
        
        lbl = QLabel("Saved Presets:")
        lbl.setStyleSheet("font-weight: bold;")
        layout.addWidget(lbl)
        
        self.list_widget = QListWidget()
        for name in self.presets.keys():
            self.list_widget.addItem(name)
        self.list_widget.itemSelectionChanged.connect(self._update_preview)
        layout.addWidget(self.list_widget)
        
        # Preview Area
        layout.addWidget(QLabel("Channels in selected preset:"))
        self.preview = QTextEdit()
        self.preview.setReadOnly(True)
        layout.addWidget(self.preview)
        
        # Buttons
        h_btns = QHBoxLayout()
        btn_del = QPushButton("Delete Selected")
        btn_del.clicked.connect(self._delete_preset)
        h_btns.addWidget(btn_del)
        
        h_btns.addStretch()
        
        btn_close = QPushButton("Done")
        btn_close.clicked.connect(lambda: self.done(1))
        h_btns.addWidget(btn_close)
        
        layout.addLayout(h_btns)

    def _update_preview(self):
        item = self.list_widget.currentItem()
        if item:
            name = item.text()
            if callable(name): name = name()
            channels = self.presets.get(name, [])
            self.preview.setText(", ".join(channels))
        else:
            self.preview.clear()

    def _delete_preset(self):
        item = self.list_widget.currentItem()
        if not item:
            return
        
        name = item.text()
        if callable(name): name = name()
        
        res = QMessageBox.question(self, "Delete Preset", f"Are you sure you want to delete preset '{name}'?",
                                 QMessageBox.Yes | QMessageBox.No)
        
        if res == QMessageBox.Yes:
            if name in self.presets:
                del self.presets[name]
                self.list_widget.takeItem(self.list_widget.row(item))
                self.preview.clear()


class DataConfigDialog(QDialog):
    """
    Unified dialog for data configuration:
    - For TOF/Vitesse: Channel selection checkboxes + global dwell time
    - For Quad/Sector: Just dwell time configuration (per-channel or global)
    """
    def __init__(self, channels, detected_at=None, is_tof=False, parent=None):
        super().__init__(parent)
        self.is_tof = is_tof
        # Filter out 'TotalBeam'
        channels = [c for c in channels if str(c).strip().lower() != "totalbeam"]
        self.channels = channels
        self.result_channels = channels.copy()  # Default: all channels selected
        self.result_dwells = {}
        
        if is_tof:
            self.setWindowTitle("Import TOF Data")
            self.resize(400, 550)
        else:
            self.setWindowTitle("Configure Dwell Times")
            self.resize(400, 520)
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(15, 15, 15, 15)
        layout.setSpacing(5)
        
        # 0. Detected AT Header
        if detected_at:
            lbl_at = QLabel(f"Detected Acquisition Time: <b>{detected_at:.3f} ms</b>")
            lbl_at.setIndent(0)
            lbl_at.setStyleSheet("margin: 0px; padding: 0px; margin-bottom: 5px;")
            layout.addWidget(lbl_at)
        
        # --- TOF MODE: Channel Selection ---
        if is_tof:
            lbl = QLabel("Select channels to load:")
            lbl.setWordWrap(True)
            lbl.setIndent(0)
            lbl.setStyleSheet("font-weight: bold; margin: 0px; padding: 0px;")
            layout.addWidget(lbl)
            
            # Select All / None
            h_sel = QHBoxLayout()
            btn_all = QPushButton("Select All")
            btn_all.setAutoDefault(False)
            btn_all.setDefault(False)
            btn_none = QPushButton("Select None")
            btn_none.setAutoDefault(False)
            btn_none.setDefault(False)
            btn_all.clicked.connect(self.select_all)
            btn_none.clicked.connect(self.select_none)
            h_sel.addWidget(btn_all)
            h_sel.addWidget(btn_none)
            h_sel.addStretch()
            
            # Channel counter
            self.lbl_count = QLabel(f"0 / {len(channels)} selected")
            h_sel.addWidget(self.lbl_count)
            layout.addLayout(h_sel)
            
            # Scroll Area for Checkboxes
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            self.scroll_content = QWidget()
            self.scroll_layout = QVBoxLayout(self.scroll_content)
            self.scroll_layout.setContentsMargins(5,5,5,5)
            self.scroll_layout.setSpacing(2)
            
            self.checks = []
            for ch in channels:
                chk = QCheckBox(ch)
                chk.setChecked(False)  # Default off for TOF
                chk.toggled.connect(self._update_count)
                self.scroll_layout.addWidget(chk)
                self.checks.append(chk)
                
            self.scroll_layout.addStretch()
            scroll.setWidget(self.scroll_content)
            layout.addWidget(scroll)
            
            # Preset Selection
            layout.addWidget(QLabel("Load Preset:"))
            self.cmb_preset = QComboBox()
            self.cmb_preset.addItem("-- Select --")
            self._load_presets()
            self.cmb_preset.currentTextChanged.connect(self._apply_preset)
            layout.addWidget(self.cmb_preset)
            
            # Global Dwell for TOF
            h_dwell = QHBoxLayout()
            h_dwell.addWidget(QLabel("Global Dwell Time (ms):"))
            self.spin_dwell = QDoubleSpinBox()
            self.spin_dwell.setRange(0.001, 10000)
            self.spin_dwell.setValue(detected_at if detected_at else 10.0)
            self.spin_dwell.setDecimals(3)
            h_dwell.addWidget(self.spin_dwell)
            h_dwell.addStretch()
            layout.addLayout(h_dwell)
            
        # --- QUAD/SECTOR MODE: Dwell Configuration ---
        else:
            self.checks = None
            
            # Toggle: Global vs Per-Channel
            # User Request: Default to Global, Label "Set individual dwell times"
            self.chk_same = QCheckBox("Set individual dwell times")
            self.chk_same.setChecked(False) # Unchecked = Global (Default)
            self.chk_same.toggled.connect(self._toggle_mode)
            layout.addWidget(self.chk_same)
            
            self.stack = QStackedWidget()
            layout.addWidget(self.stack)
            
            # Page 0: Global Spinbox (Default)
            self.page_global = QWidget()
            l_global = QVBoxLayout(self.page_global)
            l_global.setContentsMargins(0,0,0,0)
            l_global.addWidget(QLabel("Global Dwell Time (ms):"))
            self.spin_dwell = QDoubleSpinBox()
            self.spin_dwell.setRange(0.001, 10000)
            self.spin_dwell.setValue(detected_at if detected_at else 10.0)
            self.spin_dwell.setDecimals(3)
            l_global.addWidget(self.spin_dwell)
            l_global.addStretch()
            self.stack.addWidget(self.page_global) # Index 0
            
            # Page 1: Individual Table
            self.page_table = QWidget()
            l_table = QVBoxLayout(self.page_table)
            l_table.setContentsMargins(0,0,0,0)
            
            self.table = QTableWidget()
            self.table.setColumnCount(2)
            self.table.setHorizontalHeaderLabels(["Isotope", "Dwell (ms)"])
            self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.table.setRowCount(len(channels))
            
            self.spin_map = {}
            for i, ch_name in enumerate(channels):
                self.table.setItem(i, 0, QTableWidgetItem(ch_name))
                sp = QDoubleSpinBox()
                sp.setRange(0.001, 10000)
                sp.setValue(detected_at if detected_at else 10.0)
                sp.setDecimals(3)
                self.table.setCellWidget(i, 1, sp)
                self.spin_map[ch_name] = sp
                
            l_table.addWidget(self.table)
            self.stack.addWidget(self.page_table)  # Index 1
        
        # --- Buttons ---
        btns = QHBoxLayout()
        
        # Preset buttons on left (only for TOF mode)
        if is_tof:
            btn_save_preset = QPushButton("Save Preset")
            btn_save_preset.setAutoDefault(False)
            btn_save_preset.clicked.connect(self._save_preset)
            btns.addWidget(btn_save_preset)
            
            btn_manage = QPushButton("Manage Presets")
            btn_manage.setAutoDefault(False)
            btn_manage.clicked.connect(self._manage_presets)
            btns.addWidget(btn_manage)
        
        btns.addStretch()
        self.btn_ok = QPushButton("OK")
        self.btn_ok.clicked.connect(self.accept)
        btns.addWidget(self.btn_ok)
        
        self.btn_cancel = QPushButton("Cancel")
        self.btn_cancel.clicked.connect(self.reject)
        btns.addWidget(self.btn_cancel)
        layout.addLayout(btns)
    
    def _toggle_mode(self, checked):
        if hasattr(self, 'stack'):
            self.stack.setCurrentIndex(1 if checked else 0)
    
    def _update_count(self):
        if hasattr(self, 'lbl_count') and self.checks:
            count = sum(1 for chk in self.checks if chk.isChecked())
            total = len(self.checks)
            self.lbl_count.setText(f"{count} / {total} selected")
    
    def select_all(self):
        if self.checks:
            for chk in self.checks: chk.setChecked(True)
            self._update_count()
            if hasattr(self, 'cmb_preset'):
                self.cmb_preset.blockSignals(True)
                self.cmb_preset.setCurrentIndex(0)
                self.cmb_preset.blockSignals(False)
        
    def select_none(self):
        if self.checks:
            for chk in self.checks: chk.setChecked(False)
            self._update_count()
            if hasattr(self, 'cmb_preset'):
                self.cmb_preset.blockSignals(True)
                self.cmb_preset.setCurrentIndex(0)
                self.cmb_preset.blockSignals(False)
    
    def _get_preset_file(self):
        base_dir = get_settings_dir()
        return os.path.join(base_dir, 'iolite Optimiser Channel Presets.json')
    
    def _load_presets(self):
        """Load presets from JSON file into the dropdown."""
        try:
            path = self._get_preset_file()
            if os.path.exists(path):
                with open(path, 'r') as f:
                    self._presets = json.load(f)
                for name in self._presets.keys():
                    self.cmb_preset.addItem(name)
            else:
                self._presets = {}
        except Exception as e:
            IoLog.warning(f"iolite Optimiser: Could not load presets: {e}")
            self._presets = {}
    
    def _save_presets_to_file(self):
        """Save all presets to JSON file."""
        import json
        try:
            path = self._get_preset_file()
            # Ensure we are saving clean Python types
            clean_presets = {}
            for k, v in self._presets.items():
                clean_k = str(k)
                clean_v = [str(x) for x in v]
                clean_presets[clean_k] = clean_v

            with open(path, 'w') as f:
                json.dump(clean_presets, f, indent=2)
            IoLog.information(f"iolite Optimiser: Presets saved to {path} ({len(clean_presets)} presets)")
        except Exception as e:
            IoLog.error(f"iolite Optimiser: Error saving presets to file: {e}")
            IoLog.error(traceback.format_exc())
    
    def _apply_preset(self, preset_name):
        """Apply a preset to the checkboxes."""
        if not hasattr(self, '_presets'):
            return
            
        if preset_name == "-- Select --":
            # Manual return to default: Uncheck all
            if self.checks:
                for chk in self.checks: chk.setChecked(False)
            self._update_count()
            return
        
        preset_channels = self._presets.get(preset_name, [])
        if not preset_channels:
            return
        
        # Uncheck all, then check matching
        for chk in self.checks:
            t = chk.text
            if callable(t): t = t()
            chk.setChecked(t in preset_channels)
        
        self._update_count()
    
    def _save_preset(self):
        """Prompt for name and save current selection as a preset."""
        try:
            res = QInputDialog.getText(self, "Save Preset", "Enter preset name:")
            
            # Handle variable return types from different Qt wrappers (PythonQt vs PyQt)
            if isinstance(res, (list, tuple)):
                name = str(res[0]).strip()
                ok = res[1] if len(res) > 1 else bool(name)
            else:
                name = str(res).strip()
                ok = bool(name)
                
            if ok and name:
                name = str(name).strip() # Explicit Python string
                selected = []
                for chk in self.checks:
                    if chk.isChecked():
                        t = chk.text
                        if callable(t): t = t()
                        selected.append(str(t)) # Explicit Python string
                
                if not selected:
                    QMessageBox.warning(self, "Empty Preset", "No channels selected to save.")
                    return
                
                self._presets[name] = selected
                self._save_presets_to_file()
                
                # Add to dropdown if not already there
                if self.cmb_preset.findText(name) == -1:
                    self.cmb_preset.addItem(name)
                
                # Explicitly set and apply
                self.cmb_preset.setCurrentText(name)
                self._apply_preset(name) # Force apply in case setCurrentText didn't trigger signal
                
                IoLog.information(f"iolite Optimiser: Saved preset '{name}' with {len(selected)} channels")
        except Exception as e:
            IoLog.error(f"iolite Optimiser: Failed to save preset: {e}")
            IoLog.error(traceback.format_exc())
    
    def _manage_presets(self):
        """Show dialog to manage (view/delete) presets."""
        dlg = PresetManagerDialog(self._presets, self)
        if dlg.exec_():
            self._presets = dlg.presets
            self._save_presets_to_file()
            # Refresh dropdown
            self.cmb_preset.clear()
            self.cmb_preset.addItem("-- Select --")
            for name in self._presets.keys():
                self.cmb_preset.addItem(name)
    
    def _get(self, obj, attr):
        if not obj or not hasattr(obj, attr): return None
        try:
            val = getattr(obj, attr)
            return val() if callable(val) else val
        except:
            return None
        
    def accept(self):
        try:
            # TOF Mode: Collect selected channels + global dwell
            if self.is_tof:
                selected = []
                for chk in self.checks:
                    if chk.isChecked():
                        t = chk.text
                        if callable(t): t = t()
                        selected.append(t)
                self.result_channels = selected
                
                dwell_val = self._get(self.spin_dwell, 'value')
                for ch in selected:
                    self.result_dwells[ch] = dwell_val
                    
                IoLog.information(f"iolite Optimiser: TOF config - {len(selected)} channels, global dwell {dwell_val}ms")
            
            # Quad Mode: Collect dwell config
            else:
                self.result_channels = self.channels  # All channels
                
                if hasattr(self, 'chk_same') and self._get(self.chk_same, 'isChecked'):
                    # Per-channel mode (Checked = Individual)
                    for ch, sp in self.spin_map.items():
                        v = self._get(sp, 'value')
                        if v: self.result_dwells[ch] = v
                else:
                    # Global mode (Unchecked = Global)
                    dwell_val = self._get(self.spin_dwell, 'value')
                    for ch in self.channels:
                        self.result_dwells[ch] = dwell_val
                        
                IoLog.information(f"iolite Optimiser: Dwell config for {len(self.result_dwells)} channels")
            
            self.done(QDialog.Accepted)
        except Exception as e:
            IoLog.error(f"iolite Optimiser: DataConfigDialog Error: {e}")
            self.done(QDialog.Rejected)
    
    def reject(self):
        self.done(QDialog.Rejected)


class ioliteOptimiser(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.bg_times = None
        self.sig_times = None
        self.opt_df = None
        self.spr_df = None
        self.is_external_data = False
        self._opt_palette_indices = None
        self._opt_axis_limits = None
        self._suppress_margin_update = False # Flag to prevent margin shudder during rebuilds
        self.dragged_edge = None # For interactive region adjustment
        self.last_opt_rows = 1
        self.last_spr_rows = 1
        self.t_start = 0.0
        self.channel_dwells = {} # Map channel_name -> dwell_ms
        self.settings_json_path = None
        self.persistent_settings = {}
        self.is_connected = False
        self.optimum_spotsize = None
        self.override_spotsize = None
        self.isotope_configs = {} # Session-specific settings for each isotope
        self.channel_metadata = {} # {"Name": {"Element": "Pb", "Mass": "208"}}
        self.show_meta = {'Element': False, 'Mass': False}
        self.col_map = {} # Mapping of header names to current indices
        self.active_washout_level = 0.01 # Default to 1% washout
        self._last_simulated_params = None # Cache of last simulation run parameters
        
        # Plot Persistence (Optimisation Tab Only)
        self._opt_palette_indices = None
        self._opt_axis_limits = None
        
        # Initialize System Theme Early
        self.system_is_dark = self.detect_theme()
        
        # Debounce timer for auto-refresh
        self.refresh_timer = QTimer()
        self.refresh_timer.setSingleShot(True)
        self.refresh_timer.setInterval(500)
        self.refresh_timer.timeout.connect(self._perform_auto_refresh)

        # Delayed click detection for legend (ensures double-click works with large datasets)
        self._legend_click_timer = QTimer()
        self._legend_click_timer.setSingleShot(True)
        self._legend_click_timer.setInterval(300)  # 300ms window for double-click
        self._legend_click_timer.timeout.connect(self._execute_legend_single_click)
        self._pending_legend_click = None  # Stores (target_obj, target_list, canvas, ax, l_map, l_chk_rescale)
        
        # Debounce timer for SPR peak detection
        self.spr_debounce_timer = QTimer()
        self.spr_debounce_timer.setSingleShot(True)
        self.spr_debounce_timer.timeout.connect(self._run_spr_analysis_forced)

        # Debounce timer for Pulse Train Simulation
        self.pulse_debounce_timer = QTimer()
        self.pulse_debounce_timer.setSingleShot(True)
        self.pulse_debounce_timer.setInterval(500)
        self.pulse_debounce_timer.timeout.connect(self.run_pulse_simulation)

        # Debounce timer for saving persistent settings
        self.save_timer = QTimer()
        self.save_timer.setSingleShot(True)
        self.save_timer.setInterval(1000)
        self.save_timer.timeout.connect(self.save_persistent_settings)

        # Initialize hardware state
        self.icp_tech = "Quadrupole"
        self.min_dwell = 0.1
        self.precision = 0.1
        self.max_rr = 1000
        self.allowed_rr = None
        self.allowed_dwells = None
        self.max_speed = 20000
        
        try:
            ol_path = get_settings_dir()
            self.settings_json_path = os.path.join(ol_path, "iolite Optimiser Settings.json")
            self.persistent_settings = self.load_persistent_settings()
        except Exception as e:
            self.persistent_settings = {}
            IoLog.warning(f"iolite Optimiser: Error creating settings path: {e}")

        # Create Hardware Settings Dialog early
        self.settings_dlg = SettingsDialog(self)
        self.ui_initialized = False

    def _get(self, obj, attr):
        if not obj or not hasattr(obj, attr): return None
        try:
            val = getattr(obj, attr)
            return val() if callable(val) else val
        except Exception:
            return None

    def _block_signals(self, obj, block):
        if not obj: return
        try:
            if hasattr(obj, 'blockSignals'):
                func = getattr(obj, 'blockSignals')
                if callable(func):
                    func(block)
                else:
                    setattr(obj, 'blockSignals', block)
        except Exception:
            pass
        

    def on_iolite_data_changed(self, name):
        # Trigger debounce timer
        self.refresh_timer.start()

    def save_settings(self):
        try:
            if self.settings_json_path:
                with open(self.settings_json_path, 'w') as f:
                    json.dump(self.persistent_settings, f)
        except Exception as e:
            IoLog.warning(f"Failed to save settings: {e}")

    def detect_theme(self):
        # Infer theme from window background
        # Reverted to QApplication.palette() as self.palette() is unreliable on startup
        try:
            pal = QApplication.palette()
            bg_color = pal.color(QPalette.Window)
            
            # Simple lightness check (0-255)
            # &lt; 128 is usually dark
            lightness = bg_color.lightness()
            is_dark = lightness < 128
            # IoLog.information(f"Theme Detect: App Lightness={lightness}, IsDark={is_dark}")
            return is_dark
        except:
            return False

    def apply_theme(self, text=None, trigger_opt=True):
        if text is None:
            # Determine which combo sent the signal or default to main
            sender = self.sender()
            if sender == getattr(self, 'combo_theme_spr', None):
                text = self._get(self.combo_theme_spr, 'currentText')
            else:
                text = self._get(self.combo_theme, 'currentText')
        
        # Sync both combos
        for combo_attr in ['combo_theme', 'combo_theme_spr']:
            if hasattr(self, combo_attr):
                cb = getattr(self, combo_attr)
                if self._get(cb, 'currentText') != text:
                    self._block_signals(cb, True)
                    try:
                        cb.setCurrentText(text)
                    except:
                        pass
                    self._block_signals(cb, False)
        IoLog.information(f"iolite Optimiser: apply_theme called with text='{text}'")
        
        # 1. Determine System Theme (Table/UI)
        self.system_is_dark = self.detect_theme()
        
        # 2. Determine Plot Theme (Matplotlib)
        plot_is_dark = False
        if text == 'Auto':
            plot_is_dark = self.system_is_dark
        else:
            plot_is_dark = (text == 'Dark')
            
        style = 'dark_background' if plot_is_dark else 'default'
        
        try:
            # Update global style
            plt.style.use(style)
                    
            # Force update existing figure backgrounds
            # Note: We do NOT pass update_plot here for the main figure, because we handle it explicitly 
            # below with the is_theme_change flag. Calling it here without the flag wipes the zoom.
            for fig_attr, canvas_attr, update_func in [
                ('figure', 'canvas', None), # Main Optimiser Figure - Handled manually below
                ('spr_figure', 'spr_canvas', self.run_spr_analysis)
            ]:
                if hasattr(self, fig_attr):
                    fig = getattr(self, fig_attr)
                    fig.patch.set_facecolor(plt.rcParams['figure.facecolor'])
                    
                    if fig.axes:
                        for ax in fig.axes:
                             ax.set_facecolor(plt.rcParams['axes.facecolor'])
                        
                        # Trigger redraw/data re-plot
                        if update_func:
                            update_func()
                    else:
                        if hasattr(self, canvas_attr):
                            getattr(self, canvas_attr).draw()
            
            # Save setting
            self.persistent_settings['theme'] = text
            self.save_settings()
            
            # 3. Trigger Table Refresh to update highlights
            # This ensures table colours match the newly detected system theme
            if trigger_opt and hasattr(self, 'opt_df') and self.opt_df is not None:
                self.run_optimisation(refresh=False, preserve_zoom=True)
            else:
                self.update_plot(preserve_zoom=True)
            
        except Exception as e:
            IoLog.warning(f"Theme Update Error: {e}")


    def _perform_auto_refresh(self):
        # Only refresh if widget is visible
        if self.isVisible():
             self.run_optimisation(refresh=True)
        
    def showEvent(self, event):
        if not getattr(self, 'ui_initialized', False):
            self.ui_initialized = True
            self.initUI()
            self._on_strategy_changed(self._get(self.cmb_sync_strategy, 'currentText'))
        QWidget.showEvent(self, event)

    def initUI(self):
        # Resize to 80% of screen
        try:
            screen = QApplication.primaryScreen()
            rect = screen.availableGeometry()
            self.resize(int(rect.width() * 0.9), int(rect.height() * 0.9))
        except:
            self.resize(1200, 900) 

        self.setWindowTitle(f"iolite Optimiser - v{VERSION}")

        # Main Layout (Tabs)
        main_layout = QVBoxLayout()
        self.tabs = QTabWidget()
        
        # Tab 1: Determine SPR
        self.tab_spr = QWidget()
        self.init_spr_tab()
        self.tabs.addTab(self.tab_spr, "Determine SPR")
        
        # Tab 2: Method Optimiser
        self.tab_opt = QWidget()
        self.init_optimiser_tab()
        self.tabs.addTab(self.tab_opt, "Method Optimiser")
        
        # Tab 3: Pulse Train Simulator
        if SHOW_PULSE_TRAIN_SIMULATOR:
            IoLog.information("iolite Optimiser: Adding Pulse Train Tab...")
            self.tab_pulse = QWidget()
            self.init_pulse_train_tab()
            self.tabs.addTab(self.tab_pulse, "Pulse Train Simulator")
        
        # Select Optimiser by default for now (or SPR if desired)
        self.tabs.setCurrentIndex(1)
        self.tabs.currentChanged.connect(self._on_tab_changed)
        
        main_layout.addWidget(self.tabs)
        self.setLayout(main_layout)
        
    def _on_tab_changed(self, index):
        if index == 2:
            # Dynamically adjust default pulse shape based on SPR composite peak availability
            has_composite = (hasattr(self, 'spr_raw_results_df') and 
                             self.spr_raw_results_df is not None and 
                             not self.spr_raw_results_df.empty)
            default_shape = "Real Composite Peak" if has_composite else "Model Washout Peak (Lognormal)"
            self.cmb_pulse_shape.setCurrentText(default_shape)

            # 1. Determine target values (default to current values of spinboxes)
            target_rr = self.spin_pulse_rr.value
            target_at_ms = self.spin_pulse_at.value
            target_pulses = self.spin_pulse_count.value
            target_washout = self.spin_pulse_washout.value
            
            is_override = False
            if hasattr(self, 'grp_pulse_override') and self.grp_pulse_override.isChecked():
                is_override = True
                
            # If not in override, pull from last_sync (Optimiser Tab results) or Optimiser tab inputs
            if not is_override:
                if hasattr(self, 'last_sync') and self.last_sync:
                    target_rr = self.last_sync.get('rr_actual', target_rr)
                    target_at_ms = self.last_sync.get('at_actual_s', 0.1) * 1000.0
                    target_pulses = self.last_sync.get('actual_pulses', target_pulses)
                    target_washout = self.spin_wash.value
                else:
                    if hasattr(self, 'spin_init_rr'):
                        target_rr = self.spin_init_rr.value
                    if hasattr(self, 'spin_wash'):
                        target_washout = self.spin_wash.value
                    target_at_ms = getattr(self, 'detected_at_ms', 10.0) if getattr(self, 'detected_at_ms', None) is not None else 10.0
                    target_pulses = target_rr * (target_at_ms / 1000.0)
            
            target_bg_s = self.spin_pulse_bg.value
            target_sig_s = self.spin_pulse_sig.value
            try:
                target_iso = self._get(self.cmb_spr_iso, 'currentText')
            except:
                target_iso = None
            
            # 2. Check if anything has changed since the last run
            changed = False
            last_p = getattr(self, '_last_simulated_params', None)
            if last_p is None:
                changed = True
            else:
                def is_diff(k, val):
                    last_val = last_p.get(k)
                    if last_val is None or val is None:
                        return last_val != val
                    if isinstance(val, str) or isinstance(last_val, str) or isinstance(val, bool) or isinstance(last_val, bool):
                        return last_val != val
                    # Numerical comparison with tolerance
                    try:
                        return abs(float(last_val) - float(val)) > 1e-5
                    except Exception:
                        return last_val != val

                if (is_diff('rr', target_rr) or
                    is_diff('at_ms', target_at_ms) or
                    is_diff('pulses', target_pulses) or
                    is_diff('washout', target_washout) or
                    is_diff('bg_s', target_bg_s) or
                    is_diff('sig_s', target_sig_s) or
                    is_diff('iso', target_iso) or
                    is_diff('is_override', is_override)):
                    changed = True
            
            # 3. If changed (or first time), update UI and trigger recalculation
            if changed:
                self._block_signals(self.spin_pulse_rr, True)
                self._block_signals(self.spin_pulse_at, True)
                self._block_signals(self.spin_pulse_count, True)
                self._block_signals(self.spin_pulse_washout, True)
                
                self.spin_pulse_rr.setValue(target_rr)
                self.spin_pulse_at.setValue(target_at_ms)
                self.spin_pulse_count.setValue(target_pulses)
                self.spin_pulse_washout.setValue(target_washout)
                
                self._block_signals(self.spin_pulse_rr, False)
                self._block_signals(self.spin_pulse_at, False)
                self._block_signals(self.spin_pulse_count, False)
                self._block_signals(self.spin_pulse_washout, False)
                
                self.pulse_debounce_timer.start()
        
    def init_spr_tab(self):
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()

        # --- LEFT COLUMN (CONTROLS) ---
        
        # 0. Reload Button (Fixed at Top)
        h_ctrl = QHBoxLayout()
        self.btn_spr_run = QPushButton("Import SPR File")
        self.btn_spr_run.setFixedHeight(30)
        self.btn_spr_run.clicked.connect(lambda: self.import_and_unload_data(tab="spr"))
        h_ctrl.addWidget(self.btn_spr_run)
        left_layout.addLayout(h_ctrl)

        self.spr_scroll = QScrollArea()
        self.spr_scroll.setWidgetResizable(True)
        self.spr_scroll.setFrameShape(QScrollArea.NoFrame)
        self.spr_scroll.setFixedWidth(460)
        
        scroll_content = QWidget()
        l_settings = QVBoxLayout(scroll_content)

        # 1. Peak Detection Controls
        grp_det = QGroupBox("")
        v_det = QVBoxLayout()
        v_det.setSpacing(4)
        v_det.setContentsMargins(5, 5, 5, 5)
        
        lbl_det_title = QLabel("Peak Detection")
        lbl_det_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_det_title.setAlignment(Qt.AlignCenter)
        v_det.addWidget(lbl_det_title)
        
        # Use 4-column Grid Layout (Cols 2,3 empty)
        grid_det = QGridLayout()
        grid_det.setColumnStretch(3, 1) # Push content left
        
        self.cmb_spr_iso = QComboBox(self)
        self.cmb_spr_iso.currentTextChanged.connect(self._on_spr_iso_changed)
        grid_det.addWidget(QLabel("Select Channel:"), 0, 0)
        grid_det.addWidget(self.cmb_spr_iso, 0, 1)
        
        l_prom = QHBoxLayout()
        l_prom.setContentsMargins(0, 0, 0, 0)
        self.spin_spr_prom = QDoubleSpinBox(self)
        self.spin_spr_prom.setRange(0, 1e9)
        self.spin_spr_prom.setDecimals(0)
        self.spin_spr_prom.setValue(100.0)
        self.spin_spr_prom.setSingleStep(100.0)
        self.spin_spr_prom.setToolTip("Minimum peak height. Prevents detecting noise as peaks. Auto sets to 10 % of maximum data point.")
        self.spin_spr_prom.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.spin_spr_prom.valueChanged.connect(self._on_spr_prom_changed)
        l_prom.addWidget(self.spin_spr_prom)
        
        self.chk_spr_auto_prom = QCheckBox("Auto", self)
        self.chk_spr_auto_prom.setChecked(True)
        self.chk_spr_auto_prom.setToolTip("Minimum peak height. Prevents detecting noise as peaks. Auto sets to 10 % of maximum data point.")
        self.chk_spr_auto_prom.toggled.connect(self._on_spr_auto_prom_toggled)
        l_prom.addWidget(self.chk_spr_auto_prom)
        # self.spin_spr_prom.setEnabled(False) # Removed: User wants always clickable
        
        grid_det.addWidget(QLabel("Peak Cutoff:"), 1, 0)
        grid_det.addLayout(l_prom, 1, 1)
        
        self.spin_spr_dist = QDoubleSpinBox(self)
        self.spin_spr_dist.setRange(0.0001, 100000)
        self.spin_spr_dist.setDecimals(4)
        self.spin_spr_dist.setValue(float(self.persistent_settings.get('spr_min_distance', 0.1)))
        self.spin_spr_dist.setToolTip("Minimum distance between peaks in time. Prevents detecting multiple points on the same peak.")
        self.spin_spr_dist.valueChanged.connect(self._on_spr_dist_changed)
        self.spin_spr_dist.valueChanged.connect(self.save_persistent_settings)
        grid_det.addWidget(QLabel("Min. Distance:"), 2, 0)
        grid_det.addWidget(self.spin_spr_dist, 2, 1)
        
        self.spin_spr_baseline_window = QDoubleSpinBox(self)
        self.spin_spr_baseline_window.setRange(0, 100000)
        self.spin_spr_baseline_window.setDecimals(4)
        self.spin_spr_baseline_window.setValue(float(self.persistent_settings.get('spr_baseline_window', 1000.0)))
        self.spin_spr_baseline_window.setToolTip("Restricts how far the algorithm looks for a baseline minimum in time (0=Inf).")
        self.spin_spr_baseline_window.valueChanged.connect(self._run_spr_analysis_forced)
        self.spin_spr_baseline_window.valueChanged.connect(self.save_persistent_settings)
        
        lbl_baseline_window = QLabel("Baseline Window:")
        grid_det.addWidget(lbl_baseline_window, 6, 0)
        grid_det.addWidget(self.spin_spr_baseline_window, 6, 1)
        
        # Hide the baseline controls per user request, but keep the logic
        lbl_baseline_window.hide()
        self.spin_spr_baseline_window.hide()

        l_smooth = QHBoxLayout()
        l_smooth.setContentsMargins(0, 0, 0, 0)
        
        self.chk_spr_smooth = QCheckBox("Apply Smoothing", self)
        self.chk_spr_smooth.setChecked(self.persistent_settings.get('spr_apply_smooth', False))
        self.chk_spr_smooth.toggled.connect(self._run_spr_analysis_forced)
        self.chk_spr_smooth.toggled.connect(self.save_persistent_settings)
        l_smooth.addWidget(self.chk_spr_smooth)

        self.spin_spr_smooth_window = QDoubleSpinBox(self)
        self.spin_spr_smooth_window.setRange(0.0001, 100)
        self.spin_spr_smooth_window.setDecimals(4)
        self.spin_spr_smooth_window.setValue(float(self.persistent_settings.get('spr_smooth_window', 0.05)))
        self.spin_spr_smooth_window.setToolTip("Window size in time for the median filter to smooth out deep noise valleys.")
        self.spin_spr_smooth_window.valueChanged.connect(self._run_spr_analysis_forced)
        self.spin_spr_smooth_window.valueChanged.connect(self.save_persistent_settings)
        self.spin_spr_smooth_window.setEnabled(self.chk_spr_smooth.isChecked())
        self.chk_spr_smooth.toggled.connect(self.spin_spr_smooth_window.setEnabled)
        l_smooth.addWidget(self.spin_spr_smooth_window)

        lbl_baseline_smoothing = QLabel("Baseline Smoothing:")
        grid_det.addWidget(lbl_baseline_smoothing, 4, 0)
        grid_det.addLayout(l_smooth, 4, 1)
        
        # Hide the smoothing controls per user request
        lbl_baseline_smoothing.hide()
        self.chk_spr_smooth.hide()
        self.spin_spr_smooth_window.hide()
        
        self.cmb_spr_unit = QComboBox(self)
        self.cmb_spr_unit.addItems(["Seconds (s)", "Milliseconds (ms)"])
        # Load selection or default to ms
        self.cmb_spr_unit.setCurrentText(self.persistent_settings.get('spr_time_unit', "Milliseconds (ms)"))
        self.cmb_spr_unit.currentTextChanged.connect(self._update_spr_time_suffixes)
        self.cmb_spr_unit.currentTextChanged.connect(self._run_spr_analysis_forced)
        self.cmb_spr_unit.currentTextChanged.connect(self.save_persistent_settings)
        grid_det.addWidget(QLabel("Time Unit:"), 5, 0)
        
        self._update_spr_time_suffixes()  # Run once to set suffixes and decimals initially
        grid_det.addWidget(self.cmb_spr_unit, 5, 1)
        
        # Row 6: Background sub
        grid_det.addWidget(QLabel("Subtract Background:"), 6, 0)
        self.chk_spr_bg_sub = QCheckBox("", self)
        self.chk_spr_bg_sub.setChecked(self.persistent_settings.get('spr_bg_sub', True))
        self.chk_spr_bg_sub.toggled.connect(self._run_spr_analysis_forced)
        self.chk_spr_bg_sub.toggled.connect(self.save_persistent_settings)
        
        self.chk_spr_auto_bg = QCheckBox("Auto Select Region", self)
        # Default Auto to True
        self.chk_spr_auto_bg.setChecked(self.persistent_settings.get('spr_auto_bg', True))
        self.chk_spr_auto_bg.toggled.connect(self._on_spr_auto_bg_toggled)
        self.chk_spr_auto_bg.toggled.connect(self.save_persistent_settings)
        
        self.spin_spr_bg_start = QDoubleSpinBox(self)
        self.spin_spr_bg_start.setRange(0, 10000)
        self.spin_spr_bg_start.setDecimals(2)
        self.spin_spr_bg_start.setFixedWidth(60)
        self.spin_spr_bg_start.setValue(float(self.persistent_settings.get('spr_bg_start', 2.0)))
        self.spin_spr_bg_start.valueChanged.connect(self._on_spr_bg_changed)
        
        self.spin_spr_bg_end = QDoubleSpinBox(self)
        self.spin_spr_bg_end.setRange(0, 10000)
        self.spin_spr_bg_end.setDecimals(2)
        self.spin_spr_bg_end.setFixedWidth(60)
        self.spin_spr_bg_end.setValue(float(self.persistent_settings.get('spr_bg_end', 8.0)))
        self.spin_spr_bg_end.valueChanged.connect(self._on_spr_bg_changed)
        
        # We don't need enable logic, we just use visibility logic (in the right panel setup)
        
        grid_det.addWidget(self.chk_spr_bg_sub, 6, 1)
        
        h_excl = QHBoxLayout()
        h_excl.setContentsMargins(0, 0, 0, 0)
        
        self.chk_spr_auto_exclude = QCheckBox("Auto-Exclude Wide Outliers (Percentile >)", self)
        self.chk_spr_auto_exclude.setChecked(False)
        self.chk_spr_auto_exclude.toggled.connect(self._on_spr_auto_exclude_toggled)
        
        self.spin_spr_auto_exclude_pct = QDoubleSpinBox(self)
        self.spin_spr_auto_exclude_pct.setRange(50.0, 99.9)
        self.spin_spr_auto_exclude_pct.setDecimals(1)
        self.spin_spr_auto_exclude_pct.setValue(90.0)
        self.spin_spr_auto_exclude_pct.setSingleStep(1.0)
        self.spin_spr_auto_exclude_pct.setSuffix(" %")
        self.spin_spr_auto_exclude_pct.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.spin_spr_auto_exclude_pct.valueChanged.connect(self._on_spr_auto_exclude_pct_changed)
        
        h_excl.addWidget(self.chk_spr_auto_exclude)
        h_excl.addWidget(self.spin_spr_auto_exclude_pct)
        grid_det.addLayout(h_excl, 3, 0, 1, 2)
        
        v_det.addLayout(grid_det)
        grp_det.setLayout(v_det)
        l_settings.addWidget(grp_det)
        
        # 2. Results Metrics
        self.grp_spr_res = QGroupBox("")
        l_res = QVBoxLayout()
        l_res.setSpacing(4)
        l_res.setContentsMargins(5, 5, 5, 5)
        
        lbl_res_main_title = QLabel("SPR Analysis Results")
        lbl_res_main_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_res_main_title.setAlignment(Qt.AlignCenter)
        l_res.addWidget(lbl_res_main_title)
        
        # FW 0.1M (10%) Group
        grp_10 = QGroupBox("")
        v_main_10 = QVBoxLayout()
        v_main_10.setSpacing(4)
        v_main_10.setContentsMargins(5, 5, 5, 5)
        
        lbl_title_10 = QLabel("FW 0.1M (10%)")
        lbl_title_10.setStyleSheet("font-weight: bold; font-size: 9 pt; padding: 0px; margin: 0px;")
        lbl_title_10.setAlignment(Qt.AlignCenter)
        v_main_10.addWidget(lbl_title_10)
        
        grid_res_10 = QGridLayout()
        grid_res_10.setSpacing(2)
        grid_res_10.setContentsMargins(0, 0, 0, 0)
        grid_res_10.setColumnMinimumWidth(0, 75)
        grid_res_10.setColumnStretch(1, 1)
        grid_res_10.setColumnMinimumWidth(2, 75)
        grid_res_10.setColumnStretch(3, 1)
        grid_res_10.setColumnStretch(4, 0)
        
        # Row 0: Average & RSD
        lbl_fw10_label = QLabel("Average:")
        lbl_fw10_label.setStyleSheet("font-size: 8 pt;")
        lbl_fw10_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_fw10 = QLabel("- s")
        self.lbl_spr_fw10.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_fw10.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        lbl_rsd10_label = QLabel("RSD:")
        lbl_rsd10_label.setStyleSheet("font-size: 8 pt;")
        lbl_rsd10_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_rsd10 = QLabel("- %")
        self.lbl_spr_rsd10.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_rsd10.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        grid_res_10.addWidget(lbl_fw10_label, 0, 0)
        grid_res_10.addWidget(self.lbl_spr_fw10, 0, 1)
        grid_res_10.addWidget(lbl_rsd10_label, 0, 2)
        grid_res_10.addWidget(self.lbl_spr_rsd10, 0, 3)
        
        # Row 1: Area & Area RSD
        lbl_area10_label = QLabel("Area:")
        lbl_area10_label.setStyleSheet("font-size: 8 pt;")
        lbl_area10_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_area10 = QLabel("- Counts")
        self.lbl_spr_area10.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_area10.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        lbl_area_rsd10_label = QLabel("RSD:")
        lbl_area_rsd10_label.setStyleSheet("font-size: 8 pt;")
        lbl_area_rsd10_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_area_rsd10 = QLabel("- %")
        self.lbl_spr_area_rsd10.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_area_rsd10.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        grid_res_10.addWidget(lbl_area10_label, 1, 0)
        grid_res_10.addWidget(self.lbl_spr_area10, 1, 1)
        grid_res_10.addWidget(lbl_area_rsd10_label, 1, 2)
        grid_res_10.addWidget(self.lbl_spr_area_rsd10, 1, 3)
        
        # Row 2: Max (Left) & Composite (Right)
        lbl_max10_label = QLabel("Max:")
        lbl_max10_label.setStyleSheet("font-size: 8 pt;")
        lbl_max10_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_fw_max10 = QLabel("- s")
        self.lbl_spr_fw_max10.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_fw_max10.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        lbl_comp10_label = QLabel("Composite:")
        lbl_comp10_label.setStyleSheet("font-size: 8 pt;")
        lbl_comp10_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_fw_comp10 = QLabel("- s")
        self.lbl_spr_fw_comp10.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_fw_comp10.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        grid_res_10.addWidget(lbl_max10_label, 2, 0)
        grid_res_10.addWidget(self.lbl_spr_fw_max10, 2, 1)
        grid_res_10.addWidget(lbl_comp10_label, 2, 2)
        grid_res_10.addWidget(self.lbl_spr_fw_comp10, 2, 3)
        
        h_btns_10 = QHBoxLayout()
        self.btn_apply10_avg = QPushButton("Apply Average")
        self.btn_apply10_avg.clicked.connect(lambda: self._on_spr_apply_clicked("10avg"))
        self.btn_apply10_comp = QPushButton("Apply Composite")
        self.btn_apply10_comp.clicked.connect(lambda: self._on_spr_apply_clicked("10comp"))
        self.btn_apply10_max = QPushButton("Apply Max")
        self.btn_apply10_max.clicked.connect(lambda: self._on_spr_apply_clicked("10max"))
        h_btns_10.addWidget(self.btn_apply10_avg)
        h_btns_10.addWidget(self.btn_apply10_comp)
        h_btns_10.addWidget(self.btn_apply10_max)
        
        v_main_10.addLayout(grid_res_10)
        v_main_10.addSpacing(2)
        v_main_10.addLayout(h_btns_10)
        grp_10.setLayout(v_main_10)
        l_res.addWidget(grp_10)
        
        # FW 0.01M (1%) Group
        grp_1 = QGroupBox("")
        v_main_1 = QVBoxLayout()
        v_main_1.setSpacing(4)
        v_main_1.setContentsMargins(5, 5, 5, 5)
        
        lbl_title_1 = QLabel("FW 0.01M (1%)")
        lbl_title_1.setStyleSheet("font-weight: bold; font-size: 9 pt; padding: 0px; margin: 0px;")
        lbl_title_1.setAlignment(Qt.AlignCenter)
        v_main_1.addWidget(lbl_title_1)
        
        grid_res_1 = QGridLayout()
        grid_res_1.setSpacing(2)
        grid_res_1.setContentsMargins(0, 0, 0, 0)
        grid_res_1.setColumnMinimumWidth(0, 75)
        grid_res_1.setColumnStretch(1, 1)
        grid_res_1.setColumnMinimumWidth(2, 75)
        grid_res_1.setColumnStretch(3, 1)
        grid_res_1.setColumnStretch(4, 0)
        
        # Row 0: Average & RSD
        lbl_fw1_label = QLabel("Average:")
        lbl_fw1_label.setStyleSheet("font-size: 8 pt;")
        lbl_fw1_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_fw1 = QLabel("- s")
        self.lbl_spr_fw1.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_fw1.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        lbl_rsd1_label = QLabel("RSD:")
        lbl_rsd1_label.setStyleSheet("font-size: 8 pt;")
        lbl_rsd1_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_rsd1 = QLabel("- %")
        self.lbl_spr_rsd1.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_rsd1.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        grid_res_1.addWidget(lbl_fw1_label, 0, 0)
        grid_res_1.addWidget(self.lbl_spr_fw1, 0, 1)
        grid_res_1.addWidget(lbl_rsd1_label, 0, 2)
        grid_res_1.addWidget(self.lbl_spr_rsd1, 0, 3)
        
        # Row 1: Area & Area RSD
        lbl_area1_label = QLabel("Area:")
        lbl_area1_label.setStyleSheet("font-size: 8 pt;")
        lbl_area1_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_area1 = QLabel("- Counts")
        self.lbl_spr_area1.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_area1.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        lbl_area_rsd1_label = QLabel("RSD:")
        lbl_area_rsd1_label.setStyleSheet("font-size: 8 pt;")
        lbl_area_rsd1_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_area_rsd1 = QLabel("- %")
        self.lbl_spr_area_rsd1.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_area_rsd1.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        grid_res_1.addWidget(lbl_area1_label, 1, 0)
        grid_res_1.addWidget(self.lbl_spr_area1, 1, 1)
        grid_res_1.addWidget(lbl_area_rsd1_label, 1, 2)
        grid_res_1.addWidget(self.lbl_spr_area_rsd1, 1, 3)
        
        # Row 2: Max (Left) & Composite (Right)
        lbl_max1_label = QLabel("Max:")
        lbl_max1_label.setStyleSheet("font-size: 8 pt;")
        lbl_max1_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_fw_max1 = QLabel("- s")
        self.lbl_spr_fw_max1.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_fw_max1.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        lbl_comp1_label = QLabel("Composite:")
        lbl_comp1_label.setStyleSheet("font-size: 8 pt;")
        lbl_comp1_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.lbl_spr_fw_comp1 = QLabel("- s")
        self.lbl_spr_fw_comp1.setStyleSheet("font-size: 8 pt;")
        self.lbl_spr_fw_comp1.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        grid_res_1.addWidget(lbl_max1_label, 2, 0)
        grid_res_1.addWidget(self.lbl_spr_fw_max1, 2, 1)
        grid_res_1.addWidget(lbl_comp1_label, 2, 2)
        grid_res_1.addWidget(self.lbl_spr_fw_comp1, 2, 3)
        
        h_btns_1 = QHBoxLayout()
        self.btn_apply1_avg = QPushButton("Apply Average")
        self.btn_apply1_avg.clicked.connect(lambda: self._on_spr_apply_clicked("1avg"))
        self.btn_apply1_comp = QPushButton("Apply Composite")
        self.btn_apply1_comp.clicked.connect(lambda: self._on_spr_apply_clicked("1comp"))
        self.btn_apply1_max = QPushButton("Apply Max")
        self.btn_apply1_max.clicked.connect(lambda: self._on_spr_apply_clicked("1max"))
        h_btns_1.addWidget(self.btn_apply1_avg)
        h_btns_1.addWidget(self.btn_apply1_comp)
        h_btns_1.addWidget(self.btn_apply1_max)
        
        v_main_1.addLayout(grid_res_1)
        v_main_1.addSpacing(2)
        v_main_1.addLayout(h_btns_1)
        grp_1.setLayout(v_main_1)
        l_res.addWidget(grp_1)
        
        self.grp_spr_res.setLayout(l_res)
        l_settings.addWidget(self.grp_spr_res)
        
        # 3. Maximum Observed SPR Panel
        self.grp_spr_max = QGroupBox("")
        l_max = QVBoxLayout()
        l_max.setSpacing(4)
        l_max.setContentsMargins(5, 5, 5, 5)
        
        lbl_max_main_title = QLabel("Maximum Observed SPR")
        lbl_max_main_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_max_main_title.setAlignment(Qt.AlignCenter)
        l_max.addWidget(lbl_max_main_title)
        
        self.lbl_spr_max_iso = QLabel("Channel: -")
        self.lbl_spr_max_iso.setStyleSheet("font-weight: bold; font-size: 9pt;")
        self.lbl_spr_max_iso.setAlignment(Qt.AlignCenter)
        l_max.addWidget(self.lbl_spr_max_iso)
        
        # Grid for Max Stats (10% and 1%)
        grid_max_stats = QGridLayout()
        grid_max_stats.setSpacing(2)
        grid_max_stats.setContentsMargins(0, 0, 0, 0)
        
        # Helper for headers
        def _head(txt):
            lbl = QLabel(txt)
            lbl.setStyleSheet("font-weight: bold; font-size: 8pt;")
            lbl.setAlignment(Qt.AlignCenter)
            return lbl
            
        grid_max_stats.addWidget(_head("Metric"), 0, 0)
        grid_max_stats.addWidget(_head("FW 0.1M (10%)"), 0, 1)
        grid_max_stats.addWidget(_head("FW 0.01M (1%)"), 0, 2)
        
        self.lbl_spr_max_avg10 = QLabel("-"); self.lbl_spr_max_avg10.setAlignment(Qt.AlignCenter)
        self.lbl_spr_max_avg1 = QLabel("-"); self.lbl_spr_max_avg1.setAlignment(Qt.AlignCenter)
        grid_max_stats.addWidget(QLabel("Average:"), 1, 0)
        grid_max_stats.addWidget(self.lbl_spr_max_avg10, 1, 1)
        grid_max_stats.addWidget(self.lbl_spr_max_avg1, 1, 2)
        
        self.lbl_spr_max_max10 = QLabel("-"); self.lbl_spr_max_max10.setAlignment(Qt.AlignCenter)
        self.lbl_spr_max_max1 = QLabel("-"); self.lbl_spr_max_max1.setAlignment(Qt.AlignCenter)
        grid_max_stats.addWidget(QLabel("Max:"), 2, 0)
        grid_max_stats.addWidget(self.lbl_spr_max_max10, 2, 1)
        grid_max_stats.addWidget(self.lbl_spr_max_max1, 2, 2)
        
        self.lbl_spr_max_comp10 = QLabel("-"); self.lbl_spr_max_comp10.setAlignment(Qt.AlignCenter)
        self.lbl_spr_max_comp1 = QLabel("-"); self.lbl_spr_max_comp1.setAlignment(Qt.AlignCenter)
        grid_max_stats.addWidget(QLabel("Composite:"), 3, 0)
        grid_max_stats.addWidget(self.lbl_spr_max_comp10, 3, 1)
        grid_max_stats.addWidget(self.lbl_spr_max_comp1, 3, 2)
        
        l_max.addLayout(grid_max_stats)
        
        # Buttons
        h_btns_max = QHBoxLayout()
        self.btn_apply_max_avg = QPushButton("Apply Average")
        self.btn_apply_max_avg.clicked.connect(lambda: self._on_spr_apply_max_clicked("avg"))
        self.btn_apply_max_comp = QPushButton("Apply Composite")
        self.btn_apply_max_comp.clicked.connect(lambda: self._on_spr_apply_max_clicked("comp"))
        self.btn_apply_max_val = QPushButton("Apply Max")
        self.btn_apply_max_val.clicked.connect(lambda: self._on_spr_apply_max_clicked("max"))
        h_btns_max.addWidget(self.btn_apply_max_avg)
        h_btns_max.addWidget(self.btn_apply_max_comp)
        h_btns_max.addWidget(self.btn_apply_max_val)
        l_max.addLayout(h_btns_max)
        
        self.chk_spr_show_all_stats = QCheckBox("Show All Channel Stats")
        self.chk_spr_show_all_stats.setChecked(False)
        self.chk_spr_show_all_stats.toggled.connect(self._on_spr_show_all_stats_toggled)
        l_max.addWidget(self.chk_spr_show_all_stats)
        
        self.grp_spr_max.setLayout(l_max)
        l_settings.addWidget(self.grp_spr_max)
        
        l_settings.addStretch()
        self.spr_scroll.setWidget(scroll_content)
        left_layout.addWidget(self.spr_scroll)
        
        # --- RIGHT COLUMN (PLOT & TABLE) ---
        h_ctrl_spr1 = QHBoxLayout()
        h_ctrl_spr2 = QHBoxLayout()
        
        self.combo_theme_spr = QComboBox()
        self.combo_theme_spr.addItems(["Auto", "Dark", "Light"])
        self.combo_theme_spr.setCurrentText(self.persistent_settings.get('theme', 'Auto'))
        self.combo_theme_spr.currentTextChanged.connect(self.apply_theme)
        
        self.chk_y_zoom_spr = QCheckBox("Pan / Zoom Y")
        self.chk_y_zoom_spr.setChecked(self.persistent_settings.get('spr_y_zoom', False))
        self.chk_y_zoom_spr.toggled.connect(self._on_spr_y_zoom_toggled)
        
        self.chk_rescale_spr = QCheckBox("Auto-Rescale Y")
        self.chk_rescale_spr.setChecked(True) # Ephemeral: Always default to ON
        self.chk_rescale_spr.toggled.connect(self._on_spr_rescale_toggled)
        
        h_ctrl_spr1.addWidget(QLabel("Theme:"))
        h_ctrl_spr1.addWidget(self.combo_theme_spr)
        h_ctrl_spr1.addSpacing(10)
        h_ctrl_spr1.addWidget(self.chk_y_zoom_spr)
        h_ctrl_spr1.addWidget(self.chk_rescale_spr)
        h_ctrl_spr1.addStretch()
        
        h_ctrl_spr2.addWidget(self.chk_spr_auto_bg)
        self.lbl_bg_prefix = QLabel("Background (s):", self)
        h_ctrl_spr2.addWidget(self.lbl_bg_prefix)
        h_ctrl_spr2.addWidget(self.spin_spr_bg_start)
        self.lbl_bg_to = QLabel("to", self)
        h_ctrl_spr2.addWidget(self.lbl_bg_to)
        h_ctrl_spr2.addWidget(self.spin_spr_bg_end)
        h_ctrl_spr2.addStretch()
        
        def _toggle_bg_layout_vis(visible):
            self.chk_spr_auto_bg.setVisible(visible)
            self.lbl_bg_prefix.setVisible(visible)
            self.spin_spr_bg_start.setVisible(visible)
            self.lbl_bg_to.setVisible(visible)
            self.spin_spr_bg_end.setVisible(visible)
            
        self.chk_spr_bg_sub.toggled.connect(_toggle_bg_layout_vis)
        # Init
        _toggle_bg_layout_vis(self.chk_spr_bg_sub.isChecked())
        
        # Hide Summary Stats Checkbox (Right Aligned)
        self.chk_hide_stats = QCheckBox("Hide Summary Stats")
        self.chk_hide_stats.setChecked(False)
        self.chk_hide_stats.toggled.connect(self._on_stats_toggled)
        h_ctrl_spr1.addWidget(self.chk_hide_stats)
        
        right_layout.addLayout(h_ctrl_spr1)
        right_layout.addLayout(h_ctrl_spr2)

        self.spr_figure = Figure(figsize=(5, 4), dpi=96, constrained_layout=False)
        self.spr_canvas = FigureCanvas(self.spr_figure)
        self.spr_canvas.mpl_connect('resize_event', self._on_plot_resize)
        self.spr_canvas.mpl_connect('scroll_event', self.on_zoom)
        self.spr_canvas.mpl_connect('button_press_event', self.on_press)
        self.spr_canvas.mpl_connect('button_release_event', self.on_release)
        self.spr_canvas.mpl_connect('motion_notify_event', self.on_drag)
        self.spr_canvas.mpl_connect('pick_event', self.on_pick)
        right_layout.addWidget(self.spr_canvas, 1)

        self.lbl_spr_table_title = QLabel("SPR Peaks Detected - 0 | SPR Peaks Excluded - 0")
        self.lbl_spr_table_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        self.lbl_spr_table_title.setAlignment(Qt.AlignCenter)
        right_layout.addWidget(self.lbl_spr_table_title)
        
        self.spr_table_stack = QStackedWidget()
        
        self.spr_table = CopyableTableWidget()
        self.spr_table.setColumnCount(7)
        self.spr_table.setHorizontalHeaderLabels(["Peak", "Time (s)", "FW0.1M (10 %) (s)", "FW0.1M (10 %) Area", "FW0.01M (1 %) (s)", "FW0.01M (1 %) Area", "Max Intensity"])
        self.spr_table.verticalHeader().setVisible(False)
        self.spr_table.horizontalHeader().setDefaultAlignment(Qt.AlignCenter)
        self.spr_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.spr_table_stack.addWidget(self.spr_table)
        
        self.spr_all_stats_table = CopyableTableWidget()
        self.spr_all_stats_table.setColumnCount(11)
        self.spr_all_stats_table.setHorizontalHeaderLabels([
            "Channel", "Avg (10%)", "Max (10%)", "Comp (10%)", "RSD (10%)", "Area 10% (Avg)", 
            "Avg (1%)", "Max (1%)", "Comp (1%)", "RSD (1%)", "Area 1% (Avg)"
        ])
        self.spr_all_stats_table.verticalHeader().setVisible(False)
        self.spr_all_stats_table.horizontalHeader().setDefaultAlignment(Qt.AlignCenter)
        self.spr_all_stats_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.spr_table_stack.addWidget(self.spr_all_stats_table)
        
        right_layout.addWidget(self.spr_table_stack, 1)
        
        main_layout.addLayout(left_layout)
        main_layout.addLayout(right_layout)
        self.tab_spr.setLayout(main_layout)

        # Exclusion tracking
        self.spr_excluded_peaks = {} # isotope -> set of peak indices
        self.spr_peak_label_map = {} # artist -> peak index
        
        # Caching & Per-Channel States
        self.spr_channel_prominence = {} # isotope -> user's preferred prominence
        self.spr_channel_auto_prom = {} # isotope -> whether it uses auto prominence
        self.spr_auto_exclude = False
        self.spr_auto_exclude_pct = 90.0
        self.spr_all_raw_results = {} # isotope -> cached raw_results df
        self.spr_all_raw_results = {} # isotope -> cached raw_results df

    def _run_spr_analysis_forced(self, *args, **kwargs):
        self.run_spr_analysis(force_refresh=True)

    def run_spr_analysis(self, *args, rescale=None, **kwargs):
        if self.spr_df is None or find_peaks is None:
            return
        
        # Capture current limits if we don't want to rescale
        old_xlim, old_ylim = None, None
        if rescale is False and hasattr(self, 'spr_figure') and self.spr_figure.axes:
            old_ax = self.spr_figure.axes[0]
            old_xlim = old_ax.get_xlim()
            old_ylim = old_ax.get_ylim()
            
        # Auto-detect prominence if enabled
        if self._get(self.chk_spr_auto_prom, 'isChecked'):
            self._auto_detect_spr_prominence()
            
        # Auto-detect background if enabled
        if self._get(getattr(self, 'chk_spr_auto_bg', None), 'isChecked'):
            self._auto_detect_spr_bg()
            
        iso = self._get(self.cmb_spr_iso, 'currentText')
        prominence = self._get(self.spin_spr_prom, 'value')
        unit_text = self._get(self.cmb_spr_unit, 'currentText')
        distance = self._get(self.spin_spr_dist, 'value')
        wlen_val = self._get(self.spin_spr_baseline_window, 'value')
        wlen = None if wlen_val == 0 else wlen_val
        apply_smoothing = self._get(self.chk_spr_smooth, 'isChecked')
        smooth_window = self._get(self.spin_spr_smooth_window, 'value')
        
        if iso is None or prominence is None or unit_text is None:
            return

        # detect if CPS
        is_cps = True
        if hasattr(self, 'channel_is_cps') and self.channel_is_cps.get(iso, False):
            is_cps = True
        elif getattr(self, 'cached_unit_label', "Counts") == "Counts":
            is_cps = False

        # Restore missing variables
        is_ms = "ms" in unit_text
        mult = 1000.0 if is_ms else 1.0
        unit_label = "ms" if is_ms else "s"
        prec = 2 if is_ms else 4

        # Convert ui values to Seconds since analyse_washout_peaks requires seconds
        distance_s = distance / mult
        wlen_s = wlen / mult if wlen is not None else None
        smooth_window_s = smooth_window / mult

        bg_sub = self._get(self.chk_spr_bg_sub, 'isChecked')
        bg_start_s = self._get(self.spin_spr_bg_start, 'value')
        bg_end_s = self._get(self.spin_spr_bg_end, 'value')

        # 1. Detect All Peaks (Cache or find new)
        # We now run the global update first which caches raw dataframe results for all channels
        # If force_refresh is passed, we wipe the cache first to ensure new percentiles are calculated
        force_refresh = kwargs.get('force_refresh', False)
        if force_refresh and hasattr(self, 'spr_all_raw_results'):
            self.spr_all_raw_results.clear()
            
        if not getattr(self, '_is_updating_all', False) and (force_refresh or not hasattr(self, 'spr_all_raw_results') or not self.spr_all_raw_results):
            self._is_updating_all = True
            self._update_all_channel_stats()
            self._is_updating_all = False
            
        if hasattr(self, 'spr_all_raw_results') and iso in self.spr_all_raw_results:
             raw_result = self.spr_all_raw_results[iso]
        else:
             raw_result = (None, {"Error": "Channel data not available in cache."})
        
        if isinstance(raw_result, tuple):
             # Handle Error or Empty Count return (df, dict)
             self.spr_raw_results_df, stats = raw_result
             # If it's a tuple, we have no peaks to filter, so we skip step 2
             filtered_df = self.spr_raw_results_df 
             excluded = set()
        else:
             # Standard DataFrame return
             self.spr_raw_results_df = raw_result
        
             # 2. Filter Excluded Peaks for Statistics
             if iso not in self.spr_excluded_peaks:
                 self.spr_excluded_peaks[iso] = set()
                 
             excluded = self.spr_excluded_peaks[iso].copy()
             
             # Auto-exclusion filter
             auto_excl = getattr(self, 'spr_auto_exclude', False)
             if auto_excl and len(self.spr_raw_results_df) > 5:
                 pct_val = getattr(self, 'spr_auto_exclude_pct', 90.0)
                 df_temp = self.spr_raw_results_df
                 fw = df_temp['r001'] - df_temp['l001']
                 fw_p_upper = np.percentile(fw, pct_val)
                 
                 outliers = df_temp[fw > fw_p_upper]['Peak Index'].values
                 excluded.update(outliers)
                 
             filtered_df = self.spr_raw_results_df[~self.spr_raw_results_df['Peak Index'].isin(excluded)]
                 
             # 3. Summarize
             stats = Logic.summarize_peaks(filtered_df)
        
        # Update Table Title
        num_detected = len(self.spr_raw_results_df)
        num_excluded = len(excluded)
        self.lbl_spr_table_title.setText(f"SPR Peaks Detected - {num_detected} | SPR Peaks Excluded - {num_excluded}")
        
        # Calculate Time Offset for Zeroing
        t_orig = self.spr_df['Time'].values
        t0 = t_orig[0] if len(t_orig) > 0 else 0
        t_zeroed = t_orig - t0 if len(t_orig) > 0 else t_orig
        
        # Update Metrics
        # Helper for SI Formatting
        def _si(val):
            if val >= 1e9: return f"{val/1e9:.2f}", " G"
            if val >= 1e6: return f"{val/1e6:.2f}", " M"
            if val >= 1e3: return f"{val/1e3:.2f}", " k"
            return f"{val:.2f}", " "

        # Helper for Adaptive Precision (User Request)
        def _adap(v):
            if v >= 100: return f"{v:.0f}"
            if v >= 10:  return f"{v:.1f}"
            if v >= 1:   return f"{v:.2f}"
            return f"{v:.3f}"
            
        # Helper for adaptive rounding (float)
        def _adap_round(v_ms):
            if v_ms >= 100: return round(v_ms, 0)
            if v_ms >= 10:  return round(v_ms, 1)
            return round(v_ms, 2)

        has_peaks = True
        error_msg = None
        
        # Update Metrics
        if "Error" in stats:
            has_peaks = False
            error_msg = stats["Error"]
            self.lbl_spr_fw10.setText("Error")
            self.lbl_spr_fw1.setText(error_msg)
            
        elif stats.get("Count", 0) == 0:
            has_peaks = False
            self.lbl_spr_fw10.setText("-")
            self.lbl_spr_rsd10.setText("-")
            self.lbl_spr_fw_max10.setText("-")
            self.lbl_spr_fw1.setText("No peaks detected")
            self.lbl_spr_rsd1.setText("-")
            self.lbl_spr_fw_max1.setText("-")
            
            # Clear stored values
            self._last_spr_fw10_avg_ms = 0
            self._last_spr_fw10_max_ms = 0
            self._last_spr_fw1_avg_ms = 0
            self._last_spr_fw1_max_ms = 0

        if has_peaks:
            # Forced "Counts" unit for Area
            area_unit_label = " Counts"

            self.lbl_spr_fw10.setText(f"<b>{_adap(stats['FW0.1M Mean']*mult)}</b> {unit_label}")
            self.lbl_spr_rsd10.setText(f"{stats['FW0.1M RSD']:.1f} %")
            self.lbl_spr_fw_max10.setText(f"<b>{_adap(stats['FW0.1M Max']*mult)}</b> {unit_label}")
            # Use SI formatting for Area
            s_mean10, p_mean10 = _si(stats['Area 10% Mean'])
            self.lbl_spr_area10.setText(f"<b>{s_mean10}</b>{p_mean10}{area_unit_label}")
            self.lbl_spr_area_rsd10.setText(f"{stats['Area 10% RSD']:.1f} %")

            self.lbl_spr_fw1.setText(f"<b>{_adap(stats['FW0.01M Mean']*mult)}</b> {unit_label}")
            self.lbl_spr_rsd1.setText(f"{stats['FW0.01M RSD']:.1f} %")
            self.lbl_spr_fw_max1.setText(f"<b>{_adap(stats['FW0.01M Max']*mult)}</b> {unit_label}")
            # Use SI formatting for Area
            s_mean1, p_mean1 = _si(stats['Area 1% Mean'])
            self.lbl_spr_area1.setText(f"<b>{s_mean1}</b>{p_mean1}{area_unit_label}")
            self.lbl_spr_area_rsd1.setText(f"{stats['Area 1% RSD']:.1f} %")
            
            # Store for "Apply" logic (rounded to matching UI precision)
            self._last_spr_fw10_avg_ms = _adap_round(stats['FW0.1M Mean'] * 1000.0)
            self._last_spr_fw10_max_ms = _adap_round(stats['FW0.1M Max'] * 1000.0)
            self._last_spr_fw1_avg_ms = _adap_round(stats['FW0.01M Mean'] * 1000.0)
            self._last_spr_fw1_max_ms = _adap_round(stats['FW0.01M Max'] * 1000.0)
        else:
             # Ensure labels that weren't set in the error block are cleared
             if self.lbl_spr_area10.text != "-": self.lbl_spr_area10.setText("-")
             if self.lbl_spr_area_rsd10.text != "-": self.lbl_spr_area_rsd10.setText("-")
             if self.lbl_spr_area1.text != "-": self.lbl_spr_area1.setText("-")
             if self.lbl_spr_area_rsd1.text != "-": self.lbl_spr_area_rsd1.setText("-")

        # Update Table (Dynamic Units)
        headers = [
            "Peak", 
            "Time (s)",  # Always in seconds
            f"FW0.1M ({unit_label})", 
            "FW0.1M (Counts)", 
            f"FW0.01M ({unit_label})", 
            "FW0.01M (Counts)", 
            "Max Intensity"
        ]
        self.spr_table.setHorizontalHeaderLabels(headers)
        self.spr_table.setRowCount(0)
        self.spr_table.setRowCount(len(self.spr_raw_results_df))
        
        # excluded = self.persistent_settings.get('spr_excluded_peaks', {}).get(iso, []) # REMOVED: Uses runtime variable 'excluded' defined above

        for i, row in self.spr_raw_results_df.iterrows():
            p_idx = int(row['Peak Index'])
            is_excluded = p_idx in excluded
            
            # Helper to create item
            def _it(txt):
                it = QTableWidgetItem(txt)
                it.setTextAlignment(Qt.AlignCenter)
                if is_excluded:
                     it.setForeground(QColor("gray"))
                     f = it.font(); f.setStrikeOut(True); it.setFont(f)
                return it

            self.spr_table.setItem(i, 0, _it(str(p_idx)))
            self.spr_table.setItem(i, 1, _it(f"{row['Peak Time (s)'] - t0:.2f}"))
            
            # Col 2: FW0.1M (scaled & adaptive)
            self.spr_table.setItem(i, 2, _it(f"{_adap(row['FW0.1M (s)'] * mult)}"))
            
            # Col 3: FW0.1M Area (stays as is relative to key, but formatted)
            s_a01, p_a01 = _si(row['Area 10%'])
            self.spr_table.setItem(i, 3, _it(f"{s_a01}{p_a01}"))
            
            # Col 4: FW0.01M (scaled & adaptive)
            self.spr_table.setItem(i, 4, _it(f"{_adap(row['FW0.01M (s)'] * mult)}"))

            # Col 5: FW0.01M Area
            s_a001, p_a001 = _si(row['Area 1%'])
            self.spr_table.setItem(i, 5, _it(f"{s_a001}{p_a001}"))
            
            # Col 6: Max Intensity (SI)
            s_mx, p_mx = _si(row['Max Intensity'])
            self.spr_table.setItem(i, 6, _it(f"{s_mx}{p_mx}"))

        
        # Plot
        self.spr_figure.clear()
        
        # Explicitly set figure facecolor from current rcParams
        self.spr_figure.patch.set_facecolor(plt.rcParams['figure.facecolor'])
        
        # Determine layout (80:20 split if composite enabled)
        show_composite = True # Always show composite peak (User Request)
        fg_col = plt.rcParams['axes.labelcolor']
        
        if show_composite:
            gs = self.spr_figure.add_gridspec(1, 2, width_ratios=[8, 2], wspace=0.02)
            ax = self.spr_figure.add_subplot(gs[0, 0])
            self.spr_ax_comp = self.spr_figure.add_subplot(gs[0, 1])
        else:
            ax = self.spr_figure.add_subplot(111)
            self.spr_ax_comp = None
        
        # Theme-Aware Mega-Palette (Purple Index 8 Start) - Not persistent for SPR
        safe_palette = self._get_theme_safe_palette(is_opt=False)
        ax.set_prop_cycle(cycler(color=safe_palette))
        ax.set_facecolor(plt.rcParams['axes.facecolor'])

        # Get theme-aware foreground color
        fg_col = plt.rcParams['axes.labelcolor']

        # --- BACKGROUND SUBTRACTION FOR PLOT ---
        y_plot = self.spr_df[iso].values.copy()
        if bg_sub:
            bg_mask = (t_zeroed >= bg_start_s) & (t_zeroed <= bg_end_s)
            if bg_mask.any():
                bg_mean = np.mean(y_plot[bg_mask])
                y_plot = np.clip(y_plot - bg_mean, 0, None)
            
            if hasattr(self, 'cached_is_dark'):
                is_dark = self.cached_is_dark
            else:
                is_dark = plt.rcParams['axes.facecolor'] not in ['white', '#ffffff', 'w']
            bg_col = '#808080' if is_dark else '#9e9e9e' # Slate/Medium Grey (more visible)
            v_alpha = 0.35 if is_dark else 0.25

            self.spr_bg_span = ax.axvspan(bg_start_s, bg_end_s, color=bg_col, alpha=v_alpha, lw=0)
        else:
            self.spr_bg_span = None

        # Plot individual channel (Always seconds X)
        # Uses first color from palette (now Purple)
        line, = ax.plot(t_zeroed, y_plot, lw=1.5, alpha=0.9, label=iso)
        
        # Peak Index Labels (Always seconds X)
        self.spr_peak_label_map = {}
        max_y = self.spr_raw_results_df['Max Intensity'].max() if not self.spr_raw_results_df.empty else 1.0
        y_offset = max_y * 0.02
        
        for i, row in self.spr_raw_results_df.iterrows():
             p_idx = int(row['Peak Index'])
             is_excluded = p_idx in excluded
             alpha = 0.3 if is_excluded else 1.0
             t = ax.text((row['Peak Time (s)']-t0), row['Max Intensity'] + y_offset, str(p_idx), 
                     color=fg_col, ha='center', va='bottom', fontsize=8, fontweight='bold', 
                     clip_on=True, alpha=alpha, picker=5)
             self.spr_peak_label_map[t] = p_idx
        
        # Draw Peak Lines
        lines_10 = []
        lines_1 = []

        # VECTORIZED MARKER PLOTTING
        # Collect coordinates for vectorized plotting
        xs_10, ys_10 = [], []
        xs_1, ys_1 = [], []
        
        for i, row in self.spr_raw_results_df.iterrows():
            p_idx = int(row['Peak Index'])
            is_excluded = p_idx in excluded
            
            # FW0.1M (10%) Points (Always seconds X for plot alignment)
            h10, = ax.plot([(row['l01']-t0), (row['r01']-t0)], [row['h01'], row['h01']], 
                           linestyle='None', marker='s', markersize=5, alpha=0.3 if is_excluded else 0.9, 
                           color='C6', picker=5, label="FW 0.1M (10 %)" if i == 0 else "")
            lines_10.append(h10)
            self.spr_peak_label_map[h10] = p_idx
            
            # FW0.01M (1%) Points (Always seconds X for plot alignment)
            h1, = ax.plot([(row['l001']-t0), (row['r001']-t0)], [row['h001'], row['h001']], 
                          linestyle='None', marker='d', markersize=5, alpha=0.3 if is_excluded else 0.9, 
                          color='C9', picker=5, label="FW 0.01M (1 %)" if i == 0 else "")
            lines_1.append(h1)
            self.spr_peak_label_map[h1] = p_idx
            
        # Composite Peak Plotting
        if show_composite and self.spr_ax_comp:
            comp_df = Logic.generate_composite_peak(self.spr_df, iso, filtered_df, bg_sub=bg_sub, bg_start_s=bg_start_s, bg_end_s=bg_end_s)
            if comp_df is not None:
                # Use the anchor color (Index 0 of safe_palette)
                # Switch to ms for X-axis readability (User Request)
                self.spr_ax_comp.plot(comp_df['Relative Time (s)'] * 1000.0, comp_df['Normalised Intensity'], color=safe_palette[0], lw=2)
                self.spr_ax_comp.set_xlabel("Rel. Time (ms)", fontsize='medium', color=fg_col)
                # Stable composite title anchored to axis top
                self.spr_ax_comp.set_title("Composite Peak", fontsize='medium', color=fg_col, y=1.0, pad=10)
                self.spr_ax_comp.set_ylim(-0.02, 1.05) # Tighter Y limit
                self.spr_ax_comp.set_yticklabels([])    # Hide Y numbers
                self.spr_ax_comp.set_yticks([])         # Hide Y ticks
                self.spr_ax_comp.margins(x=0)         # Snap X axis to data
                self.spr_ax_comp.set_facecolor(plt.rcParams['axes.facecolor'])
                self.spr_ax_comp.tick_params(colors=fg_col)
                for spine in self.spr_ax_comp.spines.values():
                    spine.set_edgecolor(fg_col)
                self.spr_ax_comp.grid(True, alpha=0.2, color=fg_col)
                
                # Draw FW Markers for Composite
                comp_t = comp_df['Relative Time (s)'].values
                comp_y = comp_df['Normalised Intensity'].values
                
                comp_widths = {0.1: 0.0, 0.01: 0.0} # Store for stats box
                
                for lvl, mkr, col in [(0.1, 's', 'C6'), (0.01, 'd', 'C9')]:
                    above = np.where(comp_y >= lvl)[0]
                    if len(above) >= 2:
                        idx_l, idx_r = above[0], above[-1]
                        # Linear interpolation for left
                        if idx_l > 0:
                            t_l = np.interp(lvl, [comp_y[idx_l-1], comp_y[idx_l]], [comp_t[idx_l-1], comp_t[idx_l]])
                        else: t_l = comp_t[idx_l]
                        # Linear interpolation for right
                        if idx_r < len(comp_y) - 1:
                            t_r = np.interp(lvl, [comp_y[idx_r+1], comp_y[idx_r]], [comp_t[idx_r+1], comp_t[idx_r]])
                        else: t_r = comp_t[idx_r]
                        
                        # Store width
                        comp_widths[lvl] = t_r - t_l
                        
                        # Plot in ms
                        self.spr_ax_comp.plot([t_l * 1000.0, t_r * 1000.0], [lvl, lvl], linestyle='None', marker=mkr, markersize=5, color=col, alpha=0.8)

                # Summary Stats Box (Inlaid Legend)
                # Create persistent legend reference (Always create, control visibility)
                
                f10 = filtered_df['FW0.1M (s)']
                f1 = filtered_df['FW0.01M (s)']
                
                # Use calculated composite widths
                val_10 = comp_widths.get(0.1, 0.0)
                val_1 = comp_widths.get(0.01, 0.0)
                
                # Helper for adaptive rounding (Same as Avg/Max logic)
                def _adap_round_comp(v_ms):
                    if v_ms >= 100: return round(v_ms, 0)
                    if v_ms >= 10:  return round(v_ms, 1)
                    return round(v_ms, 2)
                
                # Store for Apply Buttons using exact same logic as Avg/Max (Always MS)
                self._last_spr_fw10_comp_ms = _adap_round_comp(val_10 * 1000.0)
                self._last_spr_fw1_comp_ms = _adap_round_comp(val_1 * 1000.0)
                
                # Update Labels in Results Panel (Bold Number, Normal Unit)
                if hasattr(self, 'lbl_spr_fw_comp10'):
                    self.lbl_spr_fw_comp10.setText(f"<b>{_adap(val_10*mult)}</b> {unit_label}")
                if hasattr(self, 'lbl_spr_fw_comp1'):
                    self.lbl_spr_fw_comp1.setText(f"<b>{_adap(val_1*mult)}</b> {unit_label}")
                
                label_10 = f"FW0.1M (10%)\nComp: {_adap(val_10*mult)} {unit_label}\nMax:  {_adap(f10.max()*mult)} {unit_label}"
                label_1 = f"FW0.01M (1%)\nComp: {_adap(val_1*mult)} {unit_label}\nMax:  {_adap(f1.max()*mult)} {unit_label}"
                
                h10 = Line2D([0], [0], linestyle='None', marker='s', markersize=4, color='C6', label=label_10)
                h1 = Line2D([0], [0], linestyle='None', marker='d', markersize=4, color='C9', label=label_1)
                
                self.spr_stats_legend = self.spr_ax_comp.legend(handles=[h10, h1], loc='upper right', 
                                                    fontsize=8.5, handletextpad=0.5, labelspacing=1.0,
                                                    frameon=True, facecolor=plt.rcParams['axes.facecolor'],
                                                    edgecolor=fg_col)
                self.spr_stats_legend.get_frame().set_alpha(0.6)
                for text in self.spr_stats_legend.get_texts():
                    text.set_color(fg_col)
                    
                # Set initial visibility based on checkbox
                self.spr_stats_legend.set_visible(not self.chk_hide_stats.isChecked())

        # Legend Frame Logic
        self.spr_legend_frame = None
        
        # Create Proxy Handles for Legend
        # (Must use proxies so legend doesn't inherit 'excluded' alpha from 1st point)
        legend_handles = [line]
        if lines_10:
             legend_handles.append(Line2D([0], [0], linestyle='None', marker='s', markersize=5, color='C6', alpha=0.9, label="FW 0.1M (10 %)"))
        if lines_1:
             legend_handles.append(Line2D([0], [0], linestyle='None', marker='d', markersize=5, color='C9', alpha=0.9, label="FW 0.01M (1 %)"))

        # Snap top margin to legend height exactly
        self.last_spr_rows = 1
        self._on_plot_resize(type('obj', (object,), {'canvas': self.spr_canvas}))
        
        ax.set_xlabel("Time (s)", fontsize='medium', color=fg_col)
        
        # Anchored Y-Axis Label (Fixed 0.3 inches from left edge)
        import matplotlib.transforms as mtransforms
        unit = getattr(self, 'cached_unit_label', "Counts")
        ax.set_ylabel(f"Intensity ({unit})", fontsize='medium', color=fg_col)
        # Use 0.3 inches to be safe from bezel
        w_in = self.spr_figure.get_size_inches()[0]
        ax.yaxis.set_label_coords(0.3/w_in, 0.5, 
                                  transform=mtransforms.blended_transform_factory(self.spr_figure.transFigure, ax.transAxes))
        
        # Connect dynamic margin update on zoom/pan
        ax.callbacks.connect('ylim_changed', lambda event: self._update_smart_margins(self.spr_canvas))
        
        # Engineering Notation (k, M, G)
        # Matches main Optimiser plot style
        ax.yaxis.set_major_formatter(EngFormatter(sep=" "))
        ax.yaxis.get_offset_text().set_color(fg_col)
        ax.grid(False, which='both')
        
        # Replace title with legend at top
        # Absolute top-anchored legend with fixed padding
        leg = ax.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, 1.0), 
                        borderaxespad=0.5, ncol=3, frameon=True, fontsize='medium', 
                        handlelength=1.5, handletextpad=0.7, columnspacing=1.5)
        leg.get_frame().set_alpha(0.0) # Transparent frame
        leg.get_frame().set_picker(5)
        self.spr_legend_frame = leg.get_frame()
        
        self.spr_map_legend_to_line = {}
        # Backward-compatible access to legend handles (legend_handles in Matplotlib >= 3.7)
        leg_handles = getattr(leg, 'legend_handles', getattr(leg, 'legendHandles', []))
        for legline, legtext in zip(leg_handles, leg.get_texts()):
            txt = legtext.get_text()
            legline.set_alpha(1.0) # Ensure fully visible initially
            legline.set_picker(5)
            legtext.set_picker(5)
            
            if txt == iso:
                # Main channel: cannot be hidden
                self.spr_map_legend_to_line[legline] = (line, False) # (Target, IsHidable)
                self.spr_map_legend_to_line[legtext] = (line, False)
            elif txt == "FW 0.1M (10 %)":
                self.spr_map_legend_to_line[legline] = (lines_10, True)
                self.spr_map_legend_to_line[legtext] = (lines_10, True)
            elif txt == "FW 0.01M (1 %)":
                self.spr_map_legend_to_line[legline] = (lines_1, True)
                self.spr_map_legend_to_line[legtext] = (lines_1, True)
        
        # Apply theme to spines and ticks
        ax.tick_params(colors=fg_col)
        for spine in ax.spines.values():
            spine.set_edgecolor(fg_col)

        # Consistent Grid
        # ax.grid(True, alpha=0.2, color=fg_col) # REMOVED per user request
        
        # Ensure Composite Plot also has no grid
        ax.grid(False, which='both')
        if self.spr_ax_comp:
            self.spr_ax_comp.minorticks_on()
            self.spr_ax_comp.grid(True, which='both', alpha=0.2, color=fg_col)

        # Force tight borders
        ax.margins(x=0)
        
        # Consistent layout management
        # (Removed constrained layout pad setting to preserve fixed subplots_adjust)
        
        # Consistent with fixed subplots_adjust, we do not call execute_constrained_layout


        if rescale is False and old_xlim is not None:
            ax.set_xlim(old_xlim)
            ax.set_ylim(old_ylim)
        elif rescale is True or (rescale is None and self._get(self.chk_rescale_spr, 'isChecked')):
            self.rescale_to_visible(ax=ax, canvas=None) # Defer draw

        # ALWAYS Rescale Composite Y (Autosize to 0-1 data + margin)
        # This fixes the issue where exclusion (rescale=False) left the composite plot unscaled
        if show_composite and self.spr_ax_comp:
             self.rescale_to_visible(ax=self.spr_ax_comp, rescale_y=True, canvas=None)

        # Force margin update to ensure axis doesn't collapse (especially if rescale=False)
        self._update_smart_margins(self.spr_canvas)

        self.spr_canvas.draw()

    def _auto_detect_spr_prominence(self, iso=None):
        if self.spr_df is None:
            return
            
        if iso is None:
            iso = self._get(self.cmb_spr_iso, 'currentText')
            
        if not iso or iso not in self.spr_df.columns:
            return
            
        y = self.spr_df[iso].values
        if len(y) < 2:
            return
            
        # Strategy: Use robust background noise estimation (similar to region auto-detect)
        try:
            # 1. Estimate baseline noise from the beginning (typically gas blank)
            init_window = max(10, int(len(y) * 0.05))
            bg_est = y[:init_window]
            bg_mean = np.nanmean(bg_est)
            bg_std = np.nanstd(bg_est) if len(bg_est) > 1 else 0
            
            y_max = np.nanmax(y)
            y_range = y_max - bg_mean
            
            # 2. Robust threshold: 10% of total range OR 10 sigma (to clear noise floor)
            # This prevents picking up noise spikes as valid peaks
            auto_val = max(y_range * 0.1, 10 * bg_std)
            
            auto_val = round(auto_val, 0)
                
            self._block_signals(self.spin_spr_prom, True)
            self.spin_spr_prom.setValue(auto_val)
            self._block_signals(self.spin_spr_prom, False)
            
            # Save the auto-detected baseline so it is respected as a default by the global table
            self.spr_channel_prominence[iso] = auto_val
        except Exception as e:
            # IoLog.error(f"iolite Optimiser: Error in auto-prominence detection: {e}")
            pass

    def _auto_detect_spr_bg(self, iso=None):
        if self.spr_df is None: return
        t_orig = self.spr_df['Time'].values
        if len(t_orig) == 0: return
        t0 = t_orig[0]
        
        try:
            (bg_s, bg_e), _ = Logic.auto_detect_regions(self.spr_df)
            bg_s_rel = max(0.0, float(bg_s - t0))
            bg_e_rel = max(0.0, float(bg_e - t0))
            self._block_signals(self.spin_spr_bg_start, True)
            self._block_signals(self.spin_spr_bg_end, True)
            self.spin_spr_bg_start.setValue(bg_s_rel)
            self.spin_spr_bg_end.setValue(bg_e_rel)
            self._block_signals(self.spin_spr_bg_start, False)
            self._block_signals(self.spin_spr_bg_end, False)
        except:
            pass

    def _on_spr_bg_changed(self):
        # 1. Turn off Auto if manually adjusted
        if self._get(self.chk_spr_auto_bg, 'isChecked'):
            self._block_signals(self.chk_spr_auto_bg, True)
            self.chk_spr_auto_bg.setChecked(False)
            self._block_signals(self.chk_spr_auto_bg, False)
        # 2. Debounce visual update
        self.spr_debounce_timer.start(300)

    def _on_spr_auto_bg_toggled(self, checked):
        if checked:
            self._auto_detect_spr_bg()
            self.run_spr_analysis(force_refresh=True)

    def _on_spr_iso_changed(self, text):
        if not text: return
        
        # Restore Auto Prominence state (default to True if not seen before)
        is_auto = True
        if hasattr(self, 'spr_channel_auto_prom') and text in self.spr_channel_auto_prom:
             is_auto = self.spr_channel_auto_prom[text]
             
        self._block_signals(self.chk_spr_auto_prom, True)
        self.chk_spr_auto_prom.setChecked(is_auto)
        self._block_signals(self.chk_spr_auto_prom, False)
        
        # Restore Auto-Exclude state
        self._block_signals(self.chk_spr_auto_exclude, True)
        self.chk_spr_auto_exclude.setChecked(getattr(self, 'spr_auto_exclude', False))
        self._block_signals(self.chk_spr_auto_exclude, False)
        
        # Restore Auto-Exclude Pct
        self._block_signals(self.spin_spr_auto_exclude_pct, True)
        self.spin_spr_auto_exclude_pct.setValue(getattr(self, 'spr_auto_exclude_pct', 90.0))
        self._block_signals(self.spin_spr_auto_exclude_pct, False)
        
        # Restore custom channel prominence if we're not in Auto
        if not is_auto and hasattr(self, 'spr_channel_prominence') and text in self.spr_channel_prominence:
             self._block_signals(self.spin_spr_prom, True)
             self.spin_spr_prom.setValue(self.spr_channel_prominence[text])
             self._block_signals(self.spin_spr_prom, False)
        elif is_auto:
             self._auto_detect_spr_prominence(text)
             
        self.run_spr_analysis()

    def _on_spr_prom_changed(self, *args):
        """
        Called when SPR peak cutoff (prominence) is manually adjusted.
        Auto-deselects the 'Auto' checkbox and applies 300ms debouncing.
        """
        # Save custom prominence for this channel
        curr_iso = self._get(self.cmb_spr_iso, 'currentText')
        if curr_iso and hasattr(self, 'spr_channel_prominence'):
             self.spr_channel_prominence[curr_iso] = self._get(self.spin_spr_prom, 'value')
             self.spr_channel_auto_prom[curr_iso] = False
             
        # 1. Uncheck the Auto box if it's currently checked
        if self._get(self.chk_spr_auto_prom, 'isChecked'):
            self._block_signals(self.chk_spr_auto_prom, True)
            self.chk_spr_auto_prom.setChecked(False)
            self._block_signals(self.chk_spr_auto_prom, False)
            # Ensure the spinbox remains enabled (it might have been disabled while Auto was on)
            self.spin_spr_prom.setEnabled(True)
        
        # 2. Restart/Start the 500ms debounce timer
        # This prevents the plot from flickering/re-calculating on every digit entered
        self.spr_debounce_timer.start(500)

    def _on_spr_dist_changed(self, *args):
        """
        Called when SPR Min Distance is manually adjusted.
        Applies 500ms debouncing.
        """
        # Start/Restart the 500ms debounce timer
        self.spr_debounce_timer.start(500)

    def _on_spr_auto_prom_toggled(self, checked):
        # Removed: self.spin_spr_prom.setEnabled(not checked)
        curr_iso = self._get(self.cmb_spr_iso, 'currentText')
        if curr_iso and hasattr(self, 'spr_channel_auto_prom'):
             self.spr_channel_auto_prom[curr_iso] = checked
             
        if checked:
            self.run_spr_analysis(force_refresh=True)

    def _on_spr_auto_exclude_toggled(self, checked):
        self.spr_auto_exclude = checked
        self.run_spr_analysis(force_refresh=True)
        
    def _on_spr_auto_exclude_pct_changed(self, *args):
        self.spr_auto_exclude_pct = self._get(self.spin_spr_auto_exclude_pct, 'value')
        
        # Wait for the user to finish typing 
        self.spr_debounce_timer.start(500)

    def _on_stats_toggled(self, checked):
        # Toggle visibility without re-running analysis (Preserves Zoom)
        if hasattr(self, 'spr_stats_legend') and self.spr_stats_legend:
            self.spr_stats_legend.set_visible(not checked)
            # Efficient redraw if available
            if hasattr(self.spr_canvas, 'draw_idle'):
                self.spr_canvas.draw_idle()
            else:
                self.spr_canvas.draw()

    def _update_spr_time_suffixes(self, *args):
        unit_text = self._get(self.cmb_spr_unit, 'currentText')
        if unit_text is None: return
        is_ms = "ms" in unit_text
        suffix = " ms" if is_ms else " s"
        decimals = 2 if is_ms else 4
        
        self.spin_spr_dist.setSuffix(suffix)
        self.spin_spr_dist.setDecimals(decimals)
        
        self.spin_spr_baseline_window.setSuffix(f'{suffix} (0=Inf)')
        self.spin_spr_baseline_window.setDecimals(decimals)
        
        self.spin_spr_smooth_window.setSuffix(suffix)
        self.spin_spr_smooth_window.setDecimals(decimals)

    def _on_spr_apply_clicked(self, mode):
        # mode can be '10avg', '10max', '10comp', '1avg', '1max', '1comp'
        mapping = {
            '10avg': getattr(self, '_last_spr_fw10_avg_ms', 0),
            '10max': getattr(self, '_last_spr_fw10_max_ms', 0),
            '10comp': getattr(self, '_last_spr_fw10_comp_ms', 0),
            '1avg': getattr(self, '_last_spr_fw1_avg_ms', 0),
            '1max': getattr(self, '_last_spr_fw1_max_ms', 0),
            '1comp': getattr(self, '_last_spr_fw1_comp_ms', 0)
        }
        val_ms = mapping.get(mode, 0)
        if '10' in mode:
            self.active_washout_level = 0.1
        else:
            self.active_washout_level = 0.01
        self.spin_wash.setValue(val_ms)
        self.tabs.setCurrentIndex(1) # Switch back to Optimiser
        IoLog.information(f"iolite Optimiser: Applied {mode} washout value: {val_ms:.2f} ms")

    def _on_spr_show_all_stats_toggled(self, checked):
        if checked:
            self.spr_table_stack.setCurrentIndex(1)
            self.lbl_spr_table_title.setText("All Channel Stats")
        else:
            self.spr_table_stack.setCurrentIndex(0)
            num_detected = len(getattr(self, 'spr_raw_results_df', []))
            num_excluded = len(getattr(self, 'spr_excluded_peaks', {}).get(self._get(self.cmb_spr_iso, 'currentText'), []))
            self.lbl_spr_table_title.setText(f"SPR Peaks Detected - {num_detected} | SPR Peaks Excluded - {num_excluded}")
            

    def _update_all_channel_stats(self):
        if getattr(self, 'spr_df', None) is None or find_peaks is None:
            return

        try:
            # 1. Get current detection settings
            ui_prominence = self._get(self.spin_spr_prom, 'value')
            distance = self._get(self.spin_spr_dist, 'value')
            wlen_val = self._get(self.spin_spr_baseline_window, 'value')
            wlen = None if wlen_val == 0 else wlen_val
            apply_smoothing = self._get(self.chk_spr_smooth, 'isChecked')
            smooth_window = self._get(self.spin_spr_smooth_window, 'value')
            unit_text = self._get(self.cmb_spr_unit, 'currentText')
            is_ms = "ms" in unit_text
            mult = 1000.0 if is_ms else 1.0
            unit_label = "ms" if is_ms else "s"
            
            bg_sub = self._get(self.chk_spr_bg_sub, 'isChecked')
            bg_start_s = self._get(self.spin_spr_bg_start, 'value')
            bg_end_s = self._get(self.spin_spr_bg_end, 'value')

            winning_iso = None
            max_fw1 = -1.0
            winning_stats = None
            
            table_rows = []

            # 2. Iterate through all channels in the combo boxes
            all_isos = [self.cmb_spr_iso.itemText(i) for i in range(self.cmb_spr_iso.count)]
            current_ui_iso = self._get(self.cmb_spr_iso, 'currentText')
            
            for iso in all_isos:
                if not iso or iso not in self.spr_df.columns:
                    continue
                    
                # Calculate optimal prominence for this channel
                is_auto_prom = getattr(self, 'spr_channel_auto_prom', {}).get(iso, True)
                
                if iso == current_ui_iso:
                    # Use the UI overridden value for the currently viewed channel
                    prominence = ui_prominence
                elif not is_auto_prom and hasattr(self, 'spr_channel_prominence') and iso in self.spr_channel_prominence:
                    # Use the stored manual value for other channels
                    prominence = self.spr_channel_prominence[iso]
                else:
                    y = self.spr_df[iso].values
                    if len(y) < 2:
                        prominence = 100
                    else:
                        init_window = max(10, int(len(y) * 0.05))
                        bg_est = y[:init_window]
                        bg_mean = np.nanmean(bg_est)
                        bg_std = np.nanstd(bg_est) if len(bg_est) > 1 else 0
                        y_max = np.nanmax(y)
                        y_range = y_max - bg_mean
                        prominence = max(y_range * 0.1, 10 * bg_std)
                        prominence = round(prominence, 0)
                        if hasattr(self, 'spr_channel_prominence'):
                            self.spr_channel_prominence[iso] = prominence

                # Detect CPS status for this channel
                is_cps = True
                if hasattr(self, 'channel_is_cps') and self.channel_is_cps.get(iso, False):
                    is_cps = True
                elif getattr(self, 'cached_unit_label', "Counts") == "Counts":
                    is_cps = False

                distance_s = distance / mult
                wlen_s = wlen / mult if wlen is not None else None
                smooth_window_s = smooth_window / mult

                raw_result = Logic.analyse_washout_peaks(
                    self.spr_df, iso, prominence, min_distance=distance_s, is_cps=is_cps,
                    wlen=wlen_s, apply_smoothing=apply_smoothing, smooth_window=smooth_window_s,
                    bg_sub=bg_sub, bg_start_s=bg_start_s, bg_end_s=bg_end_s
                )
                
                if isinstance(raw_result, tuple):
                    table_rows.append((iso, {'Count': 0}))
                    if hasattr(self, 'spr_all_raw_results'):
                        self.spr_all_raw_results[iso] = raw_result
                    continue
                
                if hasattr(self, 'spr_all_raw_results'):
                     self.spr_all_raw_results[iso] = raw_result
                     
                df = raw_result
                excluded = getattr(self, 'spr_excluded_peaks', {}).get(iso, set()).copy()
                
                # Auto-exclusion filter
                auto_excl = getattr(self, 'spr_auto_exclude', False)
                if auto_excl and len(df) > 5:
                    pct_val = getattr(self, 'spr_auto_exclude_pct', 90.0)
                    fw = df['r001'] - df['l001']
                    fw_p_upper = np.percentile(fw, pct_val)
                    
                    outliers = df[fw > fw_p_upper]['Peak Index'].values
                    excluded.update(outliers)
                
                if excluded:
                     df = df[~df['Peak Index'].isin(excluded)]
                
                stats = Logic.summarize_peaks(df)
                
                if "Error" in stats or stats.get("Count", 0) == 0:
                    table_rows.append((iso, {'Count': 0}))
                    continue

                comp_df = Logic.generate_composite_peak(self.spr_df, iso, df, bg_sub=bg_sub, bg_start_s=bg_start_s, bg_end_s=bg_end_s)
                comp_width_1 = 0.0
                comp_width_10 = 0.0
                if comp_df is not None:
                    comp_t = comp_df['Relative Time (s)'].values
                    comp_y = comp_df['Normalised Intensity'].values
                    for lvl in [0.1, 0.01]:
                        above = np.where(comp_y >= lvl)[0]
                        if len(above) >= 2:
                            idx_l, idx_r = above[0], above[-1]
                            t_l = np.interp(lvl, [comp_y[idx_l-1], comp_y[idx_l]], [comp_t[idx_l-1], comp_t[idx_l]]) if idx_l > 0 else comp_t[idx_l]
                            t_r = np.interp(lvl, [comp_y[idx_r+1], comp_y[idx_r]], [comp_t[idx_r+1], comp_t[idx_r]]) if idx_r < len(comp_y)-1 else comp_t[idx_r]
                            if lvl == 0.01: comp_width_1 = t_r - t_l
                            else: comp_width_10 = t_r - t_l

                stats['FW0.01M Comp'] = comp_width_1
                stats['FW0.1M Comp'] = comp_width_10
                
                table_rows.append((iso, stats))

                current_fw1 = stats['FW0.01M Max']
                if current_fw1 > max_fw1:
                    max_fw1 = current_fw1
                    winning_iso = iso
                    winning_stats = stats

            # 3. Update UI
            
            def _adap(v):
                if v >= 100: return f"{v:.0f}"
                if v >= 10:  return f"{v:.1f}"
                if v >= 1:   return f"{v:.2f}"
                return f"{v:.3f}"
                
            if winning_iso:
                self.lbl_spr_max_iso.setText(f"Channel: <b>{winning_iso}</b>")
                
                def _fmt(val_s):
                    return f"<b>{_adap(val_s * mult)}</b> {unit_label}"

                self.lbl_spr_max_avg10.setText(_fmt(winning_stats['FW0.1M Mean']))
                self.lbl_spr_max_avg1.setText(_fmt(winning_stats['FW0.01M Mean']))
                self.lbl_spr_max_max10.setText(_fmt(winning_stats['FW0.1M Max']))
                self.lbl_spr_max_max1.setText(_fmt(winning_stats['FW0.01M Max']))
                self.lbl_spr_max_comp10.setText(_fmt(winning_stats['FW0.1M Comp']))
                self.lbl_spr_max_comp1.setText(_fmt(winning_stats['FW0.01M Comp']))
                
                def _rnd(v_s):
                    v_ms = v_s * 1000.0
                    if v_ms >= 100: return round(v_ms, 0)
                    if v_ms >= 10:  return round(v_ms, 1)
                    return round(v_ms, 2)

                self._max_spr_winners = {
                    'avg': _rnd(winning_stats['FW0.01M Mean']),
                    'max': _rnd(winning_stats['FW0.01M Max']),
                    'comp': _rnd(winning_stats['FW0.01M Comp'])
                }
            else:
                self.lbl_spr_max_iso.setText("Channel: None Detected")
                self.lbl_spr_max_avg10.setText("-")
                self.lbl_spr_max_avg1.setText("-")
                self.lbl_spr_max_max10.setText("-")
                self.lbl_spr_max_max1.setText("-")
                self.lbl_spr_max_comp10.setText("-")
                self.lbl_spr_max_comp1.setText("-")
                if hasattr(self, '_max_spr_winners'):
                    del self._max_spr_winners
                    
            # 4. Populate spr_all_stats_table
            def _it(text):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                item.setTextAlignment(Qt.AlignCenter)
                return item
                
            def _si(val):
                if val >= 1e9: return f"{val/1e9:.2f}", " G"
                if val >= 1e6: return f"{val/1e6:.2f}", " M"
                if val >= 1e3: return f"{val/1e3:.2f}", " k"
                return f"{val:.2f}", " "

            self.spr_all_stats_table.setRowCount(0)
            self.spr_all_stats_table.setRowCount(len(table_rows))
            
            # ["Channel", "Avg (10%)", "Max (10%)", "Comp (10%)", "RSD (10%)", "Area 10% (Avg)", "Avg (1%)", "Max (1%)", "Comp (1%)", "RSD (1%)", "Area 1% (Avg)"]
            for i, (iso, stats) in enumerate(table_rows):
                self.spr_all_stats_table.setItem(i, 0, _it(iso))
                if stats.get('Count', 0) == 0:
                    for col in range(1, 11):
                        self.spr_all_stats_table.setItem(i, col, _it("-"))
                else:
                    self.spr_all_stats_table.setItem(i, 1, _it(f"{_adap(stats.get('FW0.1M Mean', 0) * mult)}"))
                    self.spr_all_stats_table.setItem(i, 2, _it(f"{_adap(stats.get('FW0.1M Max', 0) * mult)}"))
                    self.spr_all_stats_table.setItem(i, 3, _it(f"{_adap(stats.get('FW0.1M Comp', 0) * mult)}"))
                    self.spr_all_stats_table.setItem(i, 4, _it(f"{stats.get('FW0.1M RSD', 0):.1f} %"))
                    s_a10, p_a10 = _si(stats.get('Area 10% Mean', 0.0))
                    self.spr_all_stats_table.setItem(i, 5, _it(f"{s_a10}{p_a10}"))
                    
                    self.spr_all_stats_table.setItem(i, 6, _it(f"{_adap(stats.get('FW0.01M Mean', 0) * mult)}"))
                    self.spr_all_stats_table.setItem(i, 7, _it(f"{_adap(stats.get('FW0.01M Max', 0) * mult)}"))
                    self.spr_all_stats_table.setItem(i, 8, _it(f"{_adap(stats.get('FW0.01M Comp', 0) * mult)}"))
                    self.spr_all_stats_table.setItem(i, 9, _it(f"{stats.get('FW0.01M RSD', 0):.1f} %"))
                    s_a001, p_a001 = _si(stats.get('Area 1% Mean', 0.0))
                    self.spr_all_stats_table.setItem(i, 10, _it(f"{s_a001}{p_a001}"))

        except Exception as e:
            # IoLog.error(f"iolite Optimiser: Error in All Channel Stats update: {e}")
            pass

    def _on_spr_apply_max_clicked(self, mode):
        # mode can be 'avg', 'max', 'comp' (applying the 1% values of the winning channel)
        if not hasattr(self, '_max_spr_winners'):
            QMessageBox.warning(self, "No Results", "Please run 'All Channel Stats' first.")
            return
            
        val_ms = self._max_spr_winners.get(mode, 0)
        self.active_washout_level = 0.01
        self.spin_wash.setValue(val_ms)
        self.tabs.setCurrentIndex(1) # Switch back to Optimiser
        IoLog.information(f"iolite Optimiser: Applied Max Channel ({mode}) washout value: {val_ms:.2f} ms")

    def init_optimiser_tab(self):
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()
        
        # 0. Data Selection & Settings (Fixed at Top)
        h_ctrl = QHBoxLayout()
        self.btn_run = QPushButton("Import Optimisation File")
        self.btn_run.setFixedHeight(30)
        self.btn_run.clicked.connect(lambda: self.import_and_unload_data(tab="opt"))
        
        self.btn_settings = QPushButton("Settings")
        self.btn_settings.setFixedHeight(30)
        self.btn_settings.clicked.connect(self.show_settings_dialog)
        
        h_ctrl.addWidget(self.btn_run)
        h_ctrl.addWidget(self.btn_settings)
        left_layout.addLayout(h_ctrl)

        # --- LEFT COLUMN (SETTINGS - SCROLLABLE) ---
        self.settings_scroll_area = QScrollArea()
        self.settings_scroll_area.setWidgetResizable(True)
        self.settings_scroll_area.setFrameShape(QScrollArea.NoFrame)
        self.settings_scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        
        scroll_content = QWidget()
        l_settings = QVBoxLayout(scroll_content)
        l_settings.setContentsMargins(5, 2, 5, 2)
        l_settings.setSpacing(1) # Tightened spacing

        # 1. Hardware Configuration Summary
        self.grp_hw_summary = QGroupBox("")
        v_hw_sum = QVBoxLayout()
        v_hw_sum.setSpacing(0)
        v_hw_sum.setContentsMargins(2, 2, 2, 2)
        
        lbl_hw_title = QLabel("Hardware Configuration")
        lbl_hw_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_hw_title.setAlignment(Qt.AlignCenter)
        v_hw_sum.addWidget(lbl_hw_title)
        
        h_hw_cols = QHBoxLayout()
        h_hw_cols.setContentsMargins(2, 2, 2, 2)
        h_hw_cols.setSpacing(15)
        
        self.lbl_icp_sum = QLabel("Initializing ICP...")
        self.lbl_icp_sum.setWordWrap(True)
        self.lbl_icp_sum.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        
        self.lbl_las_sum = QLabel("Initializing Laser...")
        self.lbl_las_sum.setWordWrap(True)
        self.lbl_las_sum.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        
        h_hw_cols.addWidget(self.lbl_las_sum, 1, Qt.AlignHCenter)
        h_hw_cols.addWidget(self.lbl_icp_sum, 1, Qt.AlignHCenter)
        
        v_hw_sum.addLayout(h_hw_cols)
        self.grp_hw_summary.setLayout(v_hw_sum)
        l_settings.addWidget(self.grp_hw_summary)

        # 1. Instrument selection (Unified Form Layout for Alignment)
        grp_hw = QGroupBox("")
        v_hw = QVBoxLayout()
        v_hw.setSpacing(0)
        v_hw.setContentsMargins(5, 2, 5, 2)
        
        lbl_h_title = QLabel("ICP-MS Hardware")
        lbl_h_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_h_title.setAlignment(Qt.AlignCenter)
        v_hw.addWidget(lbl_h_title)
        
        l_hw = QVBoxLayout()
        
        self.form_icp = QFormLayout()
        self.form_icp.setRowWrapPolicy(QFormLayout.DontWrapRows)
        self.form_icp.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.form_icp.setLabelAlignment(Qt.AlignLeft)
        self.form_icp.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.form_icp.setHorizontalSpacing(20)
        self.form_icp.setVerticalSpacing(2)
        self.form_icp.setContentsMargins(0, 0, 0, 0)
        
        lbl_w = 160
        self.lbl_mfr = QLabel("Manufacturer:"); self.lbl_mfr.setFixedWidth(lbl_w)
        self.cmb_mfr = QComboBox()
        self.cmb_mfr.setFixedWidth(250)
        self.cmb_mfr.addItems(sorted(ICP_SPECS.keys()))
        saved_mfr = self.persistent_settings.get('mfr')
        if saved_mfr in ICP_SPECS:
            self.cmb_mfr.setCurrentText(saved_mfr)
        self.form_icp.addRow(self.lbl_mfr, self.cmb_mfr)
        
        self.lbl_model = QLabel("Model:"); self.lbl_model.setFixedWidth(lbl_w)
        self.cmb_model = QComboBox()
        self.cmb_model.setFixedWidth(250)
        self.form_icp.addRow(self.lbl_model, self.cmb_model)

        self.lbl_dwell_unit = QLabel("Dwell Times in:"); self.lbl_dwell_unit.setFixedWidth(lbl_w)
        self.cmb_dwell_unit = QComboBox()
        self.cmb_dwell_unit.setFixedWidth(250)
        self.cmb_dwell_unit.addItems(["Auto", "ms", "s"])
        saved_unit = self.persistent_settings.get('dwell_unit_pref', "Auto")
        self.cmb_dwell_unit.setCurrentText(saved_unit)
        self.cmb_dwell_unit.currentTextChanged.connect(lambda: self.run_optimisation(refresh=False))
        self.form_icp.addRow(self.lbl_dwell_unit, self.cmb_dwell_unit)

        # ICP Status (Always Visible)
        self.lbl_icp_status = QLabel("-")
        self.lbl_icp_status.setWordWrap(True)
        self.lbl_icp_status.setStyleSheet("font-size: 11px; font-weight: bold; font-style: italic;")
        self.form_icp.addRow(self.lbl_icp_status)
        
        # --- Custom ICP Rows ---
        self.lbl_cust_type = QLabel("Type:"); self.lbl_cust_type.setFixedWidth(lbl_w)
        self.cmb_cust_type = QComboBox()
        self.cmb_cust_type.setFixedWidth(250)
        self.cmb_cust_type.addItems(["Quadrupole", "Sector-Field", "Multi-Collector", "TOF"])
        self.cmb_cust_type.setCurrentText(self.persistent_settings.get('cust_type', "Quadrupole"))
        self.form_icp.addRow(self.lbl_cust_type, self.cmb_cust_type)
        
        self.lbl_cust_min = QLabel("Minimum Dwell Time (ms):"); self.lbl_cust_min.setFixedWidth(lbl_w)
        self.spin_cust_min = QDoubleSpinBox()
        self.spin_cust_min.setRange(0, 100); self.spin_cust_min.setFixedWidth(100); self.spin_cust_min.setDecimals(4)
        self.spin_cust_min.setValue(self.persistent_settings.get('cust_min', 0.1))
        self.form_icp.addRow(self.lbl_cust_min, self.spin_cust_min)
        
        self.lbl_cust_prec = QLabel("Dwell Resolution (ms):"); self.lbl_cust_prec.setFixedWidth(lbl_w)
        self.spin_cust_prec = QDoubleSpinBox()
        self.spin_cust_prec.setRange(0.0001, 100); self.spin_cust_prec.setFixedWidth(100); self.spin_cust_prec.setDecimals(4)
        self.spin_cust_prec.setValue(self.persistent_settings.get('cust_prec', 0.1))
        self.form_icp.addRow(self.lbl_cust_prec, self.spin_cust_prec)

        self.lbl_cust_dwells = QLabel("Allowed Dwell times (ms):"); self.lbl_cust_dwells.setFixedWidth(lbl_w)
        self.edit_cust_dwells = QLineEdit()
        self.edit_cust_dwells.setFixedWidth(250)
        self.edit_cust_dwells.setPlaceholderText("0.1, 0.2, 0.5")
        self.edit_cust_dwells.setToolTip("Comma separated list of valid dwells (ms)")
        self.edit_cust_dwells.setText(self.persistent_settings.get('cust_dwell_list', ""))
        self.form_icp.addRow(self.lbl_cust_dwells, self.edit_cust_dwells)
        
        # --- TOF Specific Rows ---
        self.lbl_tof_margin = QLabel("Washout Margin (%):"); self.lbl_tof_margin.setFixedWidth(lbl_w)
        self.spin_tof_margin = QDoubleSpinBox()
        self.spin_tof_margin.setRange(0, 1000); self.spin_tof_margin.setFixedWidth(100); self.spin_tof_margin.setDecimals(1)
        self.spin_tof_margin.setValue(self.persistent_settings.get('tof_margin', 10.0))
        self.spin_tof_margin.setToolTip("<b>Signal Capture Buffer</b><br/>A percentage buffer added to the system's baseline washout time. For TOF systems, this ensures sufficient time to measure the complete Single Pulse Response (SPR) peak, avoiding truncation due to signal jitter and extraction timings.")
        self.spin_tof_margin.valueChanged.connect(self._on_ui_change)
        self.form_icp.addRow(self.lbl_tof_margin, self.spin_tof_margin)

        l_hw.setContentsMargins(10, 2, 10, 2)
        l_hw.setSpacing(2)
        l_hw.addLayout(self.form_icp)
        
        v_hw.addLayout(l_hw)
        grp_hw.setLayout(v_hw)
        grp_hw.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        # Move to Dialog
        self.settings_dlg.main_layout.insertWidget(0, grp_hw)

        # 1.5 Laser selection (Unified Form Layout for Alignment)
        grp_laser = QGroupBox("")
        v_laser_main = QVBoxLayout()
        v_laser_main.setSpacing(0)
        v_laser_main.setContentsMargins(5, 2, 5, 2)
        
        lbl_l_title = QLabel("Laser Hardware")
        lbl_l_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_l_title.setAlignment(Qt.AlignCenter)
        v_laser_main.addWidget(lbl_l_title)
        
        l_laser = QVBoxLayout()
        
        self.form_las = QFormLayout()
        self.form_las.setRowWrapPolicy(QFormLayout.DontWrapRows)
        self.form_las.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self.form_las.setLabelAlignment(Qt.AlignLeft)
        self.form_las.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.form_las.setHorizontalSpacing(20)
        self.form_las.setVerticalSpacing(2)
        self.form_las.setContentsMargins(0, 0, 0, 0)
        
        self.lbl_plat = QLabel("Platform:"); self.lbl_plat.setFixedWidth(lbl_w)
        self.cmb_laser_mfr = QComboBox()
        self.cmb_laser_mfr.setFixedWidth(250)
        self.cmb_laser_mfr.addItems(sorted(LASER_PLATFORMS.keys()))
        saved_lmfr = self.persistent_settings.get('laser_mfr')
        if saved_lmfr in LASER_PLATFORMS:
            self.cmb_laser_mfr.setCurrentText(saved_lmfr)
        self.form_las.addRow(self.lbl_plat, self.cmb_laser_mfr)
        
        self.lbl_source = QLabel("Laser Source:"); self.lbl_source.setFixedWidth(lbl_w)
        self.cmb_laser_mod = QComboBox()
        self.cmb_laser_mod.setFixedWidth(250)
        self.form_las.addRow(self.lbl_source, self.cmb_laser_mod)

        self.lbl_cell = QLabel("Cell Type:"); self.lbl_cell.setFixedWidth(lbl_w)
        self.cmb_cell = QComboBox()
        self.cmb_cell.setFixedWidth(250)
        self.form_las.addRow(self.lbl_cell, self.cmb_cell)
        
        # Laser Status (Always Visible)
        self.lbl_laser_status = QLabel("-")
        self.lbl_laser_status.setWordWrap(True)
        self.lbl_laser_status.setStyleSheet("font-size: 11px; font-weight: bold; font-style: italic;")
        self.form_las.addRow(self.lbl_laser_status)
        
        # --- Custom Laser Rows ---
        self.lbl_cust_mode = QLabel("Mode:"); self.lbl_cust_mode.setFixedWidth(lbl_w)
        self.cmb_cust_rr_type = QComboBox()
        self.cmb_cust_rr_type.setFixedWidth(250)
        self.cmb_cust_rr_type.addItems(["Maximum Rep-Rate", "Discrete"])
        self.cmb_cust_rr_type.setCurrentText(self.persistent_settings.get('cust_rr_type', 'Maximum Rep-Rate'))
        self.form_las.addRow(self.lbl_cust_mode, self.cmb_cust_rr_type)
        
        self.lbl_cust_rr_max = QLabel("Max Rep-Rate (Hz):"); self.lbl_cust_rr_max.setFixedWidth(lbl_w)
        self.spin_cust_rr = QSpinBox()
        self.spin_cust_rr.setRange(1, 10000); self.spin_cust_rr.setFixedWidth(100)
        self.spin_cust_rr.setValue(self.persistent_settings.get('cust_rr', 100))
        self.form_las.addRow(self.lbl_cust_rr_max, self.spin_cust_rr)

        self.lbl_cust_rr_prec = QLabel("Rep-Rate Resolution (Hz):"); self.lbl_cust_rr_prec.setFixedWidth(lbl_w)
        self.spin_cust_rr_prec = QDoubleSpinBox()
        self.spin_cust_rr_prec.setRange(0.001, 100); self.spin_cust_rr_prec.setFixedWidth(100); self.spin_cust_rr_prec.setDecimals(3)
        self.spin_cust_rr_prec.setValue(self.persistent_settings.get('cust_rr_prec', 1.0))
        self.spin_cust_rr_prec.setToolTip("Rep Rate Resolution (Hz)")
        self.form_las.addRow(self.lbl_cust_rr_prec, self.spin_cust_rr_prec)

        self.lbl_cust_rr_list = QLabel("Allowed Rep-Rates (Hz):"); self.lbl_cust_rr_list.setFixedWidth(lbl_w)
        self.edit_cust_rr_list = QLineEdit()
        self.edit_cust_rr_list.setFixedWidth(250)
        self.edit_cust_rr_list.setPlaceholderText("1, 2, 5, 10, 20")
        self.edit_cust_rr_list.setToolTip("Comma separated list of valid rep rates")
        self.edit_cust_rr_list.setText(self.persistent_settings.get('cust_rr_list', ""))
        self.form_las.addRow(self.lbl_cust_rr_list, self.edit_cust_rr_list)

        self.lbl_cust_speed = QLabel("Max Speed (µm s⁻¹):"); self.lbl_cust_speed.setFixedWidth(lbl_w)
        self.spin_cust_speed = QSpinBox()
        self.spin_cust_speed.setRange(1, 100000); self.spin_cust_speed.setFixedWidth(90)
        self.spin_cust_speed.setValue(self.persistent_settings.get('cust_speed', 2000))
        self.form_las.addRow(self.lbl_cust_speed, self.spin_cust_speed)
        
        l_laser.setContentsMargins(10, 2, 10, 2)
        l_laser.setSpacing(2)
        l_laser.addLayout(self.form_las)
        

        v_laser_main.addLayout(l_laser)
        grp_laser.setLayout(v_laser_main)
        grp_laser.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        # Move to Dialog
        self.settings_dlg.main_layout.insertWidget(1, grp_laser)

        # 1.6 Plot Options (New Section)
        grp_plot = QGroupBox("")
        v_plot = QVBoxLayout()
        v_plot.setSpacing(0)
        v_plot.setContentsMargins(5, 2, 5, 2)
        
        lbl_p_title = QLabel("Plot Options")
        lbl_p_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_p_title.setAlignment(Qt.AlignCenter)
        v_plot.addWidget(lbl_p_title)
        
        l_plot = QHBoxLayout()
        l_plot.setContentsMargins(10, 5, 10, 5)
        
        self.chk_limit_vis = QCheckBox("Limit Visible Plots:")
        self.chk_limit_vis.setChecked(self.persistent_settings.get('limit_vis', False))
        self.chk_limit_vis.toggled.connect(lambda: self.update_plot(preserve_zoom=True))
        
        self.spin_max_vis = QSpinBox()
        self.spin_max_vis.setRange(1, 100)
        self.spin_max_vis.setValue(int(self.persistent_settings.get('max_vis', 5)))
        self.spin_max_vis.setFixedWidth(50)
        self.spin_max_vis.valueChanged.connect(lambda: self.update_plot(preserve_zoom=True))
        
        l_plot.addWidget(self.chk_limit_vis)
        l_plot.addWidget(self.spin_max_vis)
        l_plot.addStretch()
        
        v_plot.addLayout(l_plot)
        grp_plot.setLayout(v_plot)
        grp_plot.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        
        self.settings_dlg.main_layout.insertWidget(2, grp_plot)

        # 1.6 Connect signals for immediate update
        # (Signals will be connected centrally in initUI to avoid duplication)
        # Add stretch to bottom to force the group boxes to collapse vertically
        self.settings_dlg.main_layout.insertStretch(2)


        # Input Parameters
        grp_input = QGroupBox("")
        v_input = QVBoxLayout()
        v_input.setSpacing(0)
        v_input.setContentsMargins(2, 2, 2, 2)
        
        lbl_in_title = QLabel("Input Parameters")
        lbl_in_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_in_title.setAlignment(Qt.AlignCenter)
        v_input.addWidget(lbl_in_title)
        
        grid_input = QGridLayout()
        grid_input.setContentsMargins(2, 2, 2, 2)
        grid_input.setSpacing(5)

        lbl_w_in = 125
        spin_w_in = 70

        # Row 0
        lbl_spot = QLabel("Initial Spot Size (µm):"); lbl_spot.setFixedWidth(lbl_w_in)
        grid_input.addWidget(lbl_spot, 0, 0)
        self.spin_spot = QSpinBox()
        self.spin_spot.setRange(1, 300); self.spin_spot.setValue(int(self.persistent_settings.get('spot_size', 50)))
        self.spin_spot.setFixedWidth(spin_w_in)
        grid_input.addWidget(self.spin_spot, 0, 1)

        lbl_wash = QLabel("Washout (ms):"); lbl_wash.setFixedWidth(lbl_w_in)
        grid_input.addWidget(lbl_wash, 0, 3)
        self.spin_wash = QDoubleSpinBox()
        self.spin_wash.setRange(0.01, 50000); self.spin_wash.setDecimals(2)
        self.spin_wash.setValue(self.persistent_settings.get('washout', 500))
        self.spin_wash.setFixedWidth(spin_w_in)
        grid_input.addWidget(self.spin_wash, 0, 4)

        # Row 1
        lbl_init_rr = QLabel("Initial Rep-Rate (Hz):"); lbl_init_rr.setFixedWidth(lbl_w_in)
        grid_input.addWidget(lbl_init_rr, 1, 0)
        self.spin_init_rr = QDoubleSpinBox()
        self.spin_init_rr.setRange(1, 10000); self.spin_init_rr.setValue(self.persistent_settings.get('init_rr', 20))
        self.spin_init_rr.setFixedWidth(spin_w_in)
        grid_input.addWidget(self.spin_init_rr, 1, 1)

        lbl_mode = QLabel("Mode:"); lbl_mode.setFixedWidth(lbl_w_in)
        grid_input.addWidget(lbl_mode, 1, 3)
        self.cmb_mode = QComboBox()
        self.cmb_mode.addItems(["Spot", "Line", "Imaging"])
        self.cmb_mode.setCurrentText(self.persistent_settings.get('mode', 'Spot'))
        self.cmb_mode.setFixedWidth(spin_w_in)
        grid_input.addWidget(self.cmb_mode, 1, 4)

        # Row 2
        lbl_num_analytes = QLabel("Number of Analytes:"); lbl_num_analytes.setFixedWidth(lbl_w_in)
        grid_input.addWidget(lbl_num_analytes, 2, 0)
        self.spin_num_analytes = QSpinBox()
        self.spin_num_analytes.setRange(1, 100)
        self.spin_num_analytes.setValue(int(self.persistent_settings.get('num_analytes', 10)))
        self.spin_num_analytes.setFixedWidth(spin_w_in)
        grid_input.addWidget(self.spin_num_analytes, 2, 1)

        self.btn_suggest_prescan = QPushButton("Suggest PreScan Params")
        self.btn_suggest_prescan.clicked.connect(self.show_prescan_suggestions)
        grid_input.addWidget(self.btn_suggest_prescan, 2, 3, 1, 2)

        # Let Column 2 act as the expanding center space between Column 1 and Column 3
        grid_input.setColumnStretch(2, 1)

        v_input.addLayout(grid_input)
        grp_input.setLayout(v_input)
        l_settings.addWidget(grp_input)

        # Optimised Settings / Quality Goals
        grp_qual = QGroupBox("")
        v_qual = QVBoxLayout()
        v_qual.setSpacing(0)
        v_qual.setContentsMargins(2, 2, 2, 2)
        
        lbl_qual_title = QLabel("Optimisation Parameters")
        lbl_qual_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_qual_title.setAlignment(Qt.AlignCenter)
        v_qual.addWidget(lbl_qual_title)
        
        self.grid_qual = QGridLayout()
        self.grid_qual.setContentsMargins(2, 2, 2, 2)
        self.grid_qual.setSpacing(5)

        # Row 0: Sync Strategy
        lbl_sync_strat = QLabel("Sync Strategy:", self)
        lbl_sync_strat.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(lbl_sync_strat, 0, 0)
        
        self.cmb_sync_strategy = QComboBox(self)
        self.cmb_sync_strategy.addItems([
            "Adaptive Integer Sync (Auto)",
            "Strict Pulse Target (Manual)",
            "Oversampling (Time-Driven)",
            "Combined Integer Sync (Auto)"
        ])
        self.cmb_sync_strategy.setToolTip(
            "Select how the laser repetition rate and acquisition time are optimized:\n\n"
            "• Adaptive Integer Sync (Auto): Automatically sweeps rep-rates and acquisition times "
            "to guarantee a perfect integer pulse count sync (0% rounding error) close to your target.\n"
            "• Strict Pulse Target (Manual): Directly calculates acquisition times for your exact pulse target. "
            "Note: Hardware rounding may result in non-integer pulse counts (unsynchronized).\n"
            "• Oversampling (Time-Driven): Focuses purely on meeting the washout-time and target RSD constraints "
            "for steady-state oversampling, ignoring integer pulse snapping.\n"
            "• Combined Integer Sync (Auto): Combines both behaviors. It sweeps rep-rates above the oversampling minimum "
            "to find and guarantee a perfect integer pulse count sync."
        )
        saved_strat = self.persistent_settings.get('sync_strategy', 'Adaptive Integer Sync (Auto)')
        mapping = {
            "Synchronised (Auto)": "Adaptive Integer Sync (Auto)",
            "Synchronised (Manual)": "Strict Pulse Target (Manual)",
            "Oversampling": "Oversampling (Time-Driven)",
            "Combined (Auto)": "Combined Integer Sync (Auto)",
            "Combined (Manual)": "Strict Pulse Target (Manual)"
        }
        saved_strat = mapping.get(saved_strat, saved_strat)
        self.cmb_sync_strategy.setCurrentText(saved_strat)
        self.cmb_sync_strategy.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.cmb_sync_strategy.currentTextChanged.connect(self._on_strategy_changed)
        self.grid_qual.addWidget(self.cmb_sync_strategy, 0, 1, 1, 4)

        # Row 1: Pulses per Acq. Time and Min Duty Cycle
        self.lbl_pulses = QLabel("Sync Target (Pulses):"); self.lbl_pulses.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(self.lbl_pulses, 1, 0)
        self.spin_pulses = QDoubleSpinBox()
        self.spin_pulses.setRange(1, 100); self.spin_pulses.setValue(self.persistent_settings.get('pulses', 5)); self.spin_pulses.setDecimals(1)
        self.spin_pulses.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_pulses, 1, 1)

        # RSD Target (Oversampling Mode counterpart, placed in same grid cell)
        self.lbl_rsd = QLabel("Sync Target (RSD %):"); self.lbl_rsd.setFixedWidth(lbl_w_in)
        self.lbl_rsd.setToolTip("Target RSD % (relative standard deviation limit) for oversampling steady state")
        self.lbl_rsd.setVisible(False)
        self.grid_qual.addWidget(self.lbl_rsd, 1, 0)
        
        self.spin_rsd = QDoubleSpinBox()
        self.spin_rsd.setRange(0.1, 20.0); self.spin_rsd.setValue(self.persistent_settings.get('target_rsd', 5.0)); self.spin_rsd.setDecimals(1)
        self.spin_rsd.setFixedWidth(spin_w_in)
        self.spin_rsd.setToolTip("Target RSD % (relative standard deviation limit) for oversampling steady state")
        self.spin_rsd.setVisible(False)
        self.spin_rsd.valueChanged.connect(self._on_ui_change)
        self.grid_qual.addWidget(self.spin_rsd, 1, 1)

        lbl_duty = QLabel("Min Duty Cycle (%):"); lbl_duty.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(lbl_duty, 1, 3)
        self.spin_duty = QDoubleSpinBox()
        self.spin_duty.setRange(0, 100); self.spin_duty.setValue(self.persistent_settings.get('min_duty', 0)); self.spin_duty.setDecimals(1)
        self.spin_duty.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_duty, 1, 4)

        # Row 2: Target SNR and Minimum SNR
        lbl_sigma = QLabel("Target SNR (Sigma):"); lbl_sigma.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(lbl_sigma, 2, 0)
        self.spin_sigma = QDoubleSpinBox()
        self.spin_sigma.setRange(0.1, 100000); self.spin_sigma.setValue(self.persistent_settings.get('target_sigma', 10)); self.spin_sigma.setDecimals(1)
        self.spin_sigma.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_sigma, 2, 1)

        lbl_snr = QLabel("Minimum SNR (Sigma):"); lbl_snr.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(lbl_snr, 2, 3)
        self.spin_snr = QDoubleSpinBox()
        self.spin_snr.setRange(0, 1000); self.spin_snr.setValue(self.persistent_settings.get('min_snr', 0)); self.spin_snr.setDecimals(1)
        self.spin_snr.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_snr, 2, 4)

        # Row 3: Image Width and Height / Shots
        self.lbl_width = QLabel("Image Width (µm):"); self.lbl_width.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(self.lbl_width, 3, 0)
        self.spin_width = QDoubleSpinBox()
        self.spin_width.setRange(1, 1000000); self.spin_width.setValue(self.persistent_settings.get('img_width', 1000.0)); self.spin_width.setDecimals(0)
        self.spin_width.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_width, 3, 1)

        # Line Length widget (Line Mode counterpart)
        self.lbl_line_len = QLabel("Line Length (µm):"); self.lbl_line_len.setFixedWidth(lbl_w_in)
        self.lbl_line_len.setVisible(False)
        self.grid_qual.addWidget(self.lbl_line_len, 3, 0)
        self.spin_line_len = QDoubleSpinBox()
        self.spin_line_len.setRange(1, 1000000); self.spin_line_len.setValue(self.persistent_settings.get('line_length', 1000.0)); self.spin_line_len.setDecimals(0)
        self.spin_line_len.setFixedWidth(spin_w_in)
        self.spin_line_len.setVisible(False)
        self.spin_line_len.valueChanged.connect(self._on_ui_change)
        self.grid_qual.addWidget(self.spin_line_len, 3, 1)

        self.lbl_height = QLabel("Image Height (µm):"); self.lbl_height.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(self.lbl_height, 3, 3)
        self.spin_height = QDoubleSpinBox()
        self.spin_height.setRange(1, 1000000); self.spin_height.setValue(self.persistent_settings.get('img_height', 1000.0)); self.spin_height.setDecimals(0)
        self.spin_height.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_height, 3, 4)
        
        # Shots (Overlapping Height in Layout, visibility toggled)
        self.lbl_shots = QLabel("Number of Shots:"); self.lbl_shots.setFixedWidth(lbl_w_in)
        # self.grid_qual.addWidget(self.lbl_shots, 3, 0) # Managed by mode update
        self.spin_shots = QSpinBox()
        self.spin_shots.setRange(1, 100000000); self.spin_shots.setValue(int(self.persistent_settings.get('shots', 100)))
        self.spin_shots.setFixedWidth(spin_w_in)
        # self.grid_qual.addWidget(self.spin_shots, 3, 1) # Managed by mode update

        # Row 4: Sync Dosage Checkbox and Dosage Input (UI elements commented out but preserved)
        self.chk_sync_dosage = QCheckBox("Sync Dosage")
        self.chk_sync_dosage.setToolTip("Synchronize Dosage with Pulses per Acq. Time")
        self.chk_sync_dosage.setChecked(True) # Leave sync enabled by default in code
        # self.grid_qual.addWidget(self.chk_sync_dosage, 4, 0, 1, 2)

        lbl_dosage = QLabel("Dosage (Shots/Area):"); lbl_dosage.setFixedWidth(lbl_w_in)
        # self.grid_qual.addWidget(lbl_dosage, 4, 2)
        self.spin_dosage = QDoubleSpinBox()
        self.spin_dosage.setRange(1, 100); self.spin_dosage.setValue(self.persistent_settings.get('dosage', self.persistent_settings.get('pulses', 5))); self.spin_dosage.setDecimals(1)
        self.spin_dosage.setFixedWidth(spin_w_in)
        # self.grid_qual.addWidget(self.spin_dosage, 4, 3)

        # Row 5: Scale Signal with Rep-Rate and Avoid Gaps
        self.chk_scale = QCheckBox("Scale Signal with Rep-Rate", self)
        self.chk_scale.setToolTip("Scale sensitivity based on Optimised Rep-Rate vs the Initial Rep-Rate. An increase in Rep-Rate will generally increase sensitivity, provided there is enough material present.")
        self.chk_scale.setChecked(self.persistent_settings.get('scale_signal', False))
        self.grid_qual.addWidget(self.chk_scale, 5, 0, 1, 2)

        self.chk_avoid_gaps = QCheckBox("Avoid Gaps", self)
        self.chk_avoid_gaps.setToolTip("Checked: Always Overlap (Increase Rep-Rate). Unchecked: Allow Gaps.")
        self.chk_avoid_gaps.setChecked(self.persistent_settings.get('avoid_gaps', False))
        self.chk_avoid_gaps.toggled.connect(self._on_ui_change)
        self.grid_qual.addWidget(self.chk_avoid_gaps, 5, 3, 1, 2)

        # Row 6: Prefer Exact Pulses per Acq
        self.chk_prefer_exact_pulses = QCheckBox("Prefer Exact Pulses per Acq", self)
        self.chk_prefer_exact_pulses.setToolTip(
            "<b>Prefer Exact Pulses per Acq (Auto modes only)</b><br/>"
            "When checked, the auto-sync search prioritises finding a rep-rate that delivers "
            "exactly the requested number of pulses per acquisition cycle.<br/><br/>"
            "<b>Trade-off:</b> This may require a larger change in rep-rate and a longer "
            "acquisition time compared to the closest achievable integer sync, which can "
            "reduce throughput or increase overlap."
        )
        self.chk_prefer_exact_pulses.setChecked(
            self.persistent_settings.get('prefer_exact_pulses', False))
        self.chk_prefer_exact_pulses.toggled.connect(self._on_ui_change)
        self.grid_qual.addWidget(self.chk_prefer_exact_pulses, 5, 3, 1, 2)

        # Let Column 2 act as the expanding center space between Column 1 and Column 3
        self.grid_qual.setColumnStretch(2, 1)

        v_qual.addLayout(self.grid_qual)
        grp_qual.setLayout(v_qual)
        l_settings.addWidget(grp_qual)

        # Optimised Settings (Moved from Results)
        grp_sync = QGroupBox("")
        grp_sync.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        v_sync = QVBoxLayout()
        v_sync.setSpacing(0)
        v_sync.setContentsMargins(2, 2, 2, 2)
        
        lbl_sync_title = QLabel("Optimised Settings")
        lbl_sync_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_sync_title.setAlignment(Qt.AlignCenter)
        v_sync.addWidget(lbl_sync_title)
        
        l_sync = QHBoxLayout()
        # Column 1: Spot Size (Left)
        # Core Optimised Values (Unified Grid)
        grid_sync = QGridLayout()
        grid_sync.setContentsMargins(2, 2, 2, 2)
        grid_sync.setSpacing(5)
        grid_sync.setVerticalSpacing(1)
        
        # Force all rows to have same height to match the SpinBox in Row 0
        for r in range(4):
            grid_sync.setRowMinimumHeight(r, 22)

        # -- Column Groups (0,1), (2,3), (4,5) --

        # Row 0: Spot Size, Rep-Rate, Acq Time
        grid_sync.addWidget(QLabel("Spot Size (µm):"), 0, 0)
        self.spin_opt_spot = QSpinBox()
        self.spin_opt_spot.setRange(0, 500); self.spin_opt_spot.setFixedWidth(55)
        self.spin_opt_spot.setSpecialValueText("Auto")
        self.spin_opt_spot.valueChanged.connect(self._on_opt_spot_changed)
        grid_sync.addWidget(self.spin_opt_spot, 0, 1)

        grid_sync.addWidget(QLabel("Rep-Rate:"), 0, 2)
        self.lbl_opt_rr = QLabel("- Hz")
        grid_sync.addWidget(self.lbl_opt_rr, 0, 3)

        grid_sync.addWidget(QLabel("Acq Time:"), 0, 4)
        self.lbl_opt_at = QLabel("- ms")
        grid_sync.addWidget(self.lbl_opt_at, 0, 5)

        # Row 1: Est. Time (hms), Speed, Overhead
        grid_sync.addWidget(QLabel("Est. Time:"), 1, 0)
        self.lbl_est_time_hms = QLabel("-")
        grid_sync.addWidget(self.lbl_est_time_hms, 1, 1)

        self.lbl_opt_speed_title = QLabel("Speed:")
        grid_sync.addWidget(self.lbl_opt_speed_title, 1, 2)
        self.lbl_opt_speed = QLabel("- µm s⁻¹")
        grid_sync.addWidget(self.lbl_opt_speed, 1, 3)

        grid_sync.addWidget(QLabel("Overhead:"), 1, 4)
        self.lbl_opt_overhead = QLabel("- ms")
        grid_sync.addWidget(self.lbl_opt_overhead, 1, 5)

        # Row 2: Est. Time (seconds), Overlap, Budget
        self.lbl_est_time_sec = QLabel("-")
        grid_sync.addWidget(self.lbl_est_time_sec, 2, 1)

        self.lbl_opt_overlap_title = QLabel("Overlap:")
        grid_sync.addWidget(self.lbl_opt_overlap_title, 2, 2)
        self.lbl_opt_overlap = QLabel("- µm")
        grid_sync.addWidget(self.lbl_opt_overlap, 2, 3)

        grid_sync.addWidget(QLabel("Budget:"), 2, 4)
        self.lbl_opt_budget = QLabel("- ms")
        grid_sync.addWidget(self.lbl_opt_budget, 2, 5)

        # Row 3: Est. Data Points / Pixels Per Sec, Dosage, Duty Cycle
        self.lbl_est_datapoints_title = QLabel("Est. Data Points:")
        grid_sync.addWidget(self.lbl_est_datapoints_title, 3, 0)
        self.lbl_est_datapoints = QLabel("-")
        grid_sync.addWidget(self.lbl_est_datapoints, 3, 1)

        self.lbl_opt_dosage_title = QLabel("Dosage:")
        grid_sync.addWidget(self.lbl_opt_dosage_title, 3, 2)
        self.lbl_opt_pulses = QLabel("- Pulses")
        grid_sync.addWidget(self.lbl_opt_pulses, 3, 3)

        grid_sync.addWidget(QLabel("Duty Cycle:"), 3, 4)
        self.lbl_opt_duty = QLabel("- %")
        grid_sync.addWidget(self.lbl_opt_duty, 3, 5)

        # Set stretch and minimum widths for all columns to prevent layout shifting when toggling modes/labels
        grid_sync.setColumnMinimumWidth(0, 85)
        grid_sync.setColumnMinimumWidth(1, 60)
        grid_sync.setColumnMinimumWidth(2, 45)
        grid_sync.setColumnMinimumWidth(3, 85)
        grid_sync.setColumnMinimumWidth(4, 65)
        grid_sync.setColumnMinimumWidth(5, 60)
        grid_sync.setColumnStretch(1, 1)
        grid_sync.setColumnStretch(3, 1)
        v_sync.addLayout(grid_sync)
        grp_sync.setLayout(v_sync)
        l_settings.addWidget(grp_sync)
        
        self.lbl_result = QLabel("Ready")
        self.lbl_result.setWordWrap(True)
        self.lbl_result.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.lbl_result.setStyleSheet("padding: 5px;")
        
        self.result_scroll = QScrollArea()
        self.result_scroll.setWidget(self.lbl_result)
        self.result_scroll.setWidgetResizable(True)
        self.result_scroll.setMinimumHeight(60) # Reduced from 120
        # Using 0 for NoFrame to avoid extra imports if not present, but usually QFrame is there.
        # Let's check imports or just use QScrollArea.NoFrame if we know it exists.
        # In this codebase QScrollArea is from iolite.QtGui or PyQt5.QtWidgets.
        self.result_scroll.setFrameShape(0) # 0 is QFrame.NoFrame
        self.result_scroll.setStyleSheet("background-color: transparent;")
        
        l_settings.addWidget(self.result_scroll, 1) # Assign stretch of 1 to allow shrinking/expansion
        
        # l_settings.addStretch() # Removed static stretch to allow result_scroll to manage space
        self.settings_scroll_area.setWidget(scroll_content)
        self.settings_scroll_area.setFixedWidth(460)
        self.settings_scroll_area.setWidgetResizable(True)
        self.settings_scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff) # Main panel scrolling disabled
        left_layout.addWidget(self.settings_scroll_area)
        
        # --- RIGHT COLUMN (PLOT & RESULTS) ---
        # Plot Panel (Top)
        plot_widget = QWidget()
        plot_widget.setMinimumHeight(320) # Allow more flexible resizing
        plot_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        plot_layout = QVBoxLayout(plot_widget)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.setSpacing(2)
        
        right_layout = QVBoxLayout()
        right_layout.setContentsMargins(0, 0, 0, 0)
        
        self.figure = Figure(figsize=(5, 3), dpi=96, constrained_layout=False); self.canvas = FigureCanvas(self.figure)
        # Initialize with baseline margins (match resize logic)
        self.figure.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15)
        self.canvas.mpl_connect('resize_event', self._on_plot_resize)
        self.canvas.setMinimumHeight(200) # Allow smaller resize
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.updateGeometry()
        self.canvas.mpl_connect('scroll_event', self.on_zoom)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_drag)
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.plot_panning = False; self.press_x = None
        
        h_ctrl1 = QHBoxLayout()
        h_ctrl2 = QHBoxLayout()
        self.chk_norm = QCheckBox("Normalise")
        self.chk_norm.setChecked(self.persistent_settings.get('normalise', True))
        
        self.chk_y_zoom = QCheckBox("Pan / Zoom Y")
        self.chk_y_zoom.setChecked(self.persistent_settings.get('opt_y_zoom', False))
        self.chk_y_zoom.toggled.connect(self._on_opt_y_zoom_toggled)
        
        self.chk_rescale = QCheckBox("Auto-Rescale Y")
        self.chk_rescale.setChecked(True) # Ephemeral: Always default to ON
        self.chk_rescale.toggled.connect(self._on_opt_rescale_toggled)
        
        # Theme Override
        self.combo_theme = QComboBox()
        self.combo_theme.addItems(["Auto", "Dark", "Light"])
        self.combo_theme.setCurrentText(self.persistent_settings.get('theme', 'Auto'))
        
        # Region Controls (Moved to Plot Toolbar)
        self.chk_auto = QCheckBox("Auto Select Regions")
        self.chk_auto.setChecked(True) # Default On
        
        self.spin_bg_start = QDoubleSpinBox(); self.spin_bg_start.setRange(0, 10000); self.spin_bg_start.setDecimals(2); self.spin_bg_start.setFixedWidth(60)
        self.spin_bg_end = QDoubleSpinBox(); self.spin_bg_end.setRange(0, 10000); self.spin_bg_end.setDecimals(2); self.spin_bg_end.setFixedWidth(60)
        self.spin_sig_start = QDoubleSpinBox(); self.spin_sig_start.setRange(0, 10000); self.spin_sig_start.setDecimals(2); self.spin_sig_start.setFixedWidth(60)
        self.spin_sig_end = QDoubleSpinBox(); self.spin_sig_end.setRange(0, 10000); self.spin_sig_end.setDecimals(2); self.spin_sig_end.setFixedWidth(60)
        
        h_ctrl1.addWidget(QLabel("Theme:")); h_ctrl1.addWidget(self.combo_theme)
        h_ctrl1.addSpacing(10)
        h_ctrl1.addWidget(self.chk_norm); h_ctrl1.addWidget(self.chk_y_zoom)
        h_ctrl1.addWidget(self.chk_rescale)
        h_ctrl1.addStretch()
        
        h_ctrl2.addWidget(self.chk_auto)
        h_ctrl2.addWidget(QLabel("Background (s):"))
        h_ctrl2.addWidget(self.spin_bg_start); h_ctrl2.addWidget(QLabel("to")); h_ctrl2.addWidget(self.spin_bg_end)
        h_ctrl2.addSpacing(10)
        h_ctrl2.addWidget(QLabel("Signal (s):"))
        h_ctrl2.addWidget(self.spin_sig_start); h_ctrl2.addWidget(QLabel("to")); h_ctrl2.addWidget(self.spin_sig_end)
        h_ctrl2.addStretch()
        
        plot_layout.addLayout(h_ctrl1)
        plot_layout.addLayout(h_ctrl2)
        plot_layout.addWidget(self.canvas)
        
        # Results Panel (Bottom)
        results_widget = QWidget()
        results_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        results_layout = QVBoxLayout(results_widget)
        results_layout.setContentsMargins(0, 0, 0, 0)
        
        self.setStyleSheet("QGroupBox::title { font-weight: bold; font-size: 14px; }")

        # Optimised Dwell Time Distribution (Table)
        grp_table = QGroupBox("")
        v_table_main = QVBoxLayout()
        v_table_main.setSpacing(0)
        v_table_main.setContentsMargins(5, 2, 5, 2)
        
        lbl_table_title = QLabel("Optimised Dwell Time Distribution")
        lbl_table_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_table_title.setAlignment(Qt.AlignCenter)
        v_table_main.addWidget(lbl_table_title)
        
        l_table = QVBoxLayout()
        l_table.setContentsMargins(5, 5, 5, 5)
        
        self.table = CopyableTableWidget()
        self.table.setColumnCount(5)
        self.table.setStyleSheet("QHeaderView::section { padding-left: 10px; padding-right: 10px; }")
        # Headers will be set dynamically in run_optimisation
        self.table.setHorizontalHeaderLabels(["Isotope", "Mode", "Optimised Dwell Times (ms)", "Initial SNR", "Resultant SNR"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        l_table.addWidget(self.table)
        
        grp_table.setLayout(v_table_main)
        v_table_main.addLayout(l_table)
        results_layout.addWidget(grp_table)
        

        
        # Splitter
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(plot_widget)
        splitter.addWidget(results_widget)
        splitter.setStretchFactor(0, 1) 
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([500, 500]) # Balanced initial distribution
        
        right_layout.addWidget(splitter)
        
        main_layout.addLayout(left_layout, 1)
        main_layout.addLayout(right_layout, 2)
        self.tab_opt.setLayout(main_layout)
        
        
        # 1. Sync internal state with loaded settings (Silent population)
        IoLog.information("iolite Optimiser: Applying hardware settings...")
        self._handle_mfr_changed(self._get(self.cmb_mfr, 'currentText'))
        self._handle_laser_mfr_changed(self._get(self.cmb_laser_mfr, 'currentText'))
        self._update_hw_summary()
        self._update_ui_precisions()

        # 2. Connect signals after initial population to prevent startup save/opt spam
        IoLog.information("iolite Optimiser: Connecting interactive signals...")
        self.cmb_mfr.currentTextChanged.connect(self._handle_mfr_changed)
        self.cmb_mfr.currentTextChanged.connect(self._on_ui_change)
        
        self.cmb_model.currentTextChanged.connect(self._handle_model_changed)
        self.cmb_model.currentTextChanged.connect(self._on_ui_change)
        
        self.cmb_laser_mfr.currentTextChanged.connect(self._handle_laser_mfr_changed)
        self.cmb_laser_mfr.currentTextChanged.connect(self._on_ui_change)
        
        self.cmb_laser_mod.currentTextChanged.connect(self._handle_laser_mod_changed)
        self.cmb_laser_mod.currentTextChanged.connect(self._on_ui_change)

        self.cmb_cell.currentTextChanged.connect(self._handle_cell_changed)
        self.cmb_cell.currentTextChanged.connect(self._on_ui_change)
        
        self.spin_spot.valueChanged.connect(self._on_ui_change)
        self.spin_num_analytes.valueChanged.connect(self.save_persistent_settings)
        self.spin_height.valueChanged.connect(self._on_ui_change)
        self.spin_width.valueChanged.connect(self._on_ui_change)
        self.spin_wash.valueChanged.connect(self._on_ui_change)
        self.spin_init_rr.valueChanged.connect(self._on_ui_change)
        self.cmb_mode.currentTextChanged.connect(self._on_ui_change)
        self.spin_shots.valueChanged.connect(self._on_ui_change)
        self.chk_scale.toggled.connect(self._on_ui_change)
        self.spin_pulses.valueChanged.connect(self._on_pulses_changed)
        self.spin_dosage.valueChanged.connect(self._on_dosage_changed)
        self.chk_sync_dosage.toggled.connect(self._on_sync_dosage_toggled)
        self.spin_sigma.valueChanged.connect(self._on_ui_change)
        self.spin_snr.valueChanged.connect(self._on_ui_change)
        self.spin_duty.valueChanged.connect(self._on_ui_change)
        self.chk_norm.toggled.connect(lambda: self.update_plot(preserve_x_only=True))
        self.chk_norm.toggled.connect(self._on_ui_change)
        # self.chk_y_zoom.toggled.connect(self._on_ui_change) # Removed to prevent rescale on toggle
        self.chk_rescale.toggled.connect(self._on_ui_change)
        self.combo_theme.currentTextChanged.connect(self.apply_theme)
        self.chk_auto.toggled.connect(self.on_auto_toggled)
        self.chk_auto.toggled.connect(self._on_ui_change)
        for sb in [self.spin_bg_start, self.spin_bg_end, self.spin_sig_start, self.spin_sig_end]:
            sb.valueChanged.connect(self.on_region_edited)

        # Custom Field Signals
        self.cmb_cust_type.currentTextChanged.connect(self._handle_cust_type_changed)
        self.cmb_cust_type.currentIndexChanged.connect(self._on_ui_change)
        self.spin_cust_min.valueChanged.connect(self._on_ui_change)
        self.spin_cust_prec.valueChanged.connect(self._on_ui_change)
        self.edit_cust_dwells.textChanged.connect(self._on_ui_change)
        self.cmb_cust_rr_type.currentIndexChanged.connect(self._on_ui_change)
        self.spin_cust_rr.valueChanged.connect(self._on_ui_change)
        self.spin_cust_rr_prec.valueChanged.connect(self._on_ui_change)
        self.edit_cust_rr_list.textChanged.connect(self._on_ui_change)
        self.spin_cust_speed.valueChanged.connect(self._on_ui_change)
        
        # Trigger initial visibility check
        self._update_custom_visibility()
        self._update_mode_visibility()
        
        # Connect Result Label Link
        self.lbl_result.setTextInteractionFlags(Qt.LinksAccessibleByMouse | Qt.TextSelectableByMouse)
        self.lbl_result.setOpenExternalLinks(False)
        self.lbl_result.linkActivated.connect(self._on_result_link_clicked)
            
        # 2.5. Apply Theme EARLY (Before Opt)
        # Initializes the plot background color and system_state
        self.apply_theme(trigger_opt=False)

        # 3. Initial Data Loading & Optimization
        IoLog.information("iolite Optimiser: Starting initial data load...")
        self.refresh_data()
        
        if self.opt_df is not None:
            if self.bg_times is not None and self.sig_times is not None:
                IoLog.information("iolite Optimiser: Initial regions valid. Optimisation triggered by refresh_data.")
                # self.run_optimisation(refresh=False) -> HANDLED BY refresh_data()
            else:
                IoLog.information("iolite Optimiser: Initial regions empty (Auto-Detect required).")
        else:
            IoLog.information("iolite Optimiser: No data found on startup.")
            
        IoLog.information("iolite Optimiser: UI initialization complete.")
    
    def _on_result_link_clicked(self, link):
        if link == "#set_dwells":
            IoLog.information("iolite Optimiser: User clicked 'Set Dwell Times' link.")
            # Determine missing channels again or just force open for all
            if hasattr(self, 'opt_df') and self.opt_df is not None:
                self.resolve_dwell_times(set(self.opt_df.columns) - {'Time'})
                
                # RECALCULATE OVERHEAD with new dwell times
                if hasattr(self, 'detected_at_ms') and hasattr(self, 'channel_dwells'):
                    total_init_dwell_ms = self._get_total_dwell_time_ms()
                    self.detected_overhead_ms = max(0, self.detected_at_ms - total_init_dwell_ms)
                    IoLog.information(f"iolite Optimiser: Overhead updated to {self.detected_overhead_ms:.3f}ms")

                # Re-run optimization to see if fixed
                self.run_optimisation(refresh=False)

    def _get_total_dwell_time_ms(self):
        """
        Calculates the effective measurement time per individual analysis cycle.
        For Quads/Sector: Sum of all dwell times (sequential).
        For MC/TOF: Maximum dwell time (simultaneous).
        """
        # Ensure source exists (Opt tab priority)
        dwells_source = getattr(self, 'opt_dwells', getattr(self, 'channel_dwells', {}))

        if not dwells_source:
            return 0.0
        
        # Robust Tech Check
        # TYPE_TO_TECH_MAP = {"Quadrupole": "Quad", "Multi-Collector": "MC", "TOF": "TOF", ...}
        tech_raw = self.icp_tech
        mapped = TYPE_TO_TECH_MAP.get(tech_raw, tech_raw)
        
        # Simultaneous Systems (Max Dwell)
        if mapped in ["Multi-Collector", "TOF"]:
            return max(dwells_source.values())
        # Sequential Systems (Sum Dwells)
        else:
            return sum(dwells_source.values())

    def _recalculate_overhead(self):
        """
        Recalculates the overhead based on the CURRENT model/tech.
        Called on Load and when Model changes.
        """
        if not hasattr(self, 'detected_at_ms') or self.detected_at_ms is None:
            return

        total_init_dwell_ms = self._get_total_dwell_time_ms()
        
        # Precise calculation preserves imported file timing
        self.detected_overhead_ms = max(0, self.detected_at_ms - total_init_dwell_ms)
        
        tech_raw = getattr(self, 'icp_tech', 'Quadrupole')
        mapped = TYPE_TO_TECH_MAP.get(tech_raw, tech_raw)
        dwell_label = "Max Dwell" if mapped in ["Multi-Collector", "TOF"] else "Sum of Dwells"
        
        IoLog.information(f"iolite Optimiser: Overhead updated to {self.detected_overhead_ms:.3f}ms (AT={self.detected_at_ms:.3f}, {dwell_label}={total_init_dwell_ms:.3f})")
        
        # Re-run optimization to refresh results if needed
        # self.run_optimisation(refresh=False) # Optional: might be too aggressive on every change?

    def _handle_model_changed(self, model):
        mfr = self._get(self.cmb_mfr, 'currentText')
        is_custom = (model == "Custom Model" or mfr == "Custom")
        
        # Toggle Custom Settings Visibility (individual widgets in form)
        is_mc = (getattr(self, 'icp_tech', 'Quad') == "Multi-Collector")
        is_custom_icp = is_custom
        
        self.lbl_cust_type.setVisible(is_custom_icp)
        self.cmb_cust_type.setVisible(is_custom_icp)
        self.lbl_cust_min.setVisible(is_custom_icp)
        self.spin_cust_min.setVisible(is_custom_icp)
        self.lbl_cust_prec.setVisible(is_custom_icp)
        self.spin_cust_prec.setVisible(is_custom_icp)
        self.lbl_cust_dwells.setVisible(is_custom_icp and is_mc)
        self.edit_cust_dwells.setVisible(is_custom_icp and is_mc)
            
        if is_custom:
            self.icp_tech = self._get(self.cmb_cust_type, 'currentText')
            self.min_dwell = self._get(self.spin_cust_min, 'value')
            self.precision = self._get(self.spin_cust_prec, 'value')
            if self.icp_tech == "Multi-Collector":
                self.allowed_dwells = self._parse_rr_list(self._get(self.edit_cust_dwells, 'text'))
            else:
                self.allowed_dwells = None
        else:
            spec = ICP_SPECS.get(mfr, {}).get(model, {})
            self.min_dwell = spec.get('min_dwell', 0.1)
            self.precision = spec.get('prec', 0.1)
            
            # CRITICAL: Update Tech Type (Quad vs MC) from Spec
            raw_type = spec.get('type', 'Quadrupole')
            self.icp_tech = TYPE_TO_TECH_MAP.get(raw_type, raw_type)
            self.allowed_dwells = spec.get('allowed_dwells', None)
        
        # Always update status
        self._update_icp_status()
        self.lbl_result.setText(f"System: {mfr} {model} ({self.icp_tech})")
        
        # Persistence
        if mfr:
            self.persistent_settings[f"last_model_for_{mfr}"] = model
            
        # Trigger Overhead Recalculation
        self._recalculate_overhead()

    def _handle_cust_type_changed(self, text):
        self._handle_model_changed(self._get(self.cmb_model, 'currentText'))

    def _update_icp_status(self):
        if self.allowed_dwells:
            status = f"Allowed Dwells: {', '.join(map(str, self.allowed_dwells))} ms"
        else:
            status = f"Minimum Dwell Time: {self.min_dwell} ms | Precision: {self.precision} ms"
            
        if getattr(self, 'icp_tech', '') == "TOF":
            margin = self._get(self.spin_tof_margin, 'value') if hasattr(self, 'spin_tof_margin') else getattr(self, 'persistent_settings', {}).get('tof_margin', 10.0)
            status += f" | Margin: {margin:g} %"
            
        self.lbl_icp_status.setText(status)

    def _handle_laser_mod_changed(self, model):
        mfr = self._get(self.cmb_laser_mfr, 'currentText')
        # Only show overrides if the MODEL itself is "Custom Laser"
        is_custom_laser = (model == "Custom Laser")
        
        # Check if stage is custom
        cell = self._get(self.cmb_cell, 'currentText')
        is_custom_stage = (cell == "Custom Stage")
        
        # Toggle Custom Settings Visibility (individual widgets in form)
        is_hybrid_visible = (is_custom_laser or is_custom_stage)
        
        # Mode is only for Custom Laser
        self.lbl_cust_mode.setVisible(is_custom_laser)
        self.cmb_cust_rr_type.setVisible(is_custom_laser)
        
        if is_custom_laser:
            rr_type = self._get(self.cmb_cust_rr_type, 'currentText')
            is_discrete = (rr_type == "Discrete")
            
            self.lbl_cust_rr_max.setVisible(not is_discrete)
            self.spin_cust_rr.setVisible(not is_discrete)
            self.lbl_cust_rr_list.setVisible(is_discrete)
            self.edit_cust_rr_list.setVisible(is_discrete)
            self.lbl_cust_rr_prec.setVisible(True)
            self.spin_cust_rr_prec.setVisible(True)
            
            if not is_discrete:
                self.max_rr = self._get(self.spin_cust_rr, 'value')
                self.allowed_rr = None
                self.laser_rr_prec = self._get(self.spin_cust_rr_prec, 'value')
            else:
                self.max_rr = 10000
                self.allowed_rr = self._parse_rr_list(self._get(self.edit_cust_rr_list, 'text'))
                self.laser_rr_prec = 1.0 # Discrete lists don't use precision rounding

        else:
            # Not custom laser - Hide laser custom fields
            self.lbl_cust_rr_max.setVisible(False)
            self.spin_cust_rr.setVisible(False)
            self.lbl_cust_rr_list.setVisible(False)
            self.edit_cust_rr_list.setVisible(False)
            self.lbl_cust_rr_prec.setVisible(False)
            self.spin_cust_rr_prec.setVisible(False)
            
            # Retrieve spec from master DB
            spec = LASER_SOURCES.get(model, {})
            self.max_rr = spec.get('max_rr', 10)
            self.allowed_rr = spec.get('allowed_rr')
            self.laser_rr_prec = spec.get('rr_prec', 1.0)
            
        # Update status and check stage speed visibility
        self._handle_cell_changed(cell)
        self._update_laser_status()


    def _handle_cell_changed(self, cell):
        is_custom_stage = (cell == "Custom Stage")
        
        # Mode is only for Custom Laser
        las_mod = self._get(self.cmb_laser_mod, 'currentText')
        is_custom_laser = (las_mod == "Custom Laser")
        if hasattr(self, 'lbl_cust_mode'):
            self.lbl_cust_mode.setVisible(is_custom_laser)
            self.cmb_cust_rr_type.setVisible(is_custom_laser)

        if hasattr(self, 'lbl_cust_speed'):
            self.lbl_cust_speed.setVisible(is_custom_stage)
            self.spin_cust_speed.setVisible(is_custom_stage)
            
        if is_custom_stage:
            self.max_speed = self._get(self.spin_cust_speed, 'value')
        else:
            self.max_speed = STAGE_SPECS.get(cell, 2000)
        self._update_laser_status()
        self._update_hw_summary()

    def _update_laser_status(self):
        # Build status string
        if self.allowed_rr:
            status = f"Rep-Rates: {', '.join(map(str, self.allowed_rr))} Hz"
        else:
            status = f"Max Rep-Rate: {self.max_rr} Hz"
            
        status += f" | Precision: {getattr(self, 'laser_rr_prec', 1.0):g} Hz"
        status += f"\nMax Speed: {self.max_speed} µm s⁻¹"
        self.lbl_laser_status.setText(status)

    def show_settings_dialog(self):
        self.settings_dlg.exec_()
        # Apply any plot-related settings that changed (fresh, don't preserve visibility)
        if hasattr(self, 'opt_df') and self.opt_df is not None:
            self.update_plot()
        
    def _update_custom_visibility(self):
        # 1. Gather current state
        icp_mfr = self._get(self.cmb_mfr, 'currentText')
        icp_mod = self._get(self.cmb_model, 'currentText')
        is_custom_icp = (icp_mod == "Custom Model" or icp_mfr == "Custom")
        
        las_mfr = self._get(self.cmb_laser_mfr, 'currentText')
        las_mod = self._get(self.cmb_laser_mod, 'currentText')
        is_custom_laser = (las_mod == "Custom Laser")
        
        cell = self._get(self.cmb_cell, 'currentText')
        is_custom_stage = (cell == "Custom Stage")

        # --- ICP VISIBILITY ---
        if hasattr(self, 'cmb_cust_type'):
            icp_type = self._get(self.cmb_cust_type, 'currentText')
            is_mc = (icp_type == "Multi-Collector")
            
            # Type is only for custom
            self.lbl_cust_type.setVisible(is_custom_icp)
            self.cmb_cust_type.setVisible(is_custom_icp)
            
            # Show List if MC, Hide Min/Prec (ONLY if custom)
            self.lbl_cust_dwells.setVisible(is_custom_icp and is_mc)
            self.edit_cust_dwells.setVisible(is_custom_icp and is_mc)
            
            self.lbl_cust_min.setVisible(is_custom_icp and not is_mc)
            self.spin_cust_min.setVisible(is_custom_icp and not is_mc)
            self.lbl_cust_prec.setVisible(is_custom_icp and not is_mc)
            self.spin_cust_prec.setVisible(is_custom_icp and not is_mc)
            
        if is_custom_icp and hasattr(self, 'cmb_cust_type'):
            is_tof = (self._get(self.cmb_cust_type, 'currentText') == "TOF")
        else:
            is_tof = (getattr(self, 'icp_tech', 'Quadrupole') == "TOF")
            
        if hasattr(self, 'lbl_tof_margin'):
            self.lbl_tof_margin.setVisible(is_tof)
            self.spin_tof_margin.setVisible(is_tof)
            
        # --- LASER VISIBILITY ---
        if hasattr(self, 'cmb_cust_rr_type'):
            mode = self._get(self.cmb_cust_rr_type, 'currentText')
            is_discrete = (mode == "Discrete")
            
            # Mode is only for custom laser
            self.lbl_cust_mode.setVisible(is_custom_laser)
            self.cmb_cust_rr_type.setVisible(is_custom_laser)
            
            # Rep Rate Overrides (ONLY if custom laser)
            self.lbl_cust_rr_max.setVisible(is_custom_laser and not is_discrete)
            self.spin_cust_rr.setVisible(is_custom_laser and not is_discrete)
            self.lbl_cust_rr_list.setVisible(is_custom_laser and is_discrete)
            self.edit_cust_rr_list.setVisible(is_custom_laser and is_discrete)
            self.lbl_cust_rr_prec.setVisible(is_custom_laser)
            self.spin_cust_rr_prec.setVisible(is_custom_laser)
            
            # Stage Speed (Only if Custom Stage)
            self.lbl_cust_speed.setVisible(is_custom_stage)
            self.spin_cust_speed.setVisible(is_custom_stage)

    def _update_hw_summary(self):
        icp_mfr = self._get(self.cmb_mfr, 'currentText')
        icp_mod = self._get(self.cmb_model, 'currentText')
        las_mfr = self._get(self.cmb_laser_mfr, 'currentText')
        las_mod = self._get(self.cmb_laser_mod, 'currentText')
        cell = self._get(self.cmb_cell, 'currentText')

        # --- ICP STATE SYNC ---
        is_custom_icp = (icp_mfr == "Custom" or icp_mod == "Custom Model")
        if is_custom_icp:
            self.icp_tech = self._get(self.cmb_cust_type, 'currentText')
            self.min_dwell = self._get(self.spin_cust_min, 'value')
            self.precision = self._get(self.spin_cust_prec, 'value')
            if self.icp_tech == "Multi-Collector":
                self.allowed_dwells = self._parse_rr_list(self._get(self.edit_cust_dwells, 'text'))
            else:
                self.allowed_dwells = None
        else:
            spec = ICP_SPECS.get(icp_mfr, {}).get(icp_mod, {})
            raw_type = spec.get('type', "Quadrupole")
            self.icp_tech = TYPE_TO_TECH_MAP.get(raw_type, raw_type)
            self.min_dwell = spec.get('min_dwell', 0.1)
            self.precision = spec.get('prec', 0.1)
            self.allowed_dwells = spec.get('allowed_dwells', None)

        # --- LASER STATE SYNC ---
        is_custom_laser = (las_mod == "Custom Laser")
        if is_custom_laser:
            rr_type = self._get(self.cmb_cust_rr_type, 'currentText')
            self.laser_rr_prec = self._get(self.spin_cust_rr_prec, 'value')
            if rr_type == "Maximum Rep-Rate":
                self.max_rr = self._get(self.spin_cust_rr, 'value')
                self.allowed_rr = None
            else:
                self.max_rr = 10000
                self.allowed_rr = self._parse_rr_list(self._get(self.edit_cust_rr_list, 'text'))
        else:
            spec = LASER_SOURCES.get(las_mod, {})
            self.max_rr = spec.get('max_rr', 10)
            self.allowed_rr = spec.get('allowed_rr', None)
            self.laser_rr_prec = spec.get('rr_prec', 1)
        
        if cell == "Custom Stage":
            self.max_speed = self._get(self.spin_cust_speed, 'value')
        else:
            self.max_speed = STAGE_SPECS.get(cell, 2000)

        # --- UPDATE UI LABELS ---
        
        # 1. Dialog Labels (Detailed)
        if hasattr(self, 'lbl_icp_status'):
            icp_constraints = f"Allowed Dwell times: {', '.join(map(str, self.allowed_dwells))} ms" if self.allowed_dwells else f"Minimum Dwell Time: {self.min_dwell} ms | Precision: {self.precision} ms"
            if self.icp_tech == "TOF":
                margin = self._get(self.spin_tof_margin, 'value') if hasattr(self, 'spin_tof_margin') else getattr(self, 'persistent_settings', {}).get('tof_margin', 10.0)
                icp_constraints += f" | Margin: {margin:g} %"
            self.lbl_icp_status.setText(f"Type: {self.icp_tech} | {icp_constraints}")

        if hasattr(self, 'lbl_laser_status'):
            rr_info = f"Allowed RRs: {', '.join(map(str, self.allowed_rr))} Hz" if self.allowed_rr else f"Max RR: {self.max_rr} Hz"
            self.lbl_laser_status.setText(f"{rr_info} | Max Speed: {self.max_speed} µm/s | Prec: {self.laser_rr_prec} Hz")

        # 2. Main Window Summary (HTML formatted)
        if hasattr(self, 'lbl_icp_sum') and hasattr(self, 'lbl_las_sum'):
            icp_header = "<b style='font-size: 1.1em;'>ICP-MS</b>"
            icp_details = f"Manufacturer: <b>{icp_mfr}</b><br/>Model: <b>{icp_mod}</b><br/>Type: <b>{self.icp_tech}</b>"
            
            if self.allowed_dwells:
                icp_constraints = f"Allowed Dwell Times:<br/>{', '.join(map(str, self.allowed_dwells))} ms"
            else:
                icp_constraints = f"Minimum Dwell Time: {self.min_dwell} ms<br/>Dwell Precision: {self.precision} ms"
            
            if self.icp_tech == "TOF":
                margin = self._get(self.spin_tof_margin, 'value') if hasattr(self, 'spin_tof_margin') else getattr(self, 'persistent_settings', {}).get('tof_margin', 10.0)
                icp_constraints += f"<br/>Washout Margin: {margin:g} %"
            
            icp_full = f"{icp_header}<br/>{icp_details}<br/><span style='font-size: 11px;'><b><i>{icp_constraints}</i></b></span>"
            
            las_header = "<b style='font-size: 1.1em;'>Laser</b>"
            las_details = f"Platform: <b>{las_mfr}</b><br/>Laser Source: <b>{las_mod}</b><br/>Cell Type: <b>{cell}</b>"
            rr_info = f"Allowed Rep-Rates: {', '.join(map(str, self.allowed_rr))} Hz" if self.allowed_rr else f"Maximum Rep-Rate: {self.max_rr} Hz"
            las_constraints = f"{rr_info}<br/>Rep-Rate Precision: {self.laser_rr_prec:g} Hz<br/>Max Stage Speed: {self.max_speed} µm s⁻¹"
            las_full = f"{las_header}<br/>{las_details}<br/><span style='font-size: 11px;'><b><i>{las_constraints}</i></b></span>"
            
            self.lbl_icp_sum.setText(icp_full)
            self.lbl_las_sum.setText(las_full)

    def _on_strategy_changed(self, text):
        strategy = text
        is_oversampling_mode = "Oversampling" in strategy or "Combined" in strategy
        
        # Toggle visibility of the two independent target inputs
        if hasattr(self, 'lbl_pulses'): self.lbl_pulses.setVisible(not is_oversampling_mode)
        if hasattr(self, 'spin_pulses'): self.spin_pulses.setVisible(not is_oversampling_mode)
        if hasattr(self, 'lbl_rsd'): self.lbl_rsd.setVisible(is_oversampling_mode)
        if hasattr(self, 'spin_rsd'): self.spin_rsd.setVisible(is_oversampling_mode)
            
        # Update visibility of checkboxes based on active strategy
        is_strict = (strategy == "Strict Pulse Target (Manual)")
        is_adaptive = (strategy == "Adaptive Integer Sync (Auto)")
        if hasattr(self, 'chk_avoid_gaps'):
            self.chk_avoid_gaps.setVisible(is_strict)
        if hasattr(self, 'chk_prefer_exact_pulses'):
            self.chk_prefer_exact_pulses.setVisible(is_adaptive)

        self._on_ui_change()

    def _on_pulses_changed(self, value):
        strategy = self._get(self.cmb_sync_strategy, 'currentText')
        is_oversampling = "Oversampling" in strategy or "Combined" in strategy
        if not is_oversampling:
            if hasattr(self, 'chk_sync_dosage') and self.chk_sync_dosage.isChecked():
                self._block_signals(self.spin_dosage, True)
                self.spin_dosage.setValue(value)
                self._block_signals(self.spin_dosage, False)
        self._on_ui_change()

    def _on_dosage_changed(self, value):
        if hasattr(self, 'chk_sync_dosage') and self.chk_sync_dosage.isChecked():
            self._block_signals(self.chk_sync_dosage, True)
            self.chk_sync_dosage.setChecked(False)
            self._block_signals(self.chk_sync_dosage, False)
        self._on_ui_change()

    def _on_sync_dosage_toggled(self, checked):
        if hasattr(self, 'spin_dosage'):
            if checked:
                self._block_signals(self.spin_dosage, True)
                self.spin_dosage.setValue(self._get(self.spin_pulses, 'value'))
                self._block_signals(self.spin_dosage, False)
        self._on_ui_change()

    def _on_ui_change(self, *args):
        # Refresh visibility of custom fields
        self._update_custom_visibility()

        # Refresh the hardware dashboard and internal state sync
        self._update_hw_summary()
        
        # Update visibility keying off Mode
        self._update_mode_visibility()
        
        # Sync UI precisions across all panels
        self._update_ui_precisions()
        
        # Also force update the status labels in the dialog
        self._update_icp_status()
        self._update_laser_status()
        
        # Update persistent settings (Debounced)
        if hasattr(self, 'save_timer'):
            self.save_timer.start()
        else:
            self.save_persistent_settings()
        
        # Trigger Optimization (Preserve Zoom for parameter tweaks)
        self.run_optimisation(refresh=False, preserve_zoom=True)
        
        # Manually trigger rescale if enabled
        if hasattr(self, 'chk_rescale') and self._get(self.chk_rescale, 'isChecked'):
            self.rescale_to_visible()

    def _update_mode_visibility(self):
        mode = self._get(self.cmb_mode, 'currentText')
        
        # Guard: Don't rebuild layout if mode hasn't changed to avoid focus loss
        if hasattr(self, '_last_mode_vis') and self._last_mode_vis == mode:
            return
        self._last_mode_vis = mode

        is_spot = (mode == 'Spot')
        is_line = (mode == 'Line')
        is_imaging = (mode == 'Imaging')
        
        if not hasattr(self, 'grid_qual'): return

        # 1. Clear Row 2
        # We can't easily "clear" a row in QGridLayout without removing widgets.
        # Strategy: Hide all first, then reposition and show.
        
        # Hide everything first
        if hasattr(self, 'lbl_height'): self.lbl_height.setVisible(False); self.grid_qual.removeWidget(self.lbl_height)
        if hasattr(self, 'spin_height'): self.spin_height.setVisible(False); self.grid_qual.removeWidget(self.spin_height)
        if hasattr(self, 'lbl_width'): self.lbl_width.setVisible(False); self.grid_qual.removeWidget(self.lbl_width)
        if hasattr(self, 'spin_width'): self.spin_width.setVisible(False); self.grid_qual.removeWidget(self.spin_width)
        if hasattr(self, 'lbl_line_len'): self.lbl_line_len.setVisible(False); self.grid_qual.removeWidget(self.lbl_line_len)
        if hasattr(self, 'spin_line_len'): self.spin_line_len.setVisible(False); self.grid_qual.removeWidget(self.spin_line_len)
        if hasattr(self, 'lbl_shots'): self.lbl_shots.setVisible(False); self.grid_qual.removeWidget(self.lbl_shots)
        if hasattr(self, 'spin_shots'): self.spin_shots.setVisible(False); self.grid_qual.removeWidget(self.spin_shots)
        
        # 2. Re-populate based on Mode
        if is_spot:
            # Spot: Shots at (3,0)
            if hasattr(self, 'lbl_shots'): 
                self.grid_qual.addWidget(self.lbl_shots, 3, 0)
                self.lbl_shots.setVisible(True)
            if hasattr(self, 'spin_shots'): 
                self.grid_qual.addWidget(self.spin_shots, 3, 1)
                self.spin_shots.setVisible(True)
                
        elif is_line:
            # Line: Length at (3,0)
            if hasattr(self, 'lbl_line_len'):
                self.grid_qual.addWidget(self.lbl_line_len, 3, 0)
                self.lbl_line_len.setVisible(True)
            if hasattr(self, 'spin_line_len'):
                self.grid_qual.addWidget(self.spin_line_len, 3, 1)
                self.spin_line_len.setVisible(True)
                
        else: # Imaging
            # Imaging: Width at (3,0), Height at (3,2)
            if hasattr(self, 'lbl_width'):
                self.grid_qual.addWidget(self.lbl_width, 3, 0)
                self.lbl_width.setVisible(True)
            if hasattr(self, 'spin_width'):
                self.grid_qual.addWidget(self.spin_width, 3, 1)
                self.spin_width.setVisible(True)
                
            if hasattr(self, 'lbl_height'):
                self.grid_qual.addWidget(self.lbl_height, 3, 3)
                self.lbl_height.setVisible(True)
            if hasattr(self, 'spin_height'):
                self.grid_qual.addWidget(self.spin_height, 3, 4)
                self.spin_height.setVisible(True)

        # 3. Dynamic visibility for Optimised Settings (Results Summary)
        # Hide Speed, Overlap, and Dosage fields in Spot mode
        show_scan_settings = not is_spot
        if hasattr(self, 'lbl_opt_speed_title'): self.lbl_opt_speed_title.setVisible(show_scan_settings)
        if hasattr(self, 'lbl_opt_speed'): self.lbl_opt_speed.setVisible(show_scan_settings)
        
        if hasattr(self, 'lbl_opt_overlap_title'): self.lbl_opt_overlap_title.setVisible(show_scan_settings)
        if hasattr(self, 'lbl_opt_overlap'): self.lbl_opt_overlap.setVisible(show_scan_settings)
        
        if hasattr(self, 'lbl_opt_dosage_title'): self.lbl_opt_dosage_title.setVisible(show_scan_settings)
        if hasattr(self, 'lbl_opt_pulses'): self.lbl_opt_pulses.setVisible(show_scan_settings)



    def _handle_mfr_changed(self, mfr):
        self.cmb_model.blockSignals(True)
        self.cmb_model.clear()
        if mfr in ICP_SPECS:
            self.cmb_model.addItems(sorted(ICP_SPECS[mfr].keys()))
            
            # Smart Persistence: Try specific model for this mfr, then global, then default
            saved_model = self.persistent_settings.get(f"last_model_for_{mfr}")
            if not saved_model:
                saved_model = self.persistent_settings.get('model')
            
            if saved_model in ICP_SPECS[mfr]:
                self.cmb_model.setCurrentText(saved_model)
        self.cmb_model.blockSignals(False)
        model = self._get(self.cmb_model, 'currentText')
        self._handle_model_changed(model)
        self.lbl_result.setText(f"System: {mfr} {model} ({self.icp_tech})")
        


    def _handle_laser_mfr_changed(self, mfr):
        self.cmb_laser_mod.blockSignals(True)
        self.cmb_cell.blockSignals(True)
        
        mfr = self._get(self.cmb_laser_mfr, 'currentText')
        self.cmb_laser_mod.clear()
        self.cmb_cell.clear()
        
        if mfr in LASER_PLATFORMS:
            # Populate sources from association list
            allowed_lasers = LASER_PLATFORMS[mfr]['lasers']
            self.cmb_laser_mod.addItems(sorted(allowed_lasers))
            
            # Smart Persistence for Lasers
            saved_lmod = self.persistent_settings.get(f"last_laser_for_{mfr}")
            if not saved_lmod:
                saved_lmod = self.persistent_settings.get('laser_model')
                
            if saved_lmod in allowed_lasers:
                self.cmb_laser_mod.setCurrentText(saved_lmod)
            
            # Populate cells
            stages = LASER_PLATFORMS[mfr]['stages']
            self.cmb_cell.addItems(sorted(stages))
            
            # Smart Persistence for Cells
            saved_cell = self.persistent_settings.get(f"last_cell_for_{mfr}")
            if not saved_cell:
                saved_cell = self.persistent_settings.get('cell_type')
                
            if saved_cell in stages:
                self.cmb_cell.setCurrentText(saved_cell)
        
        self.cmb_laser_mod.blockSignals(False)
        self.cmb_cell.blockSignals(False)
        
        # Trigger updates
        self._handle_laser_mod_changed(self._get(self.cmb_laser_mod, 'currentText'))
        self._handle_cell_changed(self._get(self.cmb_cell, 'currentText'))
        
    def _handle_laser_mod_changed(self, lmod):
        mfr = self._get(self.cmb_laser_mfr, 'currentText')
        if mfr:
            self.persistent_settings[f"last_laser_for_{mfr}"] = lmod
        self._update_hw_summary()

    def _handle_cell_changed(self, cell):
        mfr = self._get(self.cmb_laser_mfr, 'currentText')
        if mfr:
            self.persistent_settings[f"last_cell_for_{mfr}"] = cell
        self._on_ui_change()

    def _on_opt_spot_changed(self, val):
        self.run_optimisation(refresh=False, fixed_spot_um=val, preserve_zoom=True)

    def _parse_rr_list(self, text):
        try:
            if not text: return None
            # Replace common delimiters with space
            clean = text.replace(',', ' ').replace(';', ' ')
            parts = clean.split()
            vals = sorted([float(p) for p in parts if p.strip()])
            return vals if vals else None
        except:
            return None

    def _get_precision_dps(self):
        """Returns number of decimal places for ICP Dwell based on current hardware precision."""
        if getattr(self, 'precision', 0.1) <= 0.001: return 3
        if getattr(self, 'precision', 0.1) <= 0.01: return 2
        return 1

    def _get_laser_precision_dps(self):
        """Returns number of decimal places for Rep-Rate based on current hardware precision."""
        prec = getattr(self, 'laser_rr_prec', 1.0)
        if prec <= 0.001: return 3
        if prec <= 0.01: return 2
        if prec <= 0.1: return 1
        return 0

    def _get_dwell_unit_info(self):
        """Returns (unit_str, scaler_to_ms, display_dps) based on hardware and user preference."""
        mfr = self._get(self.cmb_mfr, 'currentText')
        pref_unit = self._get(self.cmb_dwell_unit, 'currentText')
        dps = self._get_precision_dps()
        
        if pref_unit == "Auto":
            unit = MFR_DEFAULT_UNITS.get(mfr, "ms")
        else:
            unit = pref_unit
            
        if unit == "s":
            scaler = 1000.0
            disp_dps = dps + 3 # 0.1ms -> 0.0001s
        else:
            scaler = 1.0
            disp_dps = dps
            
        return unit, scaler, disp_dps

    def _update_ui_precisions(self):
        """Updates decimal places for all relevant UI input spinners."""
        rr_dps = self._get_laser_precision_dps()
        icp_dps = self._get_precision_dps()
        
        # Update Input panels
        if hasattr(self, 'spin_init_rr'):
            self.spin_init_rr.setDecimals(rr_dps)
            # Ensure step matches precision
            self.spin_init_rr.setSingleStep(getattr(self, 'laser_rr_prec', 1.0))
            
        # Update Pulse Train Simulator spinners to match Optimiser resolutions
        if hasattr(self, 'spin_pulse_rr'):
            self.spin_pulse_rr.setDecimals(rr_dps)
            self.spin_pulse_rr.setSingleStep(getattr(self, 'laser_rr_prec', 1.0))
        if hasattr(self, 'spin_pulse_at'):
            self.spin_pulse_at.setDecimals(icp_dps)
            self.spin_pulse_at.setSingleStep(getattr(self, 'precision', 0.1))
            
        # Update Target Panel if needed (e.g. if we had RR related there)
        # Note: Pulses and Sigma are usually not bound to hardware RR precision steps. 
        # But we ensure they are at least readable.
        
        # Update custom hardware spinners to match their own precision requirements
        if hasattr(self, 'spin_cust_min'):
            self.spin_cust_min.setDecimals(4) # Custom always high
        if hasattr(self, 'spin_cust_rr_prec'):
            self.spin_cust_rr_prec.setDecimals(3)


    def load_persistent_settings(self):
        IoLog.information(f"iolite Optimiser: Loading settings from {self.settings_json_path}")
        try:
            if os.path.exists(self.settings_json_path):
                with open(self.settings_json_path, 'r') as f:
                    s_data = json.load(f)
                    # IoLog.information("iolite Optimiser: Settings loaded successfully") # Silent
                    return s_data
            else:
                IoLog.information("iolite Optimiser: No settings file found")
        except Exception as e:
            IoLog.error(f"iolite Optimiser: Error loading settings: {e}")
        return {}

    def save_persistent_settings(self):
        try:
            # Update existing persistent settings with current UI state
            self.persistent_settings.update({
                'mfr': self._get(self.cmb_mfr, 'currentText'),
                'model': self._get(self.cmb_model, 'currentText'),
                'dwell_unit_pref': self._get(self.cmb_dwell_unit, 'currentText'),
                'mode': self._get(self.cmb_mode, 'currentText'),
                'laser_mfr': self._get(self.cmb_laser_mfr, 'currentText'),
                'laser_model': self._get(self.cmb_laser_mod, 'currentText'),
                'shots': self._get(self.spin_shots, 'value'),
                'spot_size': self._get(self.spin_spot, 'value'),
                'num_analytes': self._get(self.spin_num_analytes, 'value'),
                'washout': self._get(self.spin_wash, 'value'),
                'init_rr': self._get(self.spin_init_rr, 'value'),
                'scale_signal': self._get(self.chk_scale, 'isChecked'),
                'pulses': self._get(self.spin_pulses, 'value'),
                'target_rsd': self._get(self.spin_rsd, 'value') if hasattr(self, 'spin_rsd') else 5.0,
                'dosage': self._get(self.spin_dosage, 'value'),
                'sync_dosage': self._get(self.chk_sync_dosage, 'isChecked'),
                'target_sigma': self._get(self.spin_sigma, 'value'),
                'min_snr': self._get(self.spin_snr, 'value'),
                'min_duty': self._get(self.spin_duty, 'value'),
                'tof_margin': self._get(self.spin_tof_margin, 'value') if hasattr(self, 'spin_tof_margin') else 10.0,
                'cell_type': self._get(self.cmb_cell, 'currentText'),
                # 'auto_detect': self._get(self.chk_auto, 'isChecked'), # Not saved
                # 'normalise': self._get(self.chk_norm, 'isChecked'), # Not saved
                # 'y_zoom': self._get(self.chk_y_zoom, 'isChecked'), # Not saved
                # Add internal states that are expensive or complex to derive
                'icp_tech': self.icp_tech,
                'min_dwell': self.min_dwell,
                'precision': self.precision,
                'max_rr': self.max_rr,
                'max_speed': self.max_speed,
                'cust_min': self._get(self.spin_cust_min, 'value'),
                'cust_prec': self._get(self.spin_cust_prec, 'value'),
                'cust_type': self._get(self.cmb_cust_type, 'currentText'),
                'cust_dwell_list': self._get(self.edit_cust_dwells, 'text'),
                'cust_rr': self._get(self.spin_cust_rr, 'value'),
                'cust_rr_prec': self._get(self.spin_cust_rr_prec, 'value'),
                'cust_speed': self._get(self.spin_cust_speed, 'value'),
                'cust_rr_type': self._get(self.cmb_cust_rr_type, 'currentText'),
                'cust_rr_list': self._get(self.edit_cust_rr_list, 'text'),
                
                # Plot Options
                'limit_vis': self._get(self.chk_limit_vis, 'isChecked'),
                'max_vis': self._get(self.spin_max_vis, 'value'),
                
                'img_height': self._get(self.spin_height, 'value'),
                'img_width': self._get(self.spin_width, 'value'),
                'line_length': self._get(self.spin_line_len, 'value') if hasattr(self, 'spin_line_len') else 1000.0,
                'avoid_gaps': self._get(self.chk_avoid_gaps, 'isChecked'),
                'prefer_exact_pulses': self._get(self.chk_prefer_exact_pulses, 'isChecked') if hasattr(self, 'chk_prefer_exact_pulses') else False,
                
                # New Persistence Keys
                'theme': self._get(self.combo_theme, 'currentText'),
                'opt_y_zoom': self._get(self.chk_y_zoom, 'isChecked'),
                'normalise': self._get(self.chk_norm, 'isChecked'),
                'cust_rr_list': self._get(self.edit_cust_rr_list, 'text'),
                'img_height': self._get(self.spin_height, 'value') if hasattr(self, 'spin_height') else 1000.0,
                'img_width': self._get(self.spin_width, 'value') if hasattr(self, 'spin_width') else 1000.0,
                'avoid_gaps': self._get(self.chk_avoid_gaps, 'isChecked') if hasattr(self, 'chk_avoid_gaps') else False,
                'spr_time_unit': self._get(self.cmb_spr_unit, 'currentText'),
                'spr_min_distance': self._get(self.spin_spr_dist, 'value'),
                'spr_baseline_window': self._get(self.spin_spr_baseline_window, 'value'),
                'spr_apply_smooth': self._get(self.chk_spr_smooth, 'isChecked'),
                'spr_smooth_window': self._get(self.spin_spr_smooth_window, 'value'),
                'opt_y_zoom': self._get(self.chk_y_zoom, 'isChecked'),
                # 'opt_rescale_y': Ephemeral
                'spr_y_zoom': self._get(self.chk_y_zoom_spr, 'isChecked'),
                'sync_strategy': self._get(self.cmb_sync_strategy, 'currentText'),
            })
            if self.settings_json_path:
                with open(self.settings_json_path, 'w') as f:
                    json.dump(self.persistent_settings, f, indent=4)
                # Log success only once or rarely to avoid spam, but for now we log it:
                # IoLog.information(f"iolite Optimiser: Settings saved to {self.settings_json_path}")
        except Exception as e:
            IoLog.error(f"iolite Optimiser: Error saving settings: {e}")

    def _get_theme_safe_palette(self, is_opt=False):
        """
        Generates a 60-color Mega-Palette (b+c+std) with theme-aware lightness 
        normalization and a fixed Purple (Tab20 Index 8) start.
        
        If is_opt is True, it persists the random shuffle order so that 
        isotope-color assignments remain stable across theme updates.
        """
        # Get raw colors from chained palettes
        mega_raw = []
        for cmap_name in ['tab20b', 'tab20c', 'tab20']:
            mega_raw.extend(list(getattr(plt.cm, cmap_name).colors))
            
        # Specific Anchor: Purple from Tab20 (Index 8)
        anchor_rgb = list(plt.cm.tab20.colors)[8]
        
        # Detect Theme (is it dark?)
        face_col = plt.rcParams['axes.facecolor']
        rgb_face = mcolors.to_rgb(face_col)
        is_dark = (sum(rgb_face)/3.0) < 0.5
        
        # Contrast Bounds
        l_min, l_max = (0.55, 0.85) if is_dark else (0.2, 0.6)
        
        def normalise_rgb(c_rgb):
            h, l, s = colorsys.rgb_to_hls(*c_rgb)
            l = max(l_min, min(l_max, l))
            if is_dark: s = max(0.4, min(1.0, s * 1.1))
            c_new = colorsys.hls_to_rgb(h, l, s)
            return mcolors.to_hex(c_new).lower()

        # Generate normalised anchor
        anchor_hex = normalise_rgb(anchor_rgb)
        
        # For the pool, we want stability if is_opt
        # We store the indices of colors in mega_raw (excluding the anchor)
        anchor_idx = 48 # tab20[8] in the [b, c, std] chain
        
        if is_opt and self._opt_palette_indices is not None:
            # Reuse stored shuffle
            indices = self._opt_palette_indices
        else:
            # Create new shuffle
            indices = [i for i in range(len(mega_raw)) if i != anchor_idx]
            random.shuffle(indices)
            if is_opt:
                self._opt_palette_indices = indices
                
        # Generate final palette
        pool = [normalise_rgb(mega_raw[i]) for i in indices]
        
        # Ensure anchor is unique (unlikely but safe)
        if anchor_hex in pool:
            pool.remove(anchor_hex)
            
        return [anchor_hex] + pool

    def suggest_formats_from_timestamp(self, timestamp_str, attempted_format):
        try:
            import re
            ts_tokens = re.split(r'([^0-9a-zA-Z]+)', timestamp_str)
            ts_vals = [t for i, t in enumerate(ts_tokens) if i % 2 == 0 and t]
            
            if not ts_vals:
                return []
                
            supported_formats = [
                "dd MM yyyy hh mm ss",
                "MM dd yyyy hh mm ss",
                "yyyy MM dd hh mm ss",
                "yyyy dd MM hh mm ss",
                "dd MMMM yyyy hh mm ss",
                "MMMM dd yyyy hh mm ss",
                "dd MMM yyyy hh mm ss",
                "MMM dd yyyy hh mm ss"
            ]
            
            suggestions = []
            for fmt in supported_formats:
                fmt_tokens = re.split(r'([^a-zA-Z]+)', fmt)
                fmt_vals = [t for i, t in enumerate(fmt_tokens) if i % 2 == 0 and t]
                
                if len(fmt_vals) != len(ts_vals):
                    continue
                    
                is_valid = True
                for f_token, t_val in zip(fmt_vals, ts_vals):
                    if f_token == 'yyyy':
                        if not (len(t_val) == 4 and t_val.isdigit()):
                            is_valid = False
                            break
                    elif f_token == 'yy':
                        if not (len(t_val) == 2 and t_val.isdigit()):
                            is_valid = False
                            break
                    elif f_token == 'dd':
                        if not (t_val.isdigit() and 1 <= int(t_val) <= 31):
                            is_valid = False
                            break
                    elif f_token == 'MM':
                        if not (t_val.isdigit() and 1 <= int(t_val) <= 12):
                            is_valid = False
                            break
                    elif f_token == 'MMM':
                        if not (t_val.isalpha() and len(t_val) == 3):
                            is_valid = False
                            break
                    elif f_token == 'MMMM':
                        if not (t_val.isalpha() and len(t_val) >= 3):
                            is_valid = False
                            break
                    elif f_token == 'hh':
                        if not (t_val.isdigit() and 0 <= int(t_val) <= 23):
                            is_valid = False
                            break
                    elif f_token == 'mm' or f_token == 'ss':
                        if not (t_val.isdigit() and 0 <= int(t_val) <= 59):
                            is_valid = False
                            break
                            
                if is_valid:
                    suggestions.append(fmt)
                    
            return suggestions
        except Exception:
            return []

    def check_timestamp_error_in_log(self):
        try:
            import tempfile
            import json
            import re
            temp_dir = tempfile.gettempdir()
            temp_log_path = os.path.normpath(os.path.join(temp_dir, "iolite_optimiser_log.json"))
            
            IoLog.saveLogToFile(temp_log_path)
            
            if not os.path.exists(temp_log_path):
                return None
                
            with open(temp_log_path, 'r', encoding='utf-8') as f:
                log_data = json.load(f)
                
            entries = []
            if isinstance(log_data, dict):
                entries = log_data.get('Log Entries', [])
            elif isinstance(log_data, list):
                entries = log_data
                
            if isinstance(entries, list):
                # Search backwards for the most recent timestamp/csv import failure message
                for entry in reversed(entries):
                    msg_str = ""
                    if isinstance(entry, dict):
                        # Safely search standard fields (Message with capital M is used in iolite's JSON output)
                        msg_str = entry.get("Message", entry.get("message", entry.get("text", entry.get("msg", ""))))
                        if not msg_str:
                            msg_str = " ".join(str(v) for v in entry.values() if isinstance(v, (str, int, float)))
                    elif isinstance(entry, str):
                        msg_str = entry
                        
                    if any(p.lower() in msg_str.lower() for p in ["Could not load", "invalid time stamp", "impossible or invalid"]):
                        # Extract and parse suggestions
                        match = re.search(r'invalid time stamp "([^"]+)" using the date/time format ([^.]+)', msg_str)
                        if match:
                            ts_val = match.group(1)
                            fmt_val = match.group(2).strip().rstrip('.')
                            sugs = self.suggest_formats_from_timestamp(ts_val, fmt_val)
                            if sugs:
                                return msg_str + f"\n\nSuggested format(s) to check in preferences:\n" + "\n".join(f" - {s}" for s in sugs)
                        return msg_str
            
            try:
                os.remove(temp_log_path)
            except Exception:
                pass
        except Exception as le:
            IoLog.warning(f"iolite Optimiser: Failed to parse iolite log: {le}")
        return None

    def import_and_unload_data(self, tab):
        caption = "Import SPR File" if tab == "spr" else "Import Optimisation File"
        file_path = QFileDialog.getOpenFileName(self, caption, "", "Raw data (*.xlsx *.h5 *.io4 *.zip *.prn *.info *.lax *.FIN2 *.REP *.csv *.NC *.CDF *.xml *.vit);;iolite 3 experiment (*.pxp);;All files (*)")
        if not file_path:
            return
        
        norm_file_path = os.path.normpath(file_path)
        
        try:
            res = data.importData(norm_file_path)
            
            # Debug prints to trace import status
            print(f"iolite Optimiser DEBUG: data.importData('{norm_file_path}') returned {res} (type: {type(res)})")
            IoLog.information(f"iolite Optimiser DEBUG: data.importData('{norm_file_path}') returned {res}")
            
            # Check if the file is now in the loaded files list (covers cases where importData returns False but succeeds)
            was_imported = False
            for f in data.importedFiles():
                try:
                    f_path = os.path.normpath(f.filePath()).lower()
                    f_name = os.path.normpath(f.fileName()).lower()
                    target_name = os.path.normpath(os.path.basename(norm_file_path)).lower()
                    if f_path == norm_file_path.lower() or f_name == target_name:
                        was_imported = True
                        break
                except Exception:
                    pass
            
            if res or was_imported:
                self.refresh_data(tab=tab, file_path=norm_file_path)
                try:
                    data.unloadFile(norm_file_path)
                except Exception as ue:
                    IoLog.warning(f"iolite Optimiser: Failed to unload file: {ue}")
            else:
                IoLog.warning(f"iolite Optimiser: Import returned False and not found in imported list for {norm_file_path}")
                log_err = self.check_timestamp_error_in_log()
                msg = "No file was imported.\n\nPlease check if the time stamp format in the file is set correctly in iolite."
                if log_err:
                    msg += f"\n\nDetails from iolite log:\n{log_err}"
                QMessageBox.warning(self, "Import Failed", msg)
        except Exception as e:
            msg = f"Failed to import/load file: {e}"
            print(msg)
            print(traceback.format_exc())
            IoLog.error(f"iolite Optimiser: {msg}")
            log_err = self.check_timestamp_error_in_log()
            if log_err:
                msg += f"\n\nDetails from iolite log:\n{log_err}"
            else:
                msg += "\n\nPlease check if the time stamp format in the file is set correctly in iolite."
            QMessageBox.critical(self, "Import Error", msg)

    def refresh_data(self, tab=None, file_path=None):
        # Do NOT force tab index if None. None means "Refresh All".

        try:
            # Reset Auto-Rescale to ON for new data (Ephemeral)
            if hasattr(self, 'chk_rescale'): self.chk_rescale.setChecked(True)
            if hasattr(self, 'chk_rescale_spr'): self.chk_rescale_spr.setChecked(True)
            
            # --- SPR TAB DATA ---
            if tab == "spr" or tab is None:
                new_df_spr, new_meta_spr = self.get_input_dataframe(tab="spr", file_path=file_path)
                if new_df_spr is not None:
                    # Fully clear previous SPR state to prevent data pollution
                    if hasattr(self, 'spr_all_raw_results'): self.spr_all_raw_results.clear()
                    if hasattr(self, 'spr_excluded_peaks'): self.spr_excluded_peaks.clear()
                    if hasattr(self, 'spr_channel_prominence'): self.spr_channel_prominence.clear()
                    if hasattr(self, 'spr_channel_auto_prom'): self.spr_channel_auto_prom.clear()
                    
                    self.spr_df = new_df_spr
                    self.spr_metadata = new_meta_spr # Store specifically for SPR
                    
                    # Resolve Dwell Times & Units specifically for SPR
                    cols_all = [c for c in new_df_spr.columns if c != 'Time']
                    self.spr_dwells = self.resolve_dwell_times(cols_all, tab="spr") if cols_all else {}
                    
                    if self.spr_df is not None and 'Time' in self.spr_df.columns and len(self.spr_df['Time']) > 0:
                         cols = [c for c in self.spr_df.columns if c != 'Time']
                         if hasattr(self, 'cmb_spr_iso'):
                             self._block_signals(self.cmb_spr_iso, True)
                             self.cmb_spr_iso.clear()
                             self.cmb_spr_iso.addItems(cols)
                             self._block_signals(self.cmb_spr_iso, False)
                             
                             # Force refresh to wipe internal caches in run_spr_analysis
                             self.run_spr_analysis(force_refresh=True)

            # --- OPTIMISER TAB DATA ---
            if tab == "opt" or tab is None:
                new_df_opt, new_meta_opt = self.get_input_dataframe(tab="opt", file_path=file_path)
                if new_df_opt is not None:
                    self.opt_df = new_df_opt
                    self.opt_metadata = new_meta_opt # Store specifically for Optimiser
                    self.channel_metadata = self.opt_metadata # Backwards compat alias
                    
                    # Resolve Dwell Times & Units specifically for Optimiser
                    cols_all = [c for c in new_df_opt.columns if c != 'Time']
                    self.opt_dwells = self.resolve_dwell_times(cols_all) if cols_all else {}
                    self.channel_dwells = self.opt_dwells # Backwards compat alias for calc functions

                    if self.opt_df is not None:
                         if 'Time' in self.opt_df.columns and len(self.opt_df['Time']) > 0:
                             self.t_start = float(self.opt_df['Time'].iloc[0])
                             max_rel_t = float(self.opt_df['Time'].iloc[-1] - self.t_start)
                             
                             # Update SpinBox Ranges (Block signals to prevent inadvertent optimization triggers)
                             for sb in [self.spin_bg_start, self.spin_bg_end, self.spin_sig_start, self.spin_sig_end]:
                                 self._block_signals(sb, True)
                                 sb.setRange(-100, max_rel_t + 100) # Give some buffer
                                 self._block_signals(sb, False)
                             
                             # Calculate Actual AT via linear regression
                             time_vals = self.opt_df['Time'].values
                             slope = np.polyfit(np.arange(len(time_vals)), time_vals, 1)[0]
                             self.detected_at_ms = slope * 1000.0
        
                             # Resolve Dwells & Calculate Overhead
                             if len(self.opt_df.columns) > 1:
                                 self.lbl_result.setText(f"Loaded {len(self.opt_df.columns)-1} channels from iolite.")
                              
                             # Recalculate Overhead using the unified method
                             self._recalculate_overhead()
        
                         else:
                             self.t_start = 0.0
                             
                         # Logic for Auto-Detect vs Manual
                         if self._get(self.chk_auto, 'isChecked'):
                             self.run_auto_detect()
                         else:
                             self.update_plot()
                             # Trigger real-time optimization when data refreshed with manual regions
                             self.run_optimisation(refresh=False)
        except Exception as e:
            msg = f"Refresh Error: {e}"
            print(msg)
            IoLog.error(f"iolite Optimiser: {msg}")

    def resolve_dwell_times(self, channel_names, tab=None):
        local_dwells = {} # Local dict instead of shared self.channel_dwells
        
        try:
            # 1. Build map of existing properties
            ts_map = {ts.name: ts for ts in data.timeSeriesList(data.Input)}
            missing = []
            
            self.channel_is_cps = {} # Map col -> bool (True if CPS, False if Counts)

            for name in channel_names:
                found = False
                if name in ts_map:
                    try:
                        ts = ts_map[name]
                        prop = ts.property("Dwell Time (ms)")
                        if prop and float(prop) > 0:
                            local_dwells[name] = float(prop)
                            found = True
                        
                        # Detect Units
                        units = ts.property("Units")
                        
                        if units and "cps" in str(units).lower():
                            self.channel_is_cps[name] = True
                        else:
                            self.channel_is_cps[name] = False
                            
                    except Exception: pass
                
                if not found:
                    missing.append(name)
                    self.channel_is_cps[name] = False # Default
            
            # 2. If missing, check for preconfigured dwells (from TOF dialog) or show dialog
            
            # Log Unit Detection Summary
            n_cps = sum(1 for name, is_cps in self.channel_is_cps.items() if is_cps)
            n_counts = len(self.channel_is_cps) - n_cps
            IoLog.information(f"iolite Optimiser: Unit Detection - {n_cps} channels detected as CPS, {n_counts} as Counts")
            
            # Update Cached Unit Label
            self.cached_unit_label = "CPS" if n_cps > 0 else "Counts"
            
            if missing:
                if tab == "spr":
                    # For SPR tab, we don't need to prompt for dwell times, just default them to 10.0 ms
                    for m in missing:
                        local_dwells[m] = 10.0
                else:
                    # Check for pre-configured dwells (from TOF/Vitesse dialog)
                    preconfigured = getattr(self, '_preconfigured_dwells', {})
                    
                    still_missing = []
                    for m in missing:
                        if m in preconfigured:
                            local_dwells[m] = preconfigured[m]
                        else:
                            still_missing.append(m)
                    
                    if still_missing:
                        # Show unified config dialog for remaining channels
                        at_val = getattr(self, 'detected_at_ms', None)
                        dlg = DataConfigDialog(still_missing, detected_at=at_val, is_tof=False, parent=self)
                        if dlg.exec_():
                            local_dwells.update(dlg.result_dwells)
                            
                            # SAVE BACK TO IOLITE PROPERTIES
                            try:
                                for name, val in dlg.result_dwells.items():
                                    if name in ts_map:
                                        ts_map[name].setProperty("Dwell Time (ms)", float(val))
                                        IoLog.information(f"iolite Optimiser: Saved dwell {val}ms to channel '{name}'")
                            except Exception as e:
                                IoLog.warning(f"iolite Optimiser: Could not save property: {e}")
                        else:
                            IoLog.warning("iolite Optimiser: Dwell configuration cancelled.")
            
            IoLog.information(f"iolite Optimiser: Resolved dwells for {len(local_dwells)} channels")
            return local_dwells

        except Exception as e:
            IoLog.error(f"iolite Optimiser: Error resolving dwell times: {e}")
            return {}

    # --- PLOT INTERACTION ---
    def _on_plot_resize(self, event):
        """
        Dynamically adjusts margins to maintain constant pixel offsets for labels/legend.
        This ensures the plot scales perfectly with the splitter and window.
        """
        canvas = getattr(event, 'canvas', None)
        if not canvas:
            # Fallback for manual calls
            for c in [getattr(self, 'canvas', None), getattr(self, 'spr_canvas', None)]:
                if c: self._on_plot_resize(type('obj', (object,), {'canvas': c}))
            return

        fig = canvas.figure
        w, h = fig.get_size_inches() * fig.dpi
        if w == 0 or h == 0: return

        # Inch-Based Margins (DPI Aware)
        # Prevents scaling issues on high-DPI screens
        r_in = 0.2
        b_in = 0.55
        
        # Dynamic row-aware top padding (Inches)
        rows = getattr(self, 'last_opt_rows', 1) if canvas == getattr(self, 'canvas', None) else getattr(self, 'last_spr_rows', 1)
        # Base 0.25 in + 0.20 in per row (~24px base + 19px per row at 96 DPI) - Tighter as requested
        t_in = 0.25 + (rows * 0.20)
        
        # Convert to Pixels for margin manager
        r_px = r_in * fig.dpi
        b_px = b_in * fig.dpi
        t_px = t_in * fig.dpi

        # Smart Margins Calculation
        # Instead of fixed l_px, we calculate it!
        self._update_smart_margins(canvas, override_margins={'top': t_px, 'bottom': b_px, 'right': r_px})

    def _update_smart_margins(self, canvas, override_margins=None):
        """
        Calculates and applies the exact left margin needed to fit the Y-axis label and tick numbers.
        Keeps the right edge fixed.
        """
        if not canvas: return
        # Check suppression flag
        if getattr(self, '_suppress_margin_update', False): return

        fig = canvas.figure
        ax = fig.axes[0] if fig.axes else None
        if not ax: return

        # Calculate dimensions first
        w, h = fig.get_size_inches() * fig.dpi
        if w == 0 or h == 0: return

        # Get Renderer
        try:
            renderer = canvas.get_renderer()
        except:
            # If renderer not ready yet (before first draw), use estimate
             renderer = None

        # l_px Logic Refactor (Inch-based)
        anchor_in = 0.3
        anchor_px = anchor_in * fig.dpi
        
        l_px = 55 # Default fallback
        
        if renderer:
            try:
                # FORCE TICK UPDATE
                ax.yaxis.get_tightbbox(renderer)

                # 1. Measure Y-Axis Tick Labels (Max Visible Width)
                tick_width = 0
                ylim = ax.get_ylim()
                y_min, y_max = sorted(ylim)
                
                for tick in ax.yaxis.get_major_ticks():
                    loc = tick.get_loc()
                    if y_min <= loc <= y_max:
                        label = tick.label1 # Left label
                        if label.get_visible():
                            bbox = label.get_window_extent(renderer)
                            if bbox.width > tick_width:
                                tick_width = bbox.width
                
                # 2. Measure Y-Axis Label
                ylabel = ax.yaxis.get_label()
                label_width = ylabel.get_window_extent(renderer).width
                
                # Gap settings
                gap = 2 
                
                # Total Required Left Margin
                # AnchorPX + LabelWidth + Gap + TickWidth + Gap
                calc_l_px = anchor_px + label_width + gap + tick_width + 2
                l_px = max(anchor_px + 20, calc_l_px) 
                
            except Exception as e:
                pass
        
        # RE-ASSERT LABEL ANCHOR (Prevent drift during pan/zoom)
        if ax:
            import matplotlib.transforms as mtransforms
            # 0.3 inches from left edge
            ax.yaxis.set_label_coords(anchor_px/w, 0.5, 
                                      transform=mtransforms.blended_transform_factory(fig.transFigure, ax.transAxes))

        # Apply Margins (w, h already calculated above)
        
        # Use overrides or defaults
        # We need to preserve current top/bottom/right if not provided?
        # But subplots_adjust takes all.
        # We should store them or recalculate them.
        
        # Recalculate defaults if not passed (similar to _on_plot_resize logic)
        if not override_margins:
            rows = getattr(self, 'last_opt_rows', 1) if canvas == getattr(self, 'canvas', None) else getattr(self, 'last_spr_rows', 1)
            # Match tuned inch-based logic
            t_px = (0.25 + (rows * 0.20)) * fig.dpi
            b_px = 0.55 * fig.dpi
            r_px = 0.2 * fig.dpi
        else:
            t_px = override_margins.get('top', 42)
            b_px = override_margins.get('bottom', 50)
            r_px = override_margins.get('right', 20)

        # Apply
        fig.subplots_adjust(
            left = l_px / w,
            right = 1.0 - (r_px / w),
            top = 1.0 - (t_px / h),
            bottom = b_px / h
        )
        
        # Force redraw if needed? No, subplots_adjust usually triggers draw if in event loop
        # But update_geometry might be needed

    def rescale_to_visible(self, rescale_x=False, rescale_y=False, ax=None, canvas=None):
        if ax is None:
            if not hasattr(self, 'figure') or not self.figure.axes: return
            ax = self.figure.axes[0]
        if canvas is None:
            canvas = getattr(self, 'canvas', None)

        # Determine which settings to use
        is_spr = (hasattr(self, 'spr_figure') and ax in self.spr_figure.axes)
        chk_rescale = getattr(self, 'chk_rescale_spr', None) if is_spr else getattr(self, 'chk_rescale', None)
        
        # 1. Rescale X to full data range
        if rescale_x:
            ax.relim()
            ax.autoscale(axis='x', tight=True)
            # Ensure we don't zoom out too far if there's no data
            xlim = ax.get_xlim()
            if xlim[0] == xlim[1]:
                ax.set_xlim(0, 100)
        
        # 2. Rescale Y to visible data only
        # Force if rescale_y is True, OR if Auto-Rescale checkbox is checked
        if rescale_y or (chk_rescale and self._get(chk_rescale, 'isChecked')):
            ymin, ymax = float('inf'), float('-inf')
            found_visible = False
            
            # Lines (Signal)
            if hasattr(ax, 'get_lines'):
                for line in ax.get_lines():
                    if line.get_visible():
                        ydata = line.get_ydata()
                        if ydata is not None and len(ydata) > 0:
                            y_valid = ydata[~np.isnan(ydata)]
                            if len(y_valid) > 0:
                                ymin = min(ymin, np.min(y_valid))
                                ymax = max(ymax, np.max(y_valid))
                                found_visible = True
            
            # If nothing found, check collections (e.g. hlines markers in SPR)
            if not found_visible:
                # Default range if no data
                ax.relim()
                ax.autoscale(axis='y', tight=False)
            else:
                margin = (ymax - ymin) * 0.1 if ymax > ymin else 1.0 # Increased margin to 10%
                if margin == 0: margin = 1.0
                ax.set_ylim(ymin - margin, ymax + margin)
            
        if canvas:
            canvas.draw()

    def _on_opt_rescale_toggled(self, checked):
        if checked:
            # Auto-Rescale ON -> Pan/Zoom OFF
            self._block_signals(self.chk_y_zoom, True)
            self.chk_y_zoom.setChecked(False)
            self._block_signals(self.chk_y_zoom, False)
            
            # Explicitly clear saved zoom limits so rescale works
            self._opt_axis_limits = None
            
            # Trigger Rescale
            if hasattr(self, 'figure') and self.figure.axes:
                self.rescale_to_visible(ax=self.figure.axes[0], canvas=self.canvas)
        self.save_persistent_settings()

    def _on_spr_rescale_toggled(self, checked):
        if checked:
            # Auto-Rescale ON -> Pan/Zoom OFF
            self._block_signals(self.chk_y_zoom_spr, True)
            self.chk_y_zoom_spr.setChecked(False)
            self._block_signals(self.chk_y_zoom_spr, False)
            # Trigger Rescale
            if hasattr(self, 'spr_figure') and self.spr_figure.axes:
                self.rescale_to_visible(ax=self.spr_figure.axes[0], canvas=self.spr_canvas)
        self.save_persistent_settings()

    def _on_opt_y_zoom_toggled(self, checked):
        if checked:
            # Pan/Zoom ON -> Auto-Rescale OFF
            self._block_signals(self.chk_rescale, True)
            self.chk_rescale.setChecked(False)
            self._block_signals(self.chk_rescale, False)
        self.save_persistent_settings()

    def _on_spr_y_zoom_toggled(self, checked):
        if checked:
            # Pan/Zoom ON -> Auto-Rescale OFF
            self._block_signals(self.chk_rescale_spr, True)
            self.chk_rescale_spr.setChecked(False)
            self._block_signals(self.chk_rescale_spr, False)
        self.save_persistent_settings()

    def _on_pulse_rescale_toggled(self, checked):
        if checked:
            self._block_signals(self.chk_y_zoom_pulse, True)
            self.chk_y_zoom_pulse.setChecked(False)
            self._block_signals(self.chk_y_zoom_pulse, False)
            if hasattr(self, 'pulse_figure') and self.pulse_figure.axes:
                self.rescale_to_visible(ax=self.pulse_figure.axes[0], canvas=self.pulse_canvas)
        self.save_persistent_settings()

    def _on_pulse_y_zoom_toggled(self, checked):
        if checked:
            self._block_signals(self.chk_rescale_pulse, True)
            self.chk_rescale_pulse.setChecked(False)
            self._block_signals(self.chk_rescale_pulse, False)
        self.save_persistent_settings()

    def on_zoom(self, event):
        canvas = event.canvas
        fig = canvas.figure
        ax = event.inaxes
        
        # 1. Resolve Target AX (Handle scrolling on labels outside data area)
        if ax is None:
            for a in fig.axes:
                bbox = a.get_window_extent()
                if (bbox.y0 - 50 < event.y < bbox.y1 + 50) and (bbox.x0 - 50 < event.x < bbox.x1 + 50):
                    ax = a
                    break
        
        if ax is None or ax == getattr(self, 'spr_ax_comp', None): return
        
        # 2. Scale factor
        if event.button == 'up': scale_factor = 0.8 # Zoom IN
        elif event.button == 'down': scale_factor = 1.25 # Zoom OUT
        else: return
        
        bbox = ax.get_window_extent()
        
        # 3. Determine Zoom Region
        # Over Y-axis label/spine area (Left)
        over_y_axis = (bbox.y0 <= event.y <= bbox.y1) and (event.x < bbox.x0)
        # Over X-axis label/spine area (Bottom)
        over_x_axis = (bbox.x0 <= event.x <= bbox.x1) and (event.y < bbox.y0)
        # Over the actual plot area
        over_plot = (bbox.x0 <= event.x <= bbox.x1) and (bbox.y0 <= event.y <= bbox.y1)
        
        zoom_x = over_x_axis or over_plot
        zoom_y = over_y_axis # By default if on Y axis
        
        if over_plot:
            # Check if Y zoom is enabled by UI for central plot area
            is_spr = (hasattr(self, 'spr_canvas') and canvas == self.spr_canvas)
            is_pulse = (hasattr(self, 'pulse_canvas') and canvas == self.pulse_canvas)
            if is_spr: chk_y = getattr(self, 'chk_y_zoom_spr', None)
            elif is_pulse: chk_y = getattr(self, 'chk_y_zoom_pulse', None)
            else: chk_y = getattr(self, 'chk_y_zoom', None)
            
            if self._get(chk_y, 'isChecked'):
                zoom_y = True
        
        # 4. Apply Zoom
        if zoom_x:
            cur_xlim = ax.get_xlim()
            cur_xrange = (cur_xlim[1] - cur_xlim[0])
            
            # Use mouse position if in data area, else convert from pixels
            if event.xdata is not None:
                xdata_pivot = event.xdata
            else:
                # Convert pixel X to data X even if outside data range (e.g. over labels)
                try: xdata_pivot = ax.transData.inverted().transform((event.x, event.y))[0]
                except: xdata_pivot = (cur_xlim[0] + cur_xlim[1]) / 2.0
            
            new_width = cur_xrange * scale_factor
            relx = (cur_xlim[1] - xdata_pivot) / cur_xrange
            ax.set_xlim([xdata_pivot - new_width * (1-relx), xdata_pivot + new_width * relx])
            
        if zoom_y:
            # Uncheck Auto-Rescale if manual zooming Y
            is_spr = (hasattr(self, 'spr_canvas') and canvas == self.spr_canvas)
            is_pulse = (hasattr(self, 'pulse_canvas') and canvas == self.pulse_canvas)
            rescale_chk = getattr(self, 'chk_rescale_spr', None) if is_spr else (getattr(self, 'chk_rescale_pulse', None) if is_pulse else getattr(self, 'chk_rescale', None))
            if rescale_chk and self._get(rescale_chk, 'isChecked'):
                 self._block_signals(rescale_chk, True)
                 rescale_chk.setChecked(False)
                 self._block_signals(rescale_chk, False)

            cur_ylim = ax.get_ylim()
            cur_yrange = (cur_ylim[1] - cur_ylim[0])
            
            # Use mouse position if in data area, else convert from pixels
            if event.ydata is not None:
                ydata_pivot = event.ydata
            else:
                # Convert pixel Y to data Y even if outside data range (e.g. over labels)
                try: ydata_pivot = ax.transData.inverted().transform((event.x, event.y))[1]
                except: ydata_pivot = (cur_ylim[0] + cur_ylim[1]) / 2.0
                
            new_height = cur_yrange * scale_factor
            rely = (cur_ylim[1] - ydata_pivot) / cur_yrange
            ax.set_ylim([ydata_pivot - new_height * (1-rely), ydata_pivot + new_height * rely])
            
        canvas.draw()

    def on_press(self, event):
        canvas = event.canvas
        fig = canvas.figure
        ax = event.inaxes
        
        # Resolve Target AX (Handle clicking on labels outside data area)
        if ax is None:
            for a in fig.axes:
                bbox = a.get_window_extent()
                # Use a buffer for easier clicking
                if (bbox.y0 - 50 < event.y < bbox.y1 + 50) and (bbox.x0 - 50 < event.x < bbox.x1 + 50):
                    ax = a
                    break

        if ax is None or ax == getattr(self, 'spr_ax_comp', None): return
        
        if event.dblclick:
            # RESET VIEW (X and Y)
            self.rescale_to_visible(rescale_x=True, rescale_y=True, ax=ax, canvas=canvas)
            return
            
        if event.button == 1: # Left Click
            # REDUNDANT SHIFT DETECTION
            mod = 0
            try: mod = QApplication.instance().keyboardModifiers()
            except: 
                try: mod = QApplication.keyboardModifiers()
                except: pass
            
            is_shift = bool(mod & Qt.ShiftModifier)
            if not is_shift:
                is_shift = str(getattr(event, 'key', '')).lower() in ('shift', 's')
            
            bbox = ax.get_window_extent()
            
            # Determine Interaction Region
            # Over Y-axis label/spine area (Left)
            over_y_axis = (bbox.y0 <= event.y <= bbox.y1) and (event.x < bbox.x0)
            # Over X-axis label/spine area (Bottom)
            over_x_axis = (bbox.x0 <= event.x <= bbox.x1) and (event.y < bbox.y0)
            # Over the actual plot area
            over_plot = (bbox.x0 <= event.x <= bbox.x1) and (bbox.y0 <= event.y <= bbox.y1)
            
            self.press_region = 'y_axis' if over_y_axis else ('x_axis' if over_x_axis else 'plot')

            # 1. Handle Region Edge Picking (Opt and SPR Plots)
            is_opt_canvas = canvas == getattr(self, 'canvas', None)
            is_spr_canvas = canvas == getattr(self, 'spr_canvas', None)
            if is_shift and (is_opt_canvas or is_spr_canvas):
                if event.xdata is not None:
                    xlim = ax.get_xlim()
                    # Permissive threshold: 1% of plot width
                    threshold = (xlim[1] - xlim[0]) * 0.01
                    t_shift = getattr(self, 't_start', 0)
                    
                    found_edge = None
                    if is_opt_canvas:
                        if hasattr(self, 'bg_times') and self.bg_times:
                            s_bg, e_bg = self.bg_times[0] - t_shift, self.bg_times[1] - t_shift
                            if abs(event.xdata - s_bg) < threshold: found_edge = 'bg_start'
                            elif abs(event.xdata - e_bg) < threshold: found_edge = 'bg_end'
                            
                        if not found_edge and hasattr(self, 'sig_times') and self.sig_times:
                            s_sig, e_sig = self.sig_times[0] - t_shift, self.sig_times[1] - t_shift
                            if abs(event.xdata - s_sig) < threshold: found_edge = 'sig_start'
                            elif abs(event.xdata - e_sig) < threshold: found_edge = 'sig_end'
                    elif is_spr_canvas:
                        if self._get(getattr(self, 'chk_spr_bg_sub', None), 'isChecked'):
                            s_bg = self._get(getattr(self, 'spin_spr_bg_start', None), 'value')
                            e_bg = self._get(getattr(self, 'spin_spr_bg_end', None), 'value')
                            if s_bg is not None and e_bg is not None:
                                if abs(event.xdata - s_bg) < threshold: found_edge = 'spr_bg_start'
                                elif abs(event.xdata - e_bg) < threshold: found_edge = 'spr_bg_end'
                        
                    if found_edge:
                        self.dragged_edge = found_edge
                        if found_edge.startswith('spr_'):
                            if hasattr(self, 'chk_spr_auto_bg'):
                                self._block_signals(self.chk_spr_auto_bg, True)
                                self.chk_spr_auto_bg.setChecked(False)
                                self._block_signals(self.chk_spr_auto_bg, False)
                                try:
                                    self.spin_spr_bg_start.setEnabled(True)
                                    self.spin_spr_bg_end.setEnabled(True)
                                except: pass
                        else:
                            if hasattr(self, 'chk_auto'):
                                self._block_signals(self.chk_auto, True)
                                self.chk_auto.setChecked(False)
                                self._block_signals(self.chk_auto, False)
                        return # Edge dragging takes precedence
            
            # 2. Initiate Panning (Standard or SPR) - STORE STARTING STATE
            self.plot_panning = True
            self.press_x_pix = event.x
            self.press_y_pix = event.y
            self.start_xlim = ax.get_xlim()
            self.start_ylim = ax.get_ylim()
            self.press_ax = ax
            # Capture data width for each axis once at start
            self.start_dx_per_pix = (self.start_xlim[1] - self.start_xlim[0]) / (ax.bbox.width or 1)
            self.start_dy_per_pix = (self.start_ylim[1] - self.start_ylim[0]) / (ax.bbox.height or 1)

    def on_drag(self, event):
        # Allow panning/dragging to continue even if mouse moves outside the axis area
        is_panning = getattr(self, 'plot_panning', False)
        is_edge_drag = bool(getattr(self, 'dragged_edge', None))
        
        if event.inaxes is None and not is_panning and not is_edge_drag:
            return
            
        ax = event.inaxes if event.inaxes else getattr(self, 'press_ax', None)
        if ax is None: return

        # 1. Handle Edge Dragging (Priority)
        if self.dragged_edge:
            # For edge dragging, we need xdata. If it's None, try to compute it from pixels.
            xdata = event.xdata
            if xdata is None:
                try: xdata = ax.transData.inverted().transform((event.x, event.y))[0]
                except: return
            
            val = max(0, xdata)
            t_shift = getattr(self, 't_start', 0)
            
            # Determine which spinbox and which internal times to update
            sb = None
            if self.dragged_edge == 'bg_start':
                sb = self.spin_bg_start
                self.bg_times = (val + t_shift, self.bg_times[1])
            elif self.dragged_edge == 'bg_end':
                sb = self.spin_bg_end
                self.bg_times = (self.bg_times[0], val + t_shift)
            elif self.dragged_edge == 'sig_start':
                sb = self.spin_sig_start
                self.sig_times = (val + t_shift, self.sig_times[1])
            elif self.dragged_edge == 'sig_end':
                sb = self.spin_sig_end
                self.sig_times = (self.sig_times[0], val + t_shift)
            elif self.dragged_edge == 'spr_bg_start':
                sb = self.spin_spr_bg_start
            elif self.dragged_edge == 'spr_bg_end':
                sb = self.spin_spr_bg_end
                
            if sb:
                self._block_signals(sb, True)
                sb.setValue(val)
                self._block_signals(sb, False)
            
            # Capture current limits to preserve zoom during full redraw
            if event.inaxes:
                self._opt_axis_limits = (event.inaxes.get_xlim(), event.inaxes.get_ylim())
            
            # Update plot visuals ONLY (no optimization)
            if self.dragged_edge.startswith('spr_'):
                self.run_spr_analysis(rescale=False)
            else:
                self.update_plot(preserve_zoom=True)
            return

        # 2. Handle Panning (Fixed-reference for maximum stability)
        if self.plot_panning and hasattr(self, 'press_x_pix') and hasattr(self, 'press_ax'):
            ax = self.press_ax
            canvas = event.canvas
            if ax == getattr(self, 'spr_ax_comp', None): return
            
            # Calculate pixel displacement from START
            dx_pix = event.x - self.press_x_pix
            dy_pix = event.y - self.press_y_pix
            
            # Convert pixel displacement to data units using scaling factor from START
            dx_data = dx_pix * getattr(self, 'start_dx_per_pix', 0)
            dy_data = dy_pix * getattr(self, 'start_dy_per_pix', 0)
            
            region = getattr(self, 'press_region', 'plot')
            
            # Handle X-Panning
            if region in ('plot', 'x_axis'):
                if hasattr(self, 'start_xlim'):
                    ax.set_xlim(self.start_xlim[0] - dx_data, self.start_xlim[1] - dx_data)
            
            # Handle Y-Panning
            is_spr = (hasattr(self, 'spr_canvas') and canvas == self.spr_canvas)
            is_pulse = (hasattr(self, 'pulse_canvas') and canvas == self.pulse_canvas)
            
            if is_spr: chk_y = getattr(self, 'chk_y_zoom_spr', None)
            elif is_pulse: chk_y = getattr(self, 'chk_y_zoom_pulse', None)
            else: chk_y = getattr(self, 'chk_y_zoom', None)
            
            # Pan Y if: 1. Region is Y-axis OR 2. Region is plot AND checkbox is checked
            pan_y = (region == 'y_axis') or (region == 'plot' and self._get(chk_y, 'isChecked'))
            
            if pan_y and hasattr(self, 'start_ylim'):
                # Uncheck Auto-Rescale if manual panning Y
                rescale_chk = getattr(self, 'chk_rescale_spr', None) if is_spr else (getattr(self, 'chk_rescale_pulse', None) if is_pulse else getattr(self, 'chk_rescale', None))
                if rescale_chk and self._get(rescale_chk, 'isChecked'):
                     self._block_signals(rescale_chk, True)
                     rescale_chk.setChecked(False)
                     self._block_signals(rescale_chk, False)
                ax.set_ylim(self.start_ylim[0] - dy_data, self.start_ylim[1] - dy_data)

            canvas.draw()
            return

        # 3. Hover Cursor Feedback
        mod = 0
        try: mod = QApplication.instance().keyboardModifiers()
        except: 
            try: mod = QApplication.keyboardModifiers()
            except: pass
        is_shift = bool(mod & Qt.ShiftModifier) or str(getattr(event, 'key', '')).lower() in ('shift', 's')
                   
        is_opt_canvas = event.canvas == getattr(self, 'canvas', None)
        is_spr_canvas = event.canvas == getattr(self, 'spr_canvas', None)
        if is_shift and (is_opt_canvas or is_spr_canvas) and event.inaxes:
            xlim = event.inaxes.get_xlim()
            threshold = (xlim[1] - xlim[0]) * 0.01
            t_shift = getattr(self, 't_start', 0)
            
            near_edge = False
            if is_opt_canvas:
                for times in [self.bg_times, self.sig_times]:
                    if times:
                        if abs(event.xdata - (times[0]-t_shift)) < threshold or \
                           abs(event.xdata - (times[1]-t_shift)) < threshold:
                            near_edge = True
                            break
            elif is_spr_canvas:
                if self._get(getattr(self, 'chk_spr_bg_sub', None), 'isChecked'):
                    s_bg = self._get(getattr(self, 'spin_spr_bg_start', None), 'value')
                    e_bg = self._get(getattr(self, 'spin_spr_bg_end', None), 'value')
                    if s_bg is not None and e_bg is not None:
                        if abs(event.xdata - s_bg) < threshold or abs(event.xdata - e_bg) < threshold:
                            near_edge = True
            
            if near_edge:
                event.canvas.setCursor(Qt.SizeHorCursor)
            else:
                event.canvas.setCursor(Qt.ArrowCursor)
        elif event.canvas:
            event.canvas.setCursor(Qt.ArrowCursor)
        else:
            if getattr(event, 'canvas', None):
                event.canvas.setCursor(Qt.ArrowCursor)

    def on_release(self, event):
        was_dragging = self.dragged_edge is not None
        edged = self.dragged_edge
        self.plot_panning = False
        self.press_x_pix = None
        self.press_y_pix = None
        self.dragged_edge = None
        if event.canvas:
            event.canvas.setCursor(Qt.ArrowCursor)
            
        if was_dragging:
            # Capture limits to prevent rescale on release
            ax = getattr(event, 'inaxes', None)
            # Fallback if release happened outside axes
            if ax is None and hasattr(self, 'figure') and len(self.figure.axes) > 0:
                 ax = self.figure.axes[0]
            
            if ax:
                 self._opt_axis_limits = (ax.get_xlim(), ax.get_ylim())

            if edged and edged.startswith('spr_'):
                self.run_spr_analysis(rescale=False)
            else:
                # Trigger final optimization and sync
                self.on_region_edited()

    def on_pick(self, event):
        try:
            artist = event.artist
            mouse = event.mouseevent
            # IGNORE SCROLL WHEEL EVENTS (Zooming should take precedence)
            if mouse.button in ('up', 'down'):
                return
            
            canvas = mouse.canvas
            if canvas is None: return

            # 1. Identify context (Main or SPR)
            is_spr = (hasattr(self, 'spr_canvas') and canvas == self.spr_canvas)
            fig = getattr(self, 'spr_figure', None) if is_spr else getattr(self, 'figure', None)
            if fig is None or not fig.axes: return
            
            ax = fig.axes[0]
            l_frame = getattr(self, 'spr_legend_frame', None) if is_spr else getattr(self, 'legend_frame', None)
            l_map = getattr(self, 'spr_map_legend_to_line', {}) if is_spr else getattr(self, 'map_legend_to_line', {})
            l_chk_rescale = getattr(self, 'chk_rescale_spr', None) if is_spr else getattr(self, 'chk_rescale', None)

            # Helper to get all artists (lines + collections)
            all_artists = list(ax.lines) + list(ax.collections)

            # 2. Restore All (Double click on Legend Background/Frame)
            if l_frame is not None and artist == l_frame:
                if mouse.dblclick:
                    # Check if we hit an item anyway (prevents overriding the Isolate action)
                    hit_item = False
                    for item in l_map.keys():
                        if item.contains(mouse)[0]:
                            hit_item = True
                            break
                    
                    if not hit_item:
                        for art in all_artists:
                            art.set_visible(True)
                        for item in l_map.keys():
                            item.set_alpha(1.0)
                            # Force visual refresh by toggling visibility
                            item.set_visible(False)
                            item.set_visible(True)
                        if l_chk_rescale and self._get(l_chk_rescale, 'isChecked'):
                            self.rescale_to_visible(ax=ax, canvas=canvas)
                        canvas.draw()
                return

            # 3. Peak Picking (SPR exclusion)
            if is_spr and hasattr(self, 'spr_peak_label_map') and artist in self.spr_peak_label_map:
                p_idx = self.spr_peak_label_map[artist]
                iso = self._get(self.cmb_spr_iso, 'currentText')
                if iso not in self.spr_excluded_peaks:
                    self.spr_excluded_peaks[iso] = set()
                
                if p_idx in self.spr_excluded_peaks[iso]:
                    self.spr_excluded_peaks[iso].remove(p_idx)
                else:
                    self.spr_excluded_peaks[iso].add(p_idx)
                
                # Refresh analysis to update UI/Stats (prevent jump)
                self.run_spr_analysis(rescale=False)
                return

            # 4. Handle Legend Item (Line or Text)
            if artist not in l_map:
                return
            
            target_obj, is_hidable = l_map[artist]
            
            # Debounce: If this is the same mouse event and same target as last time, ignore it.
            # (Happens when clicking handle and text simultaneously)
            last_event = getattr(self, '_last_pick_mouse_event', None)
            last_target = getattr(self, '_last_pick_target', None)
            
            # Check for object identity equality of the mouse event (robust enough for mpl events)
            if last_event == mouse and last_target == target_obj:
                return
            
            self._last_pick_mouse_event = mouse
            self._last_pick_target = target_obj
            
            target_list = target_obj if isinstance(target_obj, list) else [target_obj]

            # Check if this is a potential double-click (second click on same target within timer window)
            pending = getattr(self, '_pending_legend_click', None)
            if pending is not None:
                pending_target = pending[1]
                # If same target, this is a double-click - execute ISOLATE
                if set(pending_target) == set(target_list):
                    self._legend_click_timer.stop()
                    self._pending_legend_click = None
                    self._execute_legend_isolate(target_list, l_map, canvas, ax, l_chk_rescale)
                    return
            
            # Not a double-click (yet) - start timer for single-click
            if not is_hidable:
                return  # Static items cannot be toggled
            
            self._pending_legend_click = (target_obj, target_list, canvas, ax, l_map, l_chk_rescale, is_hidable)
            self._legend_click_timer.start()
            return  # Don't process immediately, wait for timer

        except Exception as e:
            IoLog.error(f"Pick Error: {e}")

    def _execute_legend_single_click(self):
        """Execute toggle visibility for a single legend click (timer fired, no double-click)."""
        try:
            if self._pending_legend_click is None:
                return
            
            target_obj, target_list, canvas, ax, l_map, l_chk_rescale, is_hidable = self._pending_legend_click
            self._pending_legend_click = None
            
            if not is_hidable:
                return
            
            # TOGGLE visibility
            new_vis = not target_list[0].get_visible()
            for t in target_list:
                t.set_visible(new_vis)
            
            # Sync alpha for all legend artists mapping to THIS underlying line
            alpha = 1.0 if new_vis else 0.2
            for m_artist, m_spec in l_map.items():
                m_obj = m_spec[0]
                m_targets = m_obj if isinstance(m_obj, list) else [m_obj]
                if set(m_targets) == set(target_list):
                    m_artist.set_alpha(alpha)
                    # Force visual refresh by toggling visibility
                    m_artist.set_visible(False)
                    m_artist.set_visible(True)
            
            # Repaint and rescale
            if l_chk_rescale and self._get(l_chk_rescale, 'isChecked'):
                self.rescale_to_visible(ax=ax, canvas=canvas)
            
            canvas.draw()
            canvas.flush_events()
            if hasattr(canvas, 'repaint'):
                canvas.repaint()
                
        except Exception as e:
            IoLog.error(f"Legend Single Click Error: {e}")

    def _execute_legend_isolate(self, target_list, l_map, canvas, ax, l_chk_rescale):
        """Execute isolate (show only this channel) for a double-click."""
        try:
            # Collect lines in l_map to manage only legend-controlled lines
            manageable_hidable = []
            manageable_static = []
            for m_artist, m_spec in l_map.items():
                m_obj, m_hidable = m_spec
                m_targets = m_obj if isinstance(m_obj, list) else [m_obj]
                for t in m_targets:
                    if m_hidable:
                        if t not in manageable_hidable: manageable_hidable.append(t)
                    else:
                        if t not in manageable_static: manageable_static.append(t)
            
            # Apply Visibility
            for line in manageable_hidable:
                line.set_visible(line in target_list)
            for line in manageable_static:
                line.set_visible(True)
            
            # Sync Legend Alphas
            for m_artist, m_spec in l_map.items():
                m_obj, m_hidable = m_spec
                m_targets = m_obj if isinstance(m_obj, list) else [m_obj]
                is_active = any(t.get_visible() for t in m_targets)
                m_artist.set_alpha(1.0 if is_active else 0.2)
                # Force visual refresh by toggling visibility
                m_artist.set_visible(False)
                m_artist.set_visible(True)
            
            # Repaint and rescale
            if l_chk_rescale and self._get(l_chk_rescale, 'isChecked'):
                self.rescale_to_visible(ax=ax, canvas=canvas)
            
            canvas.draw()
            canvas.flush_events()
            if hasattr(canvas, 'repaint'):
                canvas.repaint()
                
        except Exception as e:
            IoLog.error(f"Legend Isolate Error: {e}")

    def update_plot(self, df=None, preserve_zoom=False, preserve_x_only=False):
        try:
            target_df = df if df is not None else self.opt_df
            if target_df is None: return

            # Performance: Enable path simplification for large datasets
            plt.rcParams['path.simplify'] = True
            plt.rcParams['path.simplify_threshold'] = 1.0  # Aggressive simplification

            # Reuse existing axes if possible to prevent layout shrinking
            # ADAPTING TO MATCH SPR LOGIC FOR STABILITY:
            # We explicitly clear the figure to ensure a clean slate and consistent callback connections
            # This prevents callback accumulation and ensures 'connected on creation' logic holds true.

            # Auto-Capture Limits if preserving zoom but not externally set (Robustness Fix)
            if (preserve_zoom or preserve_x_only) and self._opt_axis_limits is None:
                if hasattr(self, 'figure') and self.figure.axes:
                    ax = self.figure.axes[0]
                    if preserve_x_only:
                        # Only save X limits, let Y autoscale
                        self._opt_axis_limits = (ax.get_xlim(), None)
                    else:
                        self._opt_axis_limits = (ax.get_xlim(), ax.get_ylim())

            # Capture current visibility states to preserve user selections
            saved_visibility = {}
            if (preserve_zoom or preserve_x_only) and hasattr(self, 'figure') and self.figure.axes:
                ax = self.figure.axes[0]
                for line in ax.lines:
                    label = line.get_label()
                    if label and not label.startswith('_'):
                        saved_visibility[label] = line.get_visible()

            self.figure.clear()
            
            # Explicitly set figure facecolor from current rcParams
            self.figure.patch.set_facecolor(plt.rcParams['figure.facecolor'])
            
            ax = self.figure.add_subplot(111)
            # Connect dynamic margin update on zoom/pan (ONCE on creation)
            ax.callbacks.connect('ylim_changed', lambda event: self._update_smart_margins(self.canvas))
            
            if (preserve_zoom or preserve_x_only) and self._opt_axis_limits:
                 # We will restore these later
                 pass
            else:
                 self._opt_axis_limits = None
            
            # Suppress smart margin updates during data loading/formatting to prevent "shudder"
            self._suppress_margin_update = True
            
            # Plot Logic
            lines = []
            numeric_cols = target_df.select_dtypes(include=[np.number]).columns
            # Exclude 'Time' if present
            numeric_cols = [c for c in numeric_cols if c != 'Time']
            
            if len(numeric_cols) > 0:
                # --- THEME-SAFE RANDOMIZED MEGA-PALETTE ---
                # PERSISTENT for Optimisation Tab
                safe_palette = self._get_theme_safe_palette(is_opt=True)
                ax.set_prop_cycle(cycler(color=safe_palette))

                # Time Zeroing
                t_orig = target_df['Time'].values
                t_zeroed = t_orig - t_orig[0] if len(t_orig) > 0 else t_orig
                t_shift = t_orig[0] if len(t_orig) > 0 else 0
                
                # Update Spin Box Ranges
                if len(t_zeroed) > 0:
                    t_max = float(np.nanmax(t_zeroed))
                    for sb in [self.spin_bg_start, self.spin_bg_end, self.spin_sig_start, self.spin_sig_end]:
                        sb.setRange(0, t_max)
                
                # Normalization
                normalise = self._get(self.chk_norm, 'isChecked')
                
                plot_df = target_df.copy()
                
                if normalise:
                    for c in numeric_cols:
                        mn, mx = plot_df[c].min(), plot_df[c].max()
                        if mx > mn: plot_df[c] = (plot_df[c]-mn)/(mx-mn)
                        else: plot_df[c] = 0
                
                # Plot individual channels
                max_vis = self._get(self.spin_max_vis, 'value')
                limit_enabled = self._get(self.chk_limit_vis, 'isChecked')
                
                for i, col in enumerate(numeric_cols):
                    # Create line as visible first (affects legend handle creation)
                    l, = ax.plot(t_zeroed, plot_df[col], alpha=0.8, label=col)
                    
                    # Restore visibility from saved state if available
                    if col in saved_visibility:
                        l.set_visible(saved_visibility[col])
                    # Otherwise apply limit if enabled
                    elif limit_enabled and i >= max_vis:
                        l.set_visible(False)
                        
                    lines.append(l)
                
                # Interactive Legend
                # Constrained Layout handles resizing automatically
                
                
                # Custom Legend Sorting (Row-First filling simulation)
                n_items = len(lines)
                if n_items > 0:
                    ncols = min(10, n_items)
                    nrows = math.ceil(n_items / ncols)
                    self.last_opt_rows = nrows
                    self._on_plot_resize(type('obj', (object,), {'canvas': self.canvas}))
                    
                    grid_size = nrows * ncols
                    ordered_handles = [None] * grid_size
                    ordered_labels = [None] * grid_size
                    
                    for i, (h, lbl) in enumerate(zip(lines, numeric_cols)):
                        # Visual r, c (Row Major)
                        r = i // ncols
                        c = i % ncols
                        
                        # MPL linear index for Col Major (c * nrows + r)
                        mpl_idx = c * nrows + r
                        
                        if mpl_idx < grid_size:
                            ordered_handles[mpl_idx] = h
                            ordered_labels[mpl_idx] = lbl
                            
                    final_handles = []
                    final_labels = []
                    
                    for h, l in zip(ordered_handles, ordered_labels):
                        if h is None:
                            # Invisible proxy for gap (use simple Line2D)
                            final_handles.append(Line2D([0], [0], visible=False))
                            final_labels.append("")
                        else:
                            final_handles.append(h)
                            final_labels.append(l)

                    import matplotlib.transforms as mtransforms
                    leg = ax.legend(final_handles, final_labels, loc='lower center', 
                                    bbox_to_anchor=(0.5, 1.0), 
                                    bbox_transform=mtransforms.blended_transform_factory(ax.figure.transFigure, ax.transAxes),
                                    borderaxespad=0.5, ncol=ncols, frameon=True, fontsize='medium', 
                                    handlelength=1.5, handletextpad=0.7, columnspacing=1.5)
                else:
                    import matplotlib.transforms as mtransforms
                    leg = ax.legend(loc='lower center', 
                                    bbox_to_anchor=(0.5, 1.0), 
                                    bbox_transform=mtransforms.blended_transform_factory(ax.figure.transFigure, ax.transAxes),
                                    borderaxespad=0.5, ncol=min(10, len(lines)), frameon=True, fontsize='medium', 
                                    handlelength=1.5, handletextpad=0.7, columnspacing=1.5)

                leg.get_frame().set_alpha(0.0) # Transparent frame
                leg.get_frame().set_picker(5)
                self.legend_frame = leg.get_frame()
                
                self.map_legend_to_line = {}
                
                # Use final_handles to safely map back to original lines (skipping proxies)
                # Note: leg.get_lines() returns new artists created by legend, corresponding 1:1 with final_handles
                if n_items > 0:
                    for legline, legtext, source_h in zip(leg.get_lines(), leg.get_texts(), final_handles):
                        # Detect proxy by checking if source_h is in our original lines list
                        if source_h not in lines:
                            legline.set_visible(False)
                            legtext.set_visible(False)
                            continue

                        legline.set_picker(5)
                        legtext.set_picker(5)
                        self.map_legend_to_line[legline] = (source_h, True)
                        self.map_legend_to_line[legtext] = (source_h, True)
                        
                        # Sync initial visibility
                        vis = source_h.get_visible()
                        legline.set_alpha(1.0 if vis else 0.2)
                        legtext.set_alpha(1.0 if vis else 0.2)
                
                # Shaded regions (Adjusted for Time 0)
                # Shaded regions (Adjusted for Time 0)
                # Theme-aware colors and alpha (Visible Grey for Background, Teal/Turquoise for Signal)
                if hasattr(self, 'cached_is_dark'):
                    is_dark = self.cached_is_dark
                else:
                    is_dark = plt.rcParams['axes.facecolor'] not in ['white', '#ffffff', 'w']
                bg_col = '#808080' if is_dark else '#9e9e9e' # Slate/Medium Grey (more visible)
                sig_col = '#00ced1' if is_dark else '#008b8b' # Vibrant DarkTurquoise/Teal (Blue-Green)
                v_alpha = 0.35 if is_dark else 0.25
                
                # Force visibility (bring to front if needed, but standard z-order is fine)
                
                if self.bg_times:
                    s, e = self.bg_times
                    ax.axvspan(s - t_shift, e - t_shift, color=bg_col, alpha=v_alpha, label='_nolegend_')
                if self.sig_times:
                    s, e = self.sig_times
                    ax.axvspan(s - t_shift, e - t_shift, color=sig_col, alpha=v_alpha, label='_nolegend_')
                    
                # Use centralized unit
                unit = getattr(self, 'cached_unit_label', "Counts")

                ax.set_xlabel("Time (s)")
                ax.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
                ax.ticklabel_format(useOffset=False, axis='x')

                # Anchored Y-Axis Label (Fixed 0.3 inches from left edge)
                import matplotlib.transforms as mtransforms
                ax.set_ylabel("Norm. Intensity" if normalise else f"Intensity ({unit})")
                w_in = self.figure.get_size_inches()[0]
                ax.yaxis.set_label_coords(0.3/w_in, 0.5, 
                                          transform=mtransforms.blended_transform_factory(self.figure.transFigure, ax.transAxes))
                
                
                # Connect dynamic margin update on zoom/pan
                # MOVED TO AXIS CREATION TO PREVENT CALLBACK ACCUMULATION
                # ax.callbacks.connect('ylim_changed', lambda event: self._update_smart_margins(self.canvas))
                
                # Initial Smart Margin Calc (Post-draw simulation)
                # We defer this slightly or call it if renderer is ready
                # self._update_smart_margins(self.canvas) 

                # Apply SI Formatting to Y-Axis (if not normalised)
                if not normalise:
                    # EngFormatter with sep=" " (space before unit)
                    ax.yaxis.set_major_formatter(EngFormatter(sep=" "))
                
                # Force Colors for Labels/Ticks (Fix Dark Mode Visibility)
                # We specifically grab the 'text.color' or 'axes.labelcolor' from current params
                # which 'style.use' sets, but sometimes gets lost in resize/clear.
                fg_color = plt.rcParams['axes.labelcolor']
                ax.xaxis.label.set_color(fg_color)
                ax.yaxis.label.set_color(fg_color)
                ax.tick_params(axis='x', colors=fg_color)
                ax.tick_params(axis='y', colors=fg_color)
                ax.title.set_color(fg_color)
                for spine in ax.spines.values():
                    spine.set_edgecolor(fg_color)

                # Force tight borders
                ax.margins(x=0)
                
                # Consistent layout management (Match SPR plot stability)
                # (Removed constrained layout pad setting to preserve fixed subplots_adjust)
                
                # Consistent with fixed subplots_adjust, we do not call execute_constrained_layout

                # Use unified rescaling logic if enabled to ensure consistent 10% margins from first draw
                chk_rescale = getattr(self, 'chk_rescale', None)
                
                # --- APPLY THEME-STABLE SCALING ---
                if (preserve_zoom or preserve_x_only) and self._opt_axis_limits:
                    ax.set_xlim(self._opt_axis_limits[0])
                    if self._opt_axis_limits[1] is not None:
                        ax.set_ylim(self._opt_axis_limits[1])
                    else:
                        # Y limits were not saved (preserve_x_only) - autoscale Y
                        self.rescale_to_visible(rescale_x=False, rescale_y=True, ax=ax, canvas=None)
                    self._opt_axis_limits = None # Done
                else:
                    # FRESH LOAD, RECALC, or USER INTERACTION -> Rescale EVERYTHING (X and Y)
                    self._opt_axis_limits = None
                    
                    # Explicitly set X limits to ensure we don't stay in (0,1) default
                    if len(t_zeroed) > 0:
                        ax.set_xlim(0, np.nanmax(t_zeroed))
                    
                    
                    # Use robust Y-rescaling (Defer draw by passing canvas=None)
                    self.rescale_to_visible(rescale_x=True, rescale_y=True, ax=ax, canvas=None)
                
                # FORCE SMART MARGIN UPDATE (Ensure Formatter is accounted for before draw)
                # Re-enable smart margins now that layout is finalized
                self._suppress_margin_update = False
                self._update_smart_margins(self.canvas)
                
                self.canvas.draw()
                
                self.canvas.repaint() # Force Qt repaint
        except Exception as e:
            msg = f"Plot Error: {e}"
            print(msg)
            print(traceback.format_exc())
            self.lbl_result.setText(msg)

    def get_input_dataframe(self, tab=None, file_path=None):
        try:
            # 1. Check loaded files in iolite session
            files = data.importedFiles()
            target_file = None
            
            if file_path:
                norm_target = os.path.normpath(file_path).lower()
                for f in files:
                    try:
                        f_path = os.path.normpath(f.filePath()).lower()
                        f_name = os.path.normpath(f.fileName()).lower()
                        target_name = os.path.normpath(os.path.basename(file_path)).lower()
                        if f_path == norm_target or f_name == target_name:
                            target_file = f
                            break
                    except Exception:
                        pass

            if not target_file and files:
                target_file = files[0]
            
            # 2. Get channels corresponding to target file
            if target_file is not None:
                ch_names = target_file.channelList()
                channels = []
                for name in ch_names:
                    ch = data.timeSeries(name)
                    if ch:
                        channels.append(ch)
            else:
                channels = data.timeSeriesList(data.Input)
                
            if not channels: return None, {}
            
            # --- PRE-COMPUTE TIME SLICING AND CYCLE TIME EARLY ---
            ref_ch = channels[0]
            time_data = ref_ch.time()
            
            mask = None
            if target_file is not None:
                try:
                    t_start = target_file.startTime().toSecsSinceEpoch()
                    t_end = target_file.endTime().toSecsSinceEpoch()
                    mask = (time_data >= t_start) & (time_data <= t_end)
                    time_data = time_data[mask]
                except Exception as me:
                    IoLog.warning(f"iolite Optimiser: Time slicing failed: {me}. Loading full series.")
                    mask = None

            # Calculate average time step for the sliced timeline
            detected_at_val = None
            if len(time_data) > 1:
                if len(time_data) > 2:
                    slope = np.polyfit(np.arange(len(time_data) - 1), time_data[1:], 1)[0]
                else:
                    slope = np.polyfit(np.arange(len(time_data)), time_data, 1)[0]
                detected_at_val = slope * 1000.0
            
            if tab == "spr":
                self.spr_detected_at_ms = detected_at_val
            else:
                self.detected_at_ms = detected_at_val
                self.opt_detected_at_ms = detected_at_val
            
            # --- Channel Filtering for Vitesse/icpTOF ---
            ref_ch_0 = channels[0]
            m_name = ref_ch_0.property("Machine Name")
            
            is_tof = False
            target_machines = ["Vitesse", "icpTOF"]
            
            if m_name:
                m_str = str(m_name)
                if any(tm.lower() in m_str.lower() for tm in target_machines):
                    is_tof = True
            
            if is_tof:
                ch_names = [ch.name for ch in channels]
                
                # Unified Dialog: Channel selection + Global dwell
                det_at = self.spr_detected_at_ms if tab == "spr" else self.detected_at_ms
                
                dlg = DataConfigDialog(ch_names, detected_at=det_at, is_tof=True, parent=self)
                if dlg.exec_():
                    selected_names = set(dlg.result_channels)
                    channels = [ch for ch in channels if ch.name in selected_names]
                    
                    # Store pre-configured dwells (resolve_dwell_times will use these)
                    self._preconfigured_dwells = dlg.result_dwells.copy()
                    
                    if not channels:
                        IoLog.warning("iolite Optimiser: No channels selected. Aborting load.")
                        return None, {}
                else:
                    IoLog.information("iolite Optimiser: Configuration cancelled. Aborting load.")
                    return None, {}
            else:
                self._preconfigured_dwells = {}  # Reset for non-TOF data
            
            # --- End Filter Logic ---
            
            data_dict = {'Time': time_data}
            
            local_metadata = {}
            has_el = False
            has_mass = False
            
            for ch in channels:
                if ch.name == 'TotalBeam': continue
                
                # Extract Metadata
                el = ch.property('Element')
                ma = ch.property('Mass')
                local_metadata[ch.name] = {'Element': el, 'Mass': ma}
                if el: has_el = True
                if ma: has_mass = True
                
                # Only load data if lengths match (crucial check)
                ch_data = ch.data()
                if mask is not None:
                    if len(ch_data) == len(mask):
                        data_dict[ch.name] = ch_data[mask]
                else:
                    if len(ch_data) == len(time_data):
                        data_dict[ch.name] = ch_data
            
            self.show_meta = {'Element': has_el, 'Mass': has_mass}
            
            # Also log file name being loaded for visibility
            if target_file is not None:
                IoLog.information(f"iolite Optimiser: Loaded {len(data_dict) - 1} channels from file: {target_file.fileName()}")
                
            return pd.DataFrame(data_dict), local_metadata
        except Exception as e:
            IoLog.error(f"iolite Optimiser: Input Error: {e}")
            return None, {}


    def on_auto_toggled(self, checked):
        if checked:
            self.run_auto_detect()
            
    def on_region_edited(self):
        # User manual edit -> Uncheck Auto
        self.chk_auto.blockSignals(True)
        try:
            self.chk_auto.setChecked(False)
        finally:
            self.chk_auto.blockSignals(False)
        
        # 1. Validate and Swap Inverted Regions (Start > End)
        bg_s_raw = self._get(self.spin_bg_start, 'value')
        bg_e_raw = self._get(self.spin_bg_end, 'value')
        bg_s, bg_e = sorted([bg_s_raw, bg_e_raw])
        
        sig_s_raw = self._get(self.spin_sig_start, 'value')
        sig_e_raw = self._get(self.spin_sig_end, 'value')
        sig_s, sig_e = sorted([sig_s_raw, sig_e_raw])
        
        # Sync internal state immediately (Fixes state de-sync on swap)
        t_shift = getattr(self, 't_start', 0)
        self.bg_times = (bg_s + t_shift, bg_e + t_shift)
        self.sig_times = (sig_s + t_shift, sig_e + t_shift)
        
        # Update spinboxes (Block signals to prevent recursion)
        for sb, val in [
            (self.spin_bg_start, bg_s), (self.spin_bg_end, bg_e),
            (self.spin_sig_start, sig_s), (self.spin_sig_end, sig_e)
        ]:
            if self._get(sb, 'value') != val:
                sb.blockSignals(True)
                sb.setValue(val)
                sb.blockSignals(False)

        # Read values (Relative Time)
        try:
            rel_bg_s = self._get(self.spin_bg_start, 'value')
            rel_bg_e = self._get(self.spin_bg_end, 'value')
            rel_sig_s = self._get(self.spin_sig_start, 'value')
            rel_sig_e = self._get(self.spin_sig_end, 'value')
            
            # Check for NaNs
            if any(x is None for x in [rel_bg_s, rel_bg_e, rel_sig_s, rel_sig_e]):
                return

            status = f"Manual: BG {rel_bg_s:.1f}-{rel_bg_e:.1f}s, SIG {rel_sig_s:.1f}-{rel_sig_e:.1f}s"
            self.lbl_result.setText(status)
            
            # Redundant update removed to prevents double-consumption of axis limits
            # self.update_plot(preserve_zoom=True)
            
            # Trigger recalc dynamically without saving settings to disk
            self.run_optimisation(refresh=False, preserve_zoom=True)
        except Exception as e:
            # Silently handle transient UI read errors
            pass
 
    def run_auto_detect(self):
        try:
            if self.opt_df is None:
                self.lbl_result.setText("No Input Channels found.")
                return
            
            # 2. Run Auto Detect (Returns Absolute Times)
            try:
                (bg_start, bg_end), (sig_start, sig_end) = Logic.auto_detect_regions(self.opt_df)
            except Exception as e:
                self.lbl_result.setText(f"Detection failed: {e}")
                return
            
            # Store in self (Absolute)
            self.bg_times = (bg_start, bg_end)
            self.sig_times = (sig_start, sig_end)
            
            # Convert to Relative for UI
            rel_bg_s = float(bg_start - self.t_start)
            rel_bg_e = float(bg_end - self.t_start)
            rel_sig_s = float(sig_start - self.t_start)
            rel_sig_e = float(sig_end - self.t_start)
            
            
            # Update SpinBoxes (Block signals to prevent loops)
            s_list = [self.spin_bg_start, self.spin_bg_end, self.spin_sig_start, self.spin_sig_end]
            for s in s_list: s.blockSignals(True)
            try:
                self.spin_bg_start.setValue(rel_bg_s)
                self.spin_bg_end.setValue(rel_bg_e)
                self.spin_sig_start.setValue(rel_sig_s)
                self.spin_sig_end.setValue(rel_sig_e)
            finally:
                for s in s_list: s.blockSignals(False)
            
            self.update_plot()
            # Trigger real-time optimization after regions detected
            self.run_optimisation(refresh=False)
            
        except Exception as e:
            self.lbl_result.setText(f"Detection error: {e}")
 
    def run_optimisation(self, refresh=False, silent=False, fixed_spot_um=None, preserve_zoom=False):
        if fixed_spot_um == 0:
            fixed_spot_um = None
            
        if refresh:
            # Force reset zoom on fresh data load
            preserve_zoom = False
            
        if refresh:
            IoLog.information("iolite Optimiser: Refreshing data...")
            self.lbl_result.setText("Refreshing Data...")
            self.refresh_data(tab="opt")
            IoLog.information(f"iolite Optimiser: Data refreshed. Data frame is {'None' if self.opt_df is None else 'Valid'}")

        # Cache theme for the entire run to ensure consistency
        # Uses cached value from showEvent/apply_theme instead of re-detecting
        if not hasattr(self, 'cached_is_dark'):
            self.cached_is_dark = (self.detect_theme() == 'dark')

        if self.opt_df is None:
             self.lbl_result.setText("No data loaded.")
             # ... (axes clearing logic)
             self.figure.clear()
             self.figure.add_subplot(111)
             self.canvas.draw()
             self.canvas.flush_events()
             self.canvas.repaint()
             self.table.setRowCount(0)
             self.lbl_opt_rr.setText("- Hz")
             self.lbl_opt_speed.setText("- µm s⁻¹")
             self.lbl_opt_at.setText("- ms")
             self.lbl_opt_pulses.setText("- Pulses")
             self.lbl_est_time_hms.setText("-")
             self.lbl_est_time_sec.setText("-")
             self.lbl_est_datapoints.setText("-")
             if not silent:
                 IoLog.information("iolite Optimiser: UI Cleared (No Data)")
             QApplication.processEvents()
             return

        # Double check for regions if data is present
        if not self.bg_times or not self.sig_times:
            self.lbl_result.setText("Please select Background and Signal regions.")
            return

        self.lbl_result.setText("Optimising...")
        try:
            # ... (spec fetching logic)
            mfr = self._get(self.cmb_mfr, 'currentText')
            model = self._get(self.cmb_model, 'currentText')
            spec = ICP_SPECS.get(mfr, {}).get(model, {})
            
            # ... (overrides logic)
            
            # ... 
            
            base_washout = self._get(self.spin_wash, 'value')
            if getattr(self, 'icp_tech', 'Quadrupole') == "TOF":
                margin_pct = getattr(self, 'persistent_settings', {}).get('tof_margin', 10.0)
                if hasattr(self, 'spin_tof_margin'):
                    margin_pct = self._get(self.spin_tof_margin, 'value')
                base_washout = base_washout * (1.0 + (margin_pct / 100.0))

            pulses_val = self._get(self.spin_pulses, 'value')
            rsd_val = self._get(self.spin_rsd, 'value') if hasattr(self, 'spin_rsd') else 5.0

            c = {
                'icp_mfr': mfr,
                'icp_model': model,
                'icp_technology': self.icp_tech,
                'system_type': self.icp_tech, # Required by calculate_sigma_for_spot
                'precision_ms': self.precision,
                'min_dwell_ms': self.min_dwell,
                'max_rr_hz': self.max_rr,
                'allowed_rr': self.allowed_rr,
                'allowed_dwells': self.allowed_dwells,
                'max_speed_um_s': self.max_speed,
                'rr_prec_hz': self._get(self.spin_cust_rr_prec, 'value') if hasattr(self, 'spin_cust_rr_prec') else 1.0,
                'washout_ms': base_washout,
                'spot_size_um': fixed_spot_um if fixed_spot_um else self._get(self.spin_spot, 'value'),
                'pulses_per_pixel': pulses_val,
                'target_rsd': rsd_val,
                'dosage': self._get(self.spin_dosage, 'value'),
                'mode': self._get(self.cmb_mode, 'currentText'),
                'img_height': self._get(self.spin_height, 'value') if hasattr(self, 'spin_height') else 1000.0,
                'img_width': self._get(self.spin_width, 'value') if hasattr(self, 'spin_width') else 1000.0,
                'line_length': self._get(self.spin_line_len, 'value') if hasattr(self, 'spin_line_len') else 1000.0,
                'avoid_gaps': self._get(self.chk_avoid_gaps, 'isChecked') if hasattr(self, 'chk_avoid_gaps') else False,
                'prefer_exact_pulses': self._get(self.chk_prefer_exact_pulses, 'isChecked') if hasattr(self, 'chk_prefer_exact_pulses') else False,
                'lower_sigma_limit': self._get(self.spin_sigma, 'value'),
                'min_duty_cycle': self._get(self.spin_duty, 'value') / 100.0 if self._get(self.spin_duty, 'value') else 0.0,
                'snr_threshold': self._get(self.spin_snr, 'value'),
                'overhead_ms': getattr(self, 'detected_overhead_ms', 0.0),
                'ref_spot_size_um': self._get(self.spin_spot, 'value'),
                'scale_signal': self._get(self.chk_scale, 'isChecked'),
                'initial_rr': self._get(self.spin_init_rr, 'value'),
                'sync_strategy': self._get(self.cmb_sync_strategy, 'currentText')
            }
            
            bg_s, bg_e = self.bg_times if self.bg_times else (0, 1)
            sig_s, sig_e = self.sig_times if self.sig_times else (2, 3)
            
            df_blk = self.opt_df[(self.opt_df['Time'] >= bg_s) & (self.opt_df['Time'] <= bg_e)]
            df_sig = self.opt_df[(self.opt_df['Time'] >= sig_s) & (self.opt_df['Time'] <= sig_e)]
            
            isotope_data = []
            
            missing_dwell_channels = []
        
            # Ensure dwell source exists
            dwells_source = getattr(self, 'opt_dwells', getattr(self, 'channel_dwells', {}))
            
            # 1. Pre-Check for Dwell Metadata
            for col in self.opt_df.columns:
                if col == 'Time': continue
                if col not in dwells_source:
                     missing_dwell_channels.append(col)
                     
            if missing_dwell_channels:
                 # ... (abort logic)
                 msg = "Optimisation Aborted: Missing Dwell Times"
                 IoLog.warning(msg)
                 self.lbl_result.setText(msg)
                 return
    
            # 2. Process Data (All dwells guaranteed to exist)
            for col in self.opt_df.columns:
                if col == 'Time': continue

            
            # Fetch Resolved Dwell
                # Fetch Resolved Dwell
                ref_dt_ms = dwells_source[col]
                ref_dt_s = ref_dt_ms / 1000.0
                
                m_sig = df_sig[col].mean()
                m_blk = df_blk[col].mean()
                s_blk = df_blk[col].std()
                
                is_cps = self.channel_is_cps.get(col, False) if hasattr(self, 'channel_is_cps') else False
                
                if is_cps:
                    sig_cps = m_sig
                    blk_cps = m_blk
                    std_blk_counts = s_blk * ref_dt_s
                    
                    # For SNR scaling (Counts Math)
                    math_sig = m_sig * ref_dt_s
                    math_blk = m_blk * ref_dt_s
                    math_std_blk = s_blk * ref_dt_s
                else:
                    # Input is Counts -> Convert to CPS
                    sig_cps = m_sig / ref_dt_s
                    blk_cps = m_blk / ref_dt_s
                    std_blk_counts = s_blk
                    
                    # For SNR scaling (Counts Math)
                    math_sig = m_sig
                    math_blk = m_blk
                    math_std_blk = s_blk

                # Calculate Robust SNR (summing variances in quadrature in Counts domain)
                # Matches logic.py reference: sqrt(Poisson_Variance + Flicker_Variance)
                # Poisson Var = Mean Counts (math_blk)
                # Flicker Var = StdDev Counts^2 (math_std_blk^2)
                # Note: This follows user's logic that these sum in quadrature.
                noise = (max(0, math_blk) + math_std_blk**2)**0.5
                initial_snr = (math_sig - math_blk) / noise if noise > 0 else 0

                # Fetch session-specific config for this isotope
                iso_conf = self.isotope_configs.setdefault(col, {"status": "Auto", "custom_time_s": 0.01})
                
                isotope_data.append({
                    'name': col,
                    'sig_cps': sig_cps,
                    'blk_cps': blk_cps,
                    'stdev_blank_counts': std_blk_counts if s_blk is not None else 0,
                    'baseline_dt_s': ref_dt_s,
                    'dwell': ref_dt_s, 
                    'status': iso_conf['status'],
                    'custom_time_s': iso_conf['custom_time_s'],
                    'initial_snr': initial_snr
                })

            # --- HARMONIC SCALING PRE-OPTIMIZATION ---
            is_simultaneous = TYPE_TO_TECH_MAP.get(self.icp_tech) in ["Multi-Collector", "TOF"]
            min_dwell_total = 0.0
            if isotope_data:
                if is_simultaneous:
                    # MC/TOF: Channels are simultaneous. The budget needed is the maximum 
                    # requirement of any active channel, not the sum.
                    reqs = []
                    for iso in isotope_data:
                        if iso.get('status') == "Exclude": continue
                        if iso.get('status') == "Custom":
                            reqs.append(iso.get('custom_time_s', 0) * 1000.0)
                        else:
                            reqs.append(self.min_dwell)
                    min_dwell_total = max(reqs) if reqs else self.min_dwell
                else:
                    # Sequential systems: Sum of all active dwells
                    for iso in isotope_data:
                        stat = iso.get('status', 'Auto')
                        if stat == "Exclude": continue
                        if stat == "Custom":
                            min_dwell_total += iso.get('custom_time_s', self.min_dwell / 1000.0) * 1000.0
                        else:
                            min_dwell_total += self.min_dwell
            c['min_dwell_needed_ms'] = min_dwell_total

            # Count Even channels and calculate fixed costs for perfect split synchronization
            num_even = 0
            fixed_costs_s = 0.0
            min_dwell_s = self.min_dwell / 1000.0
            
            for iso in isotope_data:
                stat = iso.get('status', 'Auto')
                if stat == "Exclude": continue
                if stat == "Even":
                    num_even += 1
                elif stat == "Custom":
                    fixed_costs_s += iso.get('custom_time_s', min_dwell_s)
                elif stat == "Set to Min":
                    fixed_costs_s += min_dwell_s
                else: # Auto or non-optimisable
                    # For synchronization purposes, Auto pool items are treated as 
                    # minimums until the budget is optimized.
                    fixed_costs_s += min_dwell_s

            c['num_even'] = num_even
            c['fixed_costs_s'] = fixed_costs_s

            active_dwells_ms = []
            for iso in isotope_data:
                stat = iso.get('status', 'Auto')
                if stat == "Exclude":
                    continue
                elif stat == "Custom":
                    active_dwells_ms.append(iso.get('custom_time_s', min_dwell_s) * 1000.0)
                elif stat == "Set to Min":
                    active_dwells_ms.append(self.min_dwell)
                else:
                    active_dwells_ms.append(self.min_dwell)
            c['active_dwells_ms'] = active_dwells_ms

            c['log_prefix'] = "Initial Sync"
            sync = Logic.calculate_constrained_at(c)
            
            cols = [c for c in self.opt_df.columns if c != 'Time']
            if not cols:
                 self.lbl_result.setText("No valid data channels.")
                 # Clear UI if only Time remains (e.g. TotalBeam filtered out)
                 self.figure.add_subplot(111)
                 self.canvas.draw()
                 self.canvas.flush_events()
                 self.canvas.repaint()
                 self.table.setRowCount(0)
            # self.lbl_target_at.setText(f"Target: {sync.get('Target Acquisition Time (ms)', 0):.1f} ms") # REMOVED
            
            # Labels will be updated after spot size optimization below.

            # 3. Step 1: Always Calculate Minimum Required Spot Size (for reference)
            # Need to capture the sync info from the last step of the optimization if we want it perfect
            optimum_spot_raw, df_res_opt = Logic.calculate_minimum_required_spot_size(c, isotope_data)
            # To be safe, we re-run once for the final spot to get the matching sync result
            _, _, sync = Logic.calculate_sigma_for_spot(int(round(optimum_spot_raw)), c, isotope_data)
            optimum_spot = int(round(optimum_spot_raw))
            self.optimum_spotsize = optimum_spot
            self.override_spotsize = fixed_spot_um
            
            best_spot = optimum_spot
            df_res = df_res_opt
            
            # Identify limiting channel for the optimum spot
            limiting_iso = "None"
            if df_res_opt:
                # Filter for candidates that were part of the optimization (Auto or Even)
                opt_rows = [r for r in df_res_opt if r.get('Status') in ["Auto", "Even"]]
                if not opt_rows: opt_rows = df_res_opt # Fallback
                if opt_rows:
                    min_row = min(opt_rows, key=lambda x: x['Sigma Sep'])
                    limiting_iso = min_row['Isotope']
            
            final_notes = []
            if fixed_spot_um is not None:
                # Manual Override
                best_spot = fixed_spot_um
                # Re-calculate table for this specific spot
                _, df_res, sync = Logic.calculate_sigma_for_spot(best_spot, c, isotope_data)
                final_notes.append(f"Optimum Spot Size <b>{optimum_spot}</b> µm - Overridden to <b>{best_spot}</b> µm")
            else:
                # Update Spinbox (Block signals to prevent recursion)
                self.spin_opt_spot.blockSignals(True)
                self.spin_opt_spot.setValue(best_spot)
                self.spin_opt_spot.blockSignals(False)
                final_notes.append(f"Optimum Spot Size <b>{best_spot}</b> µm - Based on Minimum SNR of <b>{limiting_iso}</b>")
            # 4. Step 2: Laser Sync for that spot size
            c['spot_size_um'] = best_spot
            
            # CRITICAL: For MC systems, update the 'min_dwell_needed_ms' to reflect 
            # the ACTUAL selected dwell time (e.g. 131ms) from the optimization results.
            # Otherwise, calculate_laser_sync uses the initial minimum (66ms) and reports the wrong budget.
            if df_res is not None:
                is_mc = (self.icp_tech == "Multi-Collector")
                if is_mc:
                     max_opt_dwell = max([r.get('Final Dwell (ms)', 0) for r in df_res], default=0)
                     if max_opt_dwell > 0:
                         c['min_dwell_needed_ms'] = max_opt_dwell
                
                active_dwells = [r.get('Final Dwell (ms)', 0.0) for r in df_res if r.get('Status') != "Exclude"]
                if active_dwells:
                    c['active_dwells_ms'] = active_dwells
            
            # Reset 'pulses_per_pixel' to user input to ensure Error Calculation compares against Target
            c['pulses_per_pixel'] = pulses_val
            c['dosage'] = self._get(self.spin_dosage, 'value')
            c['log_prefix'] = "Optimised Sync"
            sync = Logic.calculate_constrained_at(c)
            
            # --- AUTO-SYNC SWEEP LOGIC ---
            # Now natively handled by Logic.calculate_constrained_at(c)
            c['override_at_s'] = sync['at_actual_s']
            c['override_rr'] = sync['rr_actual']
            
            # Re-run final pass of calculate_sigma_for_spot to update the table dwells, speed, overlap, etc.!
            _, df_res, sync = Logic.calculate_sigma_for_spot(best_spot, c, isotope_data)
            
            # STORE SYNC FOR PULSE TRAIN SIMULATOR
            self.last_sync = sync
            
            # Aggregate notes from Logic
            logic_notes = sync.get("notes", [])
            final_notes.extend(logic_notes)
            
            # Aggregate constraints from isotopes
            if df_res is not None:
                min_snr_isos = [row['Isotope'] for row in df_res if row.get('Constraint') == "Min SNR"]
                min_icp_isos = [row['Isotope'] for row in df_res if row.get('Constraint') == "Min ICP"]
                zero_bg_isos = [row['Isotope'] for row in df_res if row.get('IsZeroBG') is True]
                
                target_sigma = c.get('lower_sigma_limit', 0)
                snr_thresh = c.get('snr_threshold', 0)
                min_dwell_ms = c.get('min_dwell_ms', 0)
                
                if min_snr_isos:
                    is_mc_or_tof = (self.icp_tech in ["Multi-Collector", "TOF"])
                    if is_mc_or_tof:
                        final_notes.append(("The following channels do not meet the minimum SNR target:", ", ".join(min_snr_isos), "orange"))
                    else:
                        final_notes.append(("The following channels could not reach the minimum SNR target and have been set to the minimum dwell time:", ", ".join(min_snr_isos), "orange"))
                if min_icp_isos:
                    final_notes.append((f"The following channels cannot be set lower than the hardware minimum&nbsp;({min_dwell_ms}&nbsp;ms):", ", ".join(min_icp_isos), "blue"))
                
                if zero_bg_isos:
                    final_notes.append(("The following channels had zero counts in the background. The detection limit (L<sub>c</sub>) has been calculated using the Square Root Transform rule for variance stabilization:", ", ".join(zero_bg_isos), "gray"))
            

            
            # Construct final status summary
            summary = "<b>Optimisation Complete</b>"
            if final_notes:
                for note in final_notes:
                    if isinstance(note, tuple):
                        if len(note) == 3 and note[1] is None:
                             # Legacy support or custom tuple without detail
                             msg, _, colour_key = note
                             c_hex = "#FFFFFF" if colour_key == "white" else "#FF6D00"
                             summary += f"<div style='margin-left: 15px; text-indent: -15px;'><span style='color: {c_hex}'>•</span> {msg}</div>"
                        else:
                            hdr, det, colour_key = note
                            if colour_key == "blue": marker_colour = "#3399FF"
                            elif colour_key == "red": marker_colour = "#FF0000"
                            else: marker_colour = "#FF6D00" # Default/Orange
                            summary += f"<div style='margin-left: 15px; text-indent: -15px;'><span style='color: {marker_colour}'>•</span> {hdr}</div>"
                            summary += f"<div style='margin-left: 30px; text-indent: -15px;'><span style='color: {marker_colour}'>•</span> <b>{det}</b></div>"
                    else:
                        summary += f"<div style='margin-left: 15px; text-indent: -15px;'>• {note}</div>"
            
            # Save final status base
            self.final_status_base = summary
            
            # --- SYNCHRONIZATION ADVISOR DISPLAY ---
            advisor_html = "<br><b>Synchronization Advisor:</b><br>"
            
            # Check if perfect sync (within 1e-9 tolerance)
            F = sync['rr_actual']
            AT_s = sync['at_actual_s']
            pulses_per_cycle = F * AT_s
            is_perfect = abs(pulses_per_cycle - round(pulses_per_cycle)) < 1e-9
            
            strategy = c.get('sync_strategy', 'Adaptive Integer Sync (Auto)')
            is_oversampling_mode = "Oversampling" in strategy or "Combined" in strategy
            
            if is_oversampling_mode:
                target_rsd = c.get('target_rsd', 5.0)
                washout_ms = c.get('washout_ms', 30.0)
                
                # Retrieve active dwells (resolved or initial)
                rsd_details = []
                if df_res is not None:
                    for r in df_res:
                        if r.get('Status') != "Exclude":
                            iso_name = r.get('Isotope', 'Unknown')
                            dt_val = r.get('Final Dwell (ms)', 0.0)
                            if dt_val > 0:
                                rsd_c = Logic.calculate_lockwood_rsd(F, washout_ms, dt_val)
                                rsd_details.append((iso_name, dt_val, rsd_c))
                
                if not rsd_details:
                    at_ms = sync['at_actual_s'] * 1000.0
                    rsd_c = Logic.calculate_lockwood_rsd(F, washout_ms, at_ms)
                    rsd_details.append(('Cycle Time', at_ms, rsd_c))
                
                worst_iso, worst_dwell, worst_rsd = max(rsd_details, key=lambda x: x[2])
                
                if is_perfect:
                    advisor_html += f"<span style='color: #4CAF50;'>🟢 Synchronised ({round(pulses_per_cycle)} pulses per acq)</span><br>"
                else:
                    advisor_html += f"<span style='color: #8BC34A;'>🟢 Oversampling Mode (No integer snapping required)</span><br>"
                    
                    overlap_factor = F * (washout_ms / 1000.0)
                    if overlap_factor < 2.0:
                        advisor_html += f"<span style='color: #F44336;'>🔴 Unstable Steady State (Rep Rate {F:.1f} Hz is too low for oversampling. Must be at least {2.0 / (washout_ms/1000.0):.1f} Hz to ensure pulse overlap with {washout_ms:.1f} ms washout)</span><br>"
                    elif worst_rsd > target_rsd:
                        advisor_html += f"<span style='color: #FFC107;'>🟡 Unstable Steady State (Worst RSD is {worst_rsd:.1f}% for {worst_iso} at {worst_dwell:.1f} ms dwell - target is {target_rsd:.1f}%)</span><br>"
                        suggested_dt = Logic.get_lockwood_min_dwell(F, washout_ms, target_rsd)
                        advisor_html += f"• Increase {worst_iso} Dwell/Acquisition Time to: {suggested_dt:.1f} ms to lower RSD below {target_rsd:.1f}%<br>"
                    else:
                        advisor_html += f"<span style='color: #4CAF50;'>🟢 Stable Steady State (All channels below {target_rsd:.1f}% target RSD. Worst is {worst_rsd:.1f}% for {worst_iso})</span><br>"
            else:
                # Standard Synchronised / Snapped modes
                if is_perfect:
                    advisor_html += f"<span style='color: #4CAF50;'>🟢 Synchronised ({round(pulses_per_cycle)} pulses per acq)</span><br>"
                    target_pulses = pulses_val
                    final_N = round(pulses_per_cycle)
                    if "Auto" in strategy and abs(final_N - target_pulses) >= 1e-9:
                        advisor_html += f"• Adjusted target from {target_pulses:.1f} to {final_N:.1f} pulses per acq<br>"
                else:
                    diff_pulse = abs(pulses_per_cycle - round(pulses_per_cycle))
                    beat_s = AT_s / diff_pulse if diff_pulse > 1e-9 else 0.0
                    beat_text = f" ({beat_s:.1f}s beat pattern)" if beat_s > 0 else ""
                    
                    advisor_html += f"<span style='color: #FFC107;'>🟡 Unsynchronised ({pulses_per_cycle:.4f} pulses per acq{beat_text})</span><br>"
                    

                    
            self.final_status_base += advisor_html
            
            # Copy calculated pulses per cycle to dosage box if Sync Dosage is checked
            if hasattr(self, 'chk_sync_dosage') and self.chk_sync_dosage.isChecked():
                self._block_signals(self.spin_dosage, True)
                self.spin_dosage.setValue(pulses_per_cycle)
                self._block_signals(self.spin_dosage, False)
                
            # SNAP TO HARDWARE PRECISION at the very end
            prec_ms = self.precision
            # Local display variables derived from the 'sync' master result
            val_at_ms = round(sync['at_actual_s'] * 1000.0 / prec_ms) * prec_ms
            val_budget_ms = round((sync['at_actual_s'] * 1000 - getattr(self, 'detected_overhead_ms', 0.0)) / prec_ms) * prec_ms
            val_overhead_ms = round(getattr(self, 'detected_overhead_ms', 0.0) / prec_ms) * prec_ms
            
            
            # 5. Display Summary
            unit, scaler, disp_dps = self._get_dwell_unit_info()
            laser_dps = self._get_laser_precision_dps()

            # --- DYNAMIC TABLE HEADERS ---
            headers = ["Channel"]
            if self.show_meta.get('Element'): headers.append("Element")
            if self.show_meta.get('Mass'): headers.append("Mass")
            headers.append("Mode")
            
            dwell_hdr = f"Optimised Dwell Times ({unit})"
            headers.extend([dwell_hdr, "Initial SNR", "Resultant SNR"])
            
            self.col_map = {name: i for i, name in enumerate(headers)}
            self.table.setColumnCount(len(headers))
            self.table.setHorizontalHeaderLabels(headers)
            
            # Reset Resizing behavior for dynamic layout
            self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.table.horizontalHeader().setSectionResizeMode(self.col_map['Channel'], QHeaderView.ResizeToContents)
            if 'Element' in self.col_map:
                self.table.horizontalHeader().setSectionResizeMode(self.col_map['Element'], QHeaderView.ResizeToContents)
            if 'Mass' in self.col_map:
                self.table.horizontalHeader().setSectionResizeMode(self.col_map['Mass'], QHeaderView.ResizeToContents)
            self.table.horizontalHeader().setSectionResizeMode(self.col_map['Mode'], QHeaderView.ResizeToContents)

            self.lbl_opt_rr.setText(f"<b>{sync['rr_actual']:.{laser_dps}f}</b> Hz")
            self.lbl_opt_speed.setText(f"<b>{sync['speed']:.1f}</b> µm s⁻¹")
            
            unit_val_at = val_at_ms / scaler
            unit_val_budget = val_budget_ms / scaler
            unit_val_overhead = val_overhead_ms / scaler

            self.lbl_opt_at.setText(f"<b>{unit_val_at:.{disp_dps}f}</b> {unit}")
            self.lbl_opt_budget.setText(f"<b>{unit_val_budget:.{disp_dps}f}</b> {unit}")
            self.lbl_opt_overhead.setText(f"<b>{unit_val_overhead:.{disp_dps}f}</b> {unit}")

            # Calculate Duty Cycle (Sum of Optimised Dwells / AT)
            # Use sum of actual optimized dwells (from df_res) for sequential, max for simultaneous
            sum_opt_dwells = 0.0
            if df_res is not None:
                dwells = [r['Final Dwell (ms)'] for r in df_res]
                is_mc_or_tof = (self.icp_tech in ["Multi-Collector", "TOF"])
                
                if is_mc_or_tof:
                    # Simultaneous: Duty = Max Dwell / AT
                    sum_opt_dwells = max(dwells) if dwells else 0.0
                else:
                    # Sequential: Duty = Sum Dwells / AT
                    sum_opt_dwells = sum(dwells)
                
            raw_at = val_at_ms
            
            duty_cycle = 0.0
            if raw_at > 0:
                duty_cycle = (sum_opt_dwells / raw_at) * 100.0
            
            self.lbl_opt_duty.setText(f"<b>{duty_cycle:.1f}</b> %")



            if not self._get(self.chk_avoid_gaps, 'isChecked'):
                 # "Avoid Gaps" OFF -> Floor Strategy.
                 self.lbl_opt_pulses.setText(f"<b>{sync['n_dose']:.1f}</b> Shots")
            else:
                 # "Avoid Gaps" ON -> Ceil Strategy (Overlap).
                 self.lbl_opt_pulses.setText(f"<b>{sync['n_dose']:.1f}</b> Shots")

            self.lbl_opt_overlap.setText(f"<b>{sync['overlap_um']:.1f}</b> µm")

            # Estimated Scan Time
            # Mode Logic
            mode = self._get(self.cmb_mode, 'currentText')
            est_sec = 0.0
            
            try:
                if mode == "Spot":
                    # Spot: Time = Shots / RepRate
                    shots = float(self._get(self.spin_shots, 'value'))
                    rr = sync.get('rr_actual', 1.0)
                    if rr > 0:
                        est_sec = shots / rr
                        
                elif mode == "Line":
                    # Line: Time = Length / Speed
                    # Overhead is already included/negligible per user request
                    img_w = float(self._get(self.spin_width, 'value')) # Acts as Length
                    speed = float(sync.get("speed", 0.0))
                    
                    if speed > 1e-6:
                         est_sec = img_w / speed
                         
                else:
                    # Imaging (Default)
                    # Lines * (Length / Speed)
                    # Overhead removed per user request
                    img_h = float(self._get(self.spin_height, 'value'))
                    img_w = float(self._get(self.spin_width, 'value'))
                    spot_um = float(c.get('spot_size_um', 0.0))
                    speed_ums = float(sync.get("speed", 0.0))
                    
                    if spot_um > 1e-6 and speed_ums > 1e-6:
                        lines = math.ceil(img_h / spot_um)
                        line_time_s = img_w / speed_ums
                        est_sec = lines * line_time_s

                # Format HH:MM:SS
                # Standardize to nearest integer second for consistency
                total_sec_int = int(round(est_sec))
                
                hrs = total_sec_int // 3600
                mins = (total_sec_int % 3600) // 60
                secs = total_sec_int % 60
                
                self.lbl_est_time_hms.setText(f"<b>{hrs:02d}:{mins:02d}:{secs:02d}</b>")
                self.lbl_est_time_sec.setText(f"<b>{total_sec_int}</b> s")

                # Est # Data Points (for spot and line modes) vs Pixels Per Sec (for imaging mode)
                at_actual_s = sync.get('at_actual_s', 0.1)
                if mode in ["Spot", "Line"]:
                    self.lbl_est_datapoints_title.setText("Est. Data Points:")
                    if at_actual_s > 0:
                        datapoints = est_sec / at_actual_s
                        self.lbl_est_datapoints.setText(f"<b>{int(round(datapoints))}</b>")
                    else:
                        self.lbl_est_datapoints.setText("-")
                else:
                    # Imaging mode
                    self.lbl_est_datapoints_title.setText("Pixels Per Sec:")
                    if at_actual_s > 0:
                        self.lbl_est_datapoints.setText(f"<b>{1.0 / at_actual_s:.1f}</b> px s⁻¹")
                    else:
                        self.lbl_est_datapoints.setText("-")
                
            except Exception as e:
                self.lbl_est_time_hms.setText("Error")
                self.lbl_est_time_sec.setText(f"({str(e)})")
                self.lbl_est_datapoints.setText("-")
            
            self.lbl_result.setText(self.final_status_base)
            
            # 6. Display Table
            self.optimised_results = df_res # Store queryable results for Pulse Train Simulator
            
            # Map initial dwells for easy lookup
            init_dwell_map = {d['name']: d['dwell'] for d in isotope_data}

            self.table.blockSignals(True) # Prevent infinite loops
            
            # Smart Update: Only rebuild structure if row count changes
            # This prevents flickering of ComboBoxes
            rows_needed = len(df_res)
            # Handle rowCount as property (int) or method depending on binding
            rc = self.table.rowCount
            rows_current = rc() if callable(rc) else rc
            
            if rows_needed != rows_current:
                self.table.setRowCount(0)
                self.table.setRowCount(rows_needed)
                rebuild_widgets = True
            else:
                rebuild_widgets = False
            
            for i, row in enumerate(df_res):
                iso_name = str(row['Isotope'])
                val_final = row['Final Dwell (ms)'] / scaler
                
                # Channel Name
                item_ch = self.table.item(i, self.col_map['Channel'])
                if not item_ch:
                    item_ch = QTableWidgetItem(iso_name)
                    item_ch.setTextAlignment(Qt.AlignCenter)
                    self.table.setItem(i, self.col_map['Channel'], item_ch)
                else:
                    item_ch.setText(iso_name)
                
                # Metadata (if visible)
                # Ensure metadata dict exists
                meta_source = getattr(self, 'opt_metadata', {})
                
                if self.show_meta.get('Element'):
                    el = meta_source.get(iso_name, {}).get('Element', '')
                    item_el = self.table.item(i, self.col_map['Element'])
                    if not item_el:
                        item_el = QTableWidgetItem(str(el) if el else "")
                        item_el.setTextAlignment(Qt.AlignCenter)
                        self.table.setItem(i, self.col_map['Element'], item_el)
                    else:
                        item_el.setText(str(el) if el else "")

                if self.show_meta.get('Mass'):
                    ma = meta_source.get(iso_name, {}).get('Mass', '')
                    item_ma = self.table.item(i, self.col_map['Mass'])
                    if not item_ma:
                        item_ma = QTableWidgetItem(str(ma) if ma else "")
                        item_ma.setTextAlignment(Qt.AlignCenter)
                        self.table.setItem(i, self.col_map['Mass'], item_ma)
                    else:
                        item_ma.setText(str(ma) if ma else "")
                
                # Mode (Dropdown) - Reuse existing if possible
                current_status = self.isotope_configs.get(iso_name, {}).get("status", "Auto")
                
                # Determine allowed mode items based on hardware technology
                allowed_items = ["Auto", "Even", "Set to Min", "Exclude", "Custom"]
                tech = getattr(self, 'icp_tech', 'Quadrupole')
                if tech == "Quadrupole":
                    allowed_items.remove("Exclude")
                    if current_status == "Exclude":
                        current_status = "Auto"
                        self.isotope_configs.setdefault(iso_name, {})["status"] = "Auto"
                elif tech in ["TOF", "Multi-Collector"]:
                    allowed_items.remove("Set to Min")
                    if current_status == "Set to Min":
                        current_status = "Auto"
                        self.isotope_configs.setdefault(iso_name, {})["status"] = "Auto"

                combo = self.table.cellWidget(i, self.col_map['Mode'])
                if not combo or rebuild_widgets:
                    combo = QComboBox()
                    combo.addItems(allowed_items)
                    combo.setCurrentText(current_status)
                    combo.setProperty("iso_name", str(iso_name))
                    combo.activated.connect(lambda idx, n=iso_name, c=combo: self._on_mode_changed(n, c.itemText(idx)))
                    self.table.setCellWidget(i, self.col_map['Mode'], combo)
                else:
                    # Check if the combo items need to be refreshed
                    existing_items = [combo.itemText(j) for j in range(combo.count)]
                    if existing_items != allowed_items:
                        combo.blockSignals(True)
                        combo.clear()
                        combo.addItems(allowed_items)
                        combo.setCurrentText(current_status)
                        combo.blockSignals(False)
                    elif self._get(combo, 'currentText') != current_status:
                        combo.blockSignals(True)
                        combo.setCurrentText(current_status)
                        combo.blockSignals(False)
                
                # Optimised Dwell (Editable if Custom)
                item_final = QTableWidgetItem(f"{val_final:.{disp_dps}f}")
                item_final.setTextAlignment(Qt.AlignCenter)
                font = item_final.font()
                font.setBold(True)
                item_final.setFont(font)
                c_idx = self.col_map[dwell_hdr]
                
                constraint = row.get('Constraint')

                if current_status == "Custom":
                    item_final.setFlags(item_final.flags() | Qt.ItemIsEditable)
                    is_dark = self.system_is_dark
                    highlight_color = QColor("#554400") if is_dark else QColor("#FFF8DC")
                    item_final.setBackground(highlight_color)
                elif constraint == "Min SNR":
                    item_final.setFlags(item_final.flags() & ~Qt.ItemIsEditable)
                    # More distinct orange warning for SNR failure
                    is_dark = self.system_is_dark
                    highlight_color = QColor("#853D00") if is_dark else QColor("#FFCCBC")
                    item_final.setBackground(highlight_color)
                elif constraint == "Min ICP":
                    item_final.setFlags(item_final.flags() & ~Qt.ItemIsEditable)
                    # Blue info for hardware minimum
                    is_dark = self.system_is_dark
                    highlight_color = QColor("#1E3A5F") if is_dark else QColor("#E3F2FD")
                    item_final.setBackground(highlight_color)
                else:
                    item_final.setFlags(item_final.flags() & ~Qt.ItemIsEditable)
                
                self.table.setCellWidget(i, c_idx, None)
                self.table.setItem(i, c_idx, item_final)
                
                # Initial SNR
                item_snr = QTableWidgetItem(f"{row['Initial SNR']:.1f}")
                item_snr.setTextAlignment(Qt.AlignCenter)
                item_snr.setFlags(item_snr.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(i, self.col_map['Initial SNR'], item_snr)
                
                # Resultant SNR
                val_sigma = row['Sigma Sep']
                if constraint == "Min SNR":
                    txt_sigma = "Below Minimum SNR"
                elif val_sigma < 0:
                    txt_sigma = "Undetectable"
                else:
                    txt_sigma = f"{val_sigma:.1f}"
                
                item_sigma = QTableWidgetItem(txt_sigma)
                item_sigma.setTextAlignment(Qt.AlignCenter)
                item_sigma.setFlags(item_sigma.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(i, self.col_map['Resultant SNR'], item_sigma)
            
            self.table.blockSignals(False)
            
            # Connect itemChanged for Custom dwell edits if not already connected
            try: self.table.itemChanged.disconnect(self._on_dwell_changed)
            except: pass
            self.table.itemChanged.connect(self._on_dwell_changed)
            
            # TRIGGER PLOT UPDATE (Reflecting changes in scaling/labels)
            self.update_plot(preserve_zoom=preserve_zoom)
                
        except Exception as e:
            msg = f"Optimisation Error: {str(e)}"
            IoLog.error(msg)
            IoLog.error(traceback.format_exc())
            self.lbl_result.setText(msg)

    def _on_mode_changed(self, iso_name, text):
        IoLog.information(f"iolite Optimiser: Mode Changed: {iso_name} -> {text}")
        
        if iso_name in self.isotope_configs:
            # If switching to Custom, try to seed with the current result from the table
            if text == "Custom":
                try:
                    # Determine dynamic column for optimized dwell
                    unit, scaler, _ = self._get_dwell_unit_info()
                    dwell_hdr = f"Optimised Dwell Times ({unit})"
                    d_col = self.col_map.get(dwell_hdr, -1)

                    if d_col != -1:
                        # Find row for this isotope (Channel column is always first)
                        r_count = self._get(self.table, 'rowCount')
                        for r in range(0 if r_count is None else r_count):
                            if self._get(self.table.item(r, 0), 'text') == iso_name:
                                current_val_text = self._get(self.table.item(r, d_col), 'text')
                                val = float(current_val_text)
                                
                                # Convert displayed value back to seconds (internal storage unit)
                                self.isotope_configs[iso_name]["custom_time_s"] = (val * scaler) / 1000.0
                                break
                except Exception as e:
                    IoLog.warning(f"Failed to seed custom value for {iso_name}: {e}")

            self.isotope_configs[iso_name]["status"] = text
            
            # Synchronize Auto/Even pool
            if text in ["Auto", "Even"]:
                target_mode = text
                for name, conf in self.isotope_configs.items():
                    if conf.get('status') in ["Auto", "Even"]:
                        conf['status'] = target_mode

            # Trigger real-time optimization
            self.run_optimisation(refresh=False, preserve_zoom=True)

    def _on_dwell_changed(self, item):
        column = item.column()
        row = item.row()
        
        # Determine dynamic column for optimized dwell
        unit, scaler, _ = self._get_dwell_unit_info()
        dwell_hdr = f"Optimised Dwell Times ({unit})"
        d_col = self.col_map.get(dwell_hdr, -1)

        # We only care about edits in the "Optimised Dwell" column
        if d_col == -1 or column != d_col: return
        
        # Get the isotope name from the first column
        iso_item = self.table.item(row, 0)
        if not iso_item: return
        iso_name = self._get(iso_item, 'text')
        
        # Ensure it's in "Custom" mode
        status = self.isotope_configs.get(iso_name, {}).get("status")
        if status != "Custom": return
        
        try:
            val_text = self._get(item, 'text')
            new_val = float(val_text)
            
            unit, scaler, _ = self._get_dwell_unit_info()
            
            IoLog.information(f"iolite Optimiser: Custom Dwell Adjusted (Text): {iso_name} -> {new_val} {unit}")
            
            # Convert UI value to seconds (internal storage unit)
            self.isotope_configs[iso_name]["custom_time_s"] = (new_val * scaler) / 1000.0
                
            # Trigger real-time optimization
            self.run_optimisation(refresh=False, preserve_zoom=True)
        except ValueError:
            # Re-run to reset the cell to its previously calculated valid value
            self.run_optimisation(refresh=False, preserve_zoom=True)
        except Exception as e:
            IoLog.error(f"Dwell Edit Error: {e}")

    def closeEvent(self, event):
        """Override close event to save settings."""
        self.save_persistent_settings()
        event.accept()

    def show_prescan_suggestions(self):
        try:
            num_analytes = self.spin_num_analytes.value
            washout_ms = self.spin_wash.value
            spot_size = self.spin_spot.value

            # Dwell / Cycle time resolution
            min_dwell = getattr(self, 'min_dwell', 0.1)
            precision = getattr(self, 'precision', 0.1)
            
            dt_raw = washout_ms / num_analytes
            dt_ms = max(min_dwell, round(dt_raw / precision) * precision)
            target_cycle_time_ms = dt_ms * num_analytes
            
            # Rep-Rate dynamic sweep at 5.0% RSD oversampling level
            target_rsd = 5.0
            rr_prec = getattr(self, 'laser_rr_prec', 1.0)
            step = rr_prec if rr_prec > 0 else 1.0
            
            # Start from f_start = 2.0 / (washout_ms / 1000.0)
            f_start = 2.0 / (washout_ms / 1000.0) if washout_ms > 0 else 1.0
            max_allowed_rr = getattr(self, 'max_rr', 10000.0)
            if hasattr(self, 'spin_init_rr'):
                max_allowed_rr = self.spin_init_rr.maximum
                
            curr_f = f_start
            found_f = None
            while curr_f <= max_allowed_rr + 1e-9:
                rsd = Logic.calculate_lockwood_rsd(curr_f, washout_ms, dt_ms)
                if rsd <= target_rsd:
                    found_f = curr_f
                    break
                curr_f += step
                
            if found_f is None:
                found_f = max_allowed_rr
                
            # Stage Speed (µm/s)
            cycle_time_s = target_cycle_time_ms / 1000.0
            suggested_speed = spot_size / cycle_time_s if cycle_time_s > 0 else 0.0
            
            # Dosage
            suggested_dosage = found_f * cycle_time_s
            
            # Overlap
            suggested_overlap_um = spot_size - (suggested_speed / found_f) if found_f > 0 else spot_size
            suggested_overlap_pct = (suggested_overlap_um / spot_size) * 100.0 if spot_size > 0 else 0.0
            
            # Display dialog
            dlg = QDialog(self)
            dlg.setWindowTitle("Suggested PreScan Parameters")
            dlg.resize(300, 200)
            dlg_layout = QVBoxLayout(dlg)
            
            form = QFormLayout()
            
            lbl_rr = QLabel(f"<b>{found_f:.3f} Hz</b>")
            lbl_speed = QLabel(f"<b>{suggested_speed:.3f} µm/s</b>")
            lbl_overlap = QLabel(f"<b>{suggested_overlap_um:.3f} µm</b>")
            lbl_dosage = QLabel(f"<b>{suggested_dosage:.2f} pulses/pixel</b>")
            lbl_cycle = QLabel(f"<b>{target_cycle_time_ms:.3f} ms</b>")
            lbl_dwell = QLabel(f"<b>{dt_ms:.3f} ms</b>")
            
            form.addRow("Suggested Rep-Rate:", lbl_rr)
            form.addRow("Suggested Stage Speed:", lbl_speed)
            form.addRow("Suggested Overlap:", lbl_overlap)
            form.addRow("Suggested Dosage:", lbl_dosage)
            form.addRow("Target Cycle Time:", lbl_cycle)
            form.addRow("Dwell Time per Analyte:", lbl_dwell)
            
            dlg_layout.addLayout(form)
            
            # Close button
            btn_box = QHBoxLayout()
            btn_box.addStretch()
            btn_close = QPushButton("Close")
            btn_close.clicked.connect(lambda: dlg.accept())
            btn_box.addWidget(btn_close)
            dlg_layout.addLayout(btn_box)
            
            dlg.exec_()
        except Exception as e:
            IoLog.error(f"iolite Optimiser: Error calculating PreScan suggestions: {e}")
            QMessageBox.critical(self, "Calculation Error", f"Failed to compute suggestions: {e}")

    # --- Pulse Train Simulator Methods ---
    def init_pulse_train_tab(self):
        IoLog.information("iolite Optimiser: init_pulse_train_tab starting...")
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()
        
        # --- LEFT COLUMN (CONTROLS) ---
        self.pulse_scroll = QScrollArea()
        self.pulse_scroll.setWidgetResizable(True)
        self.pulse_scroll.setFrameShape(QScrollArea.NoFrame)
        self.pulse_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        
        scroll_content = QWidget()
        l_pulse = QVBoxLayout(scroll_content)
        l_pulse.setContentsMargins(5, 2, 5, 2)
        l_pulse.setSpacing(5)
        
        # Settings Group
        grp_settings = QGroupBox("")
        v_settings = QVBoxLayout()
        v_settings.setContentsMargins(5, 5, 5, 5)
        v_settings.setSpacing(5)
        
        lbl_settings_title = QLabel("Simulation Settings")
        lbl_settings_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_settings_title.setAlignment(Qt.AlignCenter)
        v_settings.addWidget(lbl_settings_title)
        
        grid_settings = QGridLayout()
        grid_settings.setSpacing(5)
        
        self.spin_pulse_bg = QDoubleSpinBox()
        self.spin_pulse_bg.setRange(0, 9999)
        self.spin_pulse_bg.setValue(5.0)
        self.spin_pulse_bg.setSuffix(" s")
        grid_settings.addWidget(QLabel("BG Time:"), 0, 0)
        grid_settings.addWidget(self.spin_pulse_bg, 0, 1)
        
        self.spin_pulse_sig = QDoubleSpinBox()
        self.spin_pulse_sig.setRange(0, 9999)
        self.spin_pulse_sig.setValue(30.0)
        self.spin_pulse_sig.setSuffix(" s")
        grid_settings.addWidget(QLabel("Signal Duration:"), 0, 2)
        grid_settings.addWidget(self.spin_pulse_sig, 0, 3)
        
        self.cmb_pulse_shape = QComboBox()
        self.cmb_pulse_shape.addItems([
            "Real Composite Peak",
            "Model Washout Peak (Lognormal)",
            "Sawtooth Approximation"
        ])
        has_composite = (hasattr(self, 'spr_raw_results_df') and 
                         self.spr_raw_results_df is not None and 
                         not self.spr_raw_results_df.empty)
        default_shape = "Real Composite Peak" if has_composite else "Model Washout Peak (Lognormal)"
        self.cmb_pulse_shape.setCurrentText(default_shape)
        self.cmb_pulse_shape.currentIndexChanged.connect(lambda: self.pulse_debounce_timer.start())
        
        grid_settings.addWidget(QLabel("Pulse Shape:"), 1, 0)
        grid_settings.addWidget(self.cmb_pulse_shape, 1, 1, 1, 3)
        
        v_settings.addLayout(grid_settings)
        grp_settings.setLayout(v_settings)
        l_pulse.addWidget(grp_settings)
        
        # Overrides Group
        self.grp_pulse_override = QGroupBox("")
        v_override = QVBoxLayout()
        v_override.setContentsMargins(5, 5, 5, 5)
        v_override.setSpacing(5)
        
        lbl_override_title = QLabel("Overrides / Manual Adjustment")
        lbl_override_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_override_title.setAlignment(Qt.AlignCenter)
        v_override.addWidget(lbl_override_title)
        
        grid_override = QGridLayout()
        grid_override.setSpacing(5)
        
        self.btn_pulse_reset = QPushButton("Reset to Optimum")
        self.btn_pulse_reset.clicked.connect(self.reset_pulse_params_to_optimum)
        grid_override.addWidget(self.btn_pulse_reset, 0, 0, 1, 4)
        
        self.spin_pulse_rr = QDoubleSpinBox()
        self.spin_pulse_rr.setRange(0.0001, 100000)
        self.spin_pulse_rr.setDecimals(4)
        self.spin_pulse_rr.setSuffix(" Hz")
        grid_override.addWidget(QLabel("Rep Rate:"), 1, 0)
        grid_override.addWidget(self.spin_pulse_rr, 1, 1)
        
        self.spin_pulse_at = QDoubleSpinBox()
        self.spin_pulse_at.setRange(0.001, 10000)
        self.spin_pulse_at.setDecimals(4)
        self.spin_pulse_at.setSuffix(" ms")
        grid_override.addWidget(QLabel("Acq Time:"), 1, 2)
        grid_override.addWidget(self.spin_pulse_at, 1, 3)
        
        self.spin_pulse_count = QDoubleSpinBox()
        self.spin_pulse_count.setRange(0.001, 100000)
        self.spin_pulse_count.setDecimals(3)
        grid_override.addWidget(QLabel("Pulses / Acq:"), 2, 0)
        grid_override.addWidget(self.spin_pulse_count, 2, 1)
        
        self.spin_pulse_washout = QDoubleSpinBox(self)
        self.spin_pulse_washout.setRange(0.01, 50000.0)
        self.spin_pulse_washout.setDecimals(2)
        self.spin_pulse_washout.setSuffix(" ms")
        grid_override.addWidget(QLabel("Washout:"), 2, 2)
        grid_override.addWidget(self.spin_pulse_washout, 2, 3)
        
        # Connect signals for auto-calc and debounce timer
        self.spin_pulse_bg.valueChanged.connect(lambda: self.pulse_debounce_timer.start())
        self.spin_pulse_sig.valueChanged.connect(lambda: self.pulse_debounce_timer.start())
        
        self.spin_pulse_rr.valueChanged.connect(self._calc_pulse_params_from_rr)
        self.spin_pulse_at.valueChanged.connect(self._calc_pulse_params_from_at)
        self.spin_pulse_count.valueChanged.connect(self._calc_pulse_params_from_count)
        self.spin_pulse_washout.valueChanged.connect(lambda: self.pulse_debounce_timer.start())
        
        self.spin_pulse_rr.valueChanged.connect(lambda: self.pulse_debounce_timer.start())
        self.spin_pulse_at.valueChanged.connect(lambda: self.pulse_debounce_timer.start())
        self.spin_pulse_count.valueChanged.connect(lambda: self.pulse_debounce_timer.start())
        
        v_override.addLayout(grid_override)
        self.grp_pulse_override.setLayout(v_override)
        l_pulse.addWidget(self.grp_pulse_override)
        
        l_pulse.addStretch()
        
        self.pulse_scroll.setWidget(scroll_content)
        self.pulse_scroll.setFixedWidth(460)
        left_layout.addWidget(self.pulse_scroll)
        
        # --- RIGHT COLUMN (PLOT) ---
        
        # Plot Controls
        h_ctrl = QHBoxLayout()
        
        self.combo_theme_pulse = QComboBox()
        self.combo_theme_pulse.addItems(["Auto", "Dark", "Light"])
        self.combo_theme_pulse.setCurrentText(self.persistent_settings.get('theme', 'Auto'))
        self.combo_theme_pulse.currentTextChanged.connect(self.apply_theme)
        
        self.chk_norm_pulse = QCheckBox("Normalize")
        self.chk_norm_pulse.setChecked(False)
        self.chk_norm_pulse.toggled.connect(lambda: self.update_pulse_plot(preserve_zoom=True))
        
        self.chk_y_zoom_pulse = QCheckBox("Pan / Zoom Y")
        self.chk_y_zoom_pulse.toggled.connect(self._on_pulse_y_zoom_toggled)
        
        self.chk_rescale_pulse = QCheckBox("Auto-Rescale Y")
        self.chk_rescale_pulse.setChecked(True)
        self.chk_rescale_pulse.toggled.connect(self._on_pulse_rescale_toggled)
        
        self.chk_show_bg = QCheckBox("Show Background")
        self.chk_show_bg.setChecked(True)
        self.chk_show_bg.toggled.connect(lambda: self.update_pulse_plot(preserve_zoom=False))
        
        h_ctrl.addWidget(QLabel("Theme:"))
        h_ctrl.addWidget(self.combo_theme_pulse)
        h_ctrl.addWidget(self.chk_norm_pulse)
        h_ctrl.addWidget(self.chk_y_zoom_pulse)
        h_ctrl.addWidget(self.chk_rescale_pulse)
        h_ctrl.addWidget(self.chk_show_bg)
        h_ctrl.addStretch()
        
        right_layout.addLayout(h_ctrl)
        
        # Canvas
        self.pulse_figure = Figure(figsize=(5, 4), dpi=100)
        self.pulse_canvas = FigureCanvas(self.pulse_figure)
        self.pulse_canvas.mpl_connect('resize_event', self._on_plot_resize)
        self.pulse_canvas.mpl_connect('scroll_event', self.on_zoom)
        self.pulse_canvas.mpl_connect('button_press_event', self.on_press)
        self.pulse_canvas.mpl_connect('button_release_event', self.on_release)
        self.pulse_canvas.mpl_connect('motion_notify_event', self.on_drag)
        self.pulse_canvas.mpl_connect('pick_event', self.on_pulse_pick)
        
        self.pulse_legend_map = {} # Mapping for interactive legend
        
        right_layout.addWidget(self.pulse_canvas, 1) # Set stretch to 1 to give plot more space
        
        # RSD Stats Table
        self.lbl_pulse_table_title = QLabel("<b>Pulse Simulation Stability (Signal Area ± 1s crop)</b>")
        self.lbl_pulse_table_title.setAlignment(Qt.AlignCenter)
        right_layout.addWidget(self.lbl_pulse_table_title)
        
        # Use the custom CopyableTableWidget if available, or fallback to standard
        try:
            self.pulse_table = CopyableTableWidget()
        except:
            self.pulse_table = QTableWidget()
            
        self.pulse_table.setColumnCount(6)
        self.pulse_table.setHorizontalHeaderLabels([
            "Channel", "Mean Intensity", "RSD (%)", 
            "Ideal Match (%)", "Max Error (%)", "Duty Cycle (%)"
        ])
        self.pulse_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        # self.pulse_table.setFixedHeight(120) # REMOVED: Allow it to stretch
        self.pulse_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.pulse_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.pulse_table.setAlternatingRowColors(True)
        right_layout.addWidget(self.pulse_table, 1) # Give it equal stretch with the canvas
        
        main_layout.addLayout(left_layout, 1)
        main_layout.addLayout(right_layout, 4)
        
        self.tab_pulse.setLayout(main_layout)

    def reset_pulse_params_to_optimum(self):
        # Pull from last_sync (Optimiser Tab results) or Optimiser tab inputs
        if hasattr(self, 'last_sync') and self.last_sync:
            rr = self.last_sync.get('rr_actual', 10.0)
            at_ms = self.last_sync.get('at_actual_s', 0.1) * 1000.0
            pulses = self.last_sync.get('actual_pulses', 1.0)
        else:
            rr = self.spin_init_rr.value if hasattr(self, 'spin_init_rr') else 10.0
            at_ms = getattr(self, 'detected_at_ms', 10.0) if getattr(self, 'detected_at_ms', None) is not None else 10.0
            pulses = rr * (at_ms / 1000.0)
            
        # Update spinboxes (Block signals to avoid double-triggers, then manually start timer)
        self.spin_pulse_rr.blockSignals(True)
        self.spin_pulse_at.blockSignals(True)
        self.spin_pulse_count.blockSignals(True)
        self.spin_pulse_washout.blockSignals(True)
        
        self.spin_pulse_rr.setValue(rr)
        self.spin_pulse_at.setValue(at_ms)
        self.spin_pulse_count.setValue(pulses)
        self.spin_pulse_washout.setValue(self.spin_wash.value)
        
        self.spin_pulse_rr.blockSignals(False)
        self.spin_pulse_at.blockSignals(False)
        self.spin_pulse_count.blockSignals(False)
        self.spin_pulse_washout.blockSignals(False)
        
        # Trigger simulation refresh
        self.pulse_debounce_timer.start()

    def _calc_pulse_params_from_rr(self):
        self.spin_pulse_count.blockSignals(True)
        try:
            rr = self.spin_pulse_rr.value
            at_s = self.spin_pulse_at.value / 1000.0
            self.spin_pulse_count.setValue(rr * at_s)
        except: pass
        self.spin_pulse_count.blockSignals(False)

    def _calc_pulse_params_from_at(self):
        self.spin_pulse_count.blockSignals(True)
        try:
            rr = self.spin_pulse_rr.value
            at_s = self.spin_pulse_at.value / 1000.0
            self.spin_pulse_count.setValue(rr * at_s)
        except: pass
        self.spin_pulse_count.blockSignals(False)

    def _calc_pulse_params_from_count(self):
        self.spin_pulse_rr.blockSignals(True)
        try:
            p = self.spin_pulse_count.value
            at_s = self.spin_pulse_at.value / 1000.0
            if at_s > 0:
                rr = p / at_s
                self.spin_pulse_rr.setValue(rr)
        except: pass
        self.spin_pulse_rr.blockSignals(False)

    def run_pulse_simulation(self):
        try:
            IoLog.information("iolite Optimiser: run_pulse_simulation called")
            # 1. Gather Inputs
            if not hasattr(self, 'last_sync') or self.last_sync is None:
                IoLog.warning("Pulse Train: No sync data found yet. Using current or default parameters.")
                # If override checked, we can proceed with just defaults if logic handles None sync

            # Get Selected Isotope (for composite peak)
            try:
                iso = self._get(self.cmb_spr_iso, 'currentText')
                IoLog.information(f"iolite Optimiser: Selected ISO: {iso}")
            except: iso = None
            
            shape_pref = self._get(self.cmb_pulse_shape, 'currentText')
            
            if shape_pref == "Real Composite Peak":
                if not iso:
                    IoLog.warning("iolite Optimiser: No isotope selected")
                    return
                # Get Analyzed Isotopes
                if self.spr_raw_results_df is None or self.spr_raw_results_df.empty:
                     IoLog.warning("iolite Optimiser: No SPR results available")
                     if hasattr(self, 'lbl_result'): self.lbl_result.setText("No SPR analysis available. Please run SPR Analysis first.")
                     return
            else:
                if not iso:
                    iso = "Model Isotope"

            # Get Timing (Always use spinbox values)
            rr = self.spin_pulse_rr.value
            at_ms = self.spin_pulse_at.value
            at_s = at_ms / 1000.0
            pulses = self.spin_pulse_count.value
            
            IoLog.information(f"iolite Optimiser: Simulation Params - RR: {rr}, AT: {at_ms}ms, Pulses: {pulses}")

            try:
                bg_s = self.spin_pulse_bg.value
                sig_s = self.spin_pulse_sig.value
            except Exception as e:
                 IoLog.error(f"iolite Optimiser: Failed to get spinbox values: {e}")
                 return

            # 2. Generate Data with Sequential Sampling (V3)
            self.pulse_results = {} # Dict of {iso: (theo_df, meas_df)}
            
            # --- Build Channel Specs ---
            channel_specs = []
            
            if hasattr(self, 'optimised_results') and self.optimised_results:
                # Use optimized order
                is_simultaneous = (self.icp_tech in ["Multi-Collector", "TOF"])
                current_offset_s = 0.0
                
                for res in self.optimised_results:
                    ch_name = str(res.get('Isotope', 'Unknown'))
                    dwell_ms = res.get('Final Dwell (ms)', 10.0)
                    
                    try:
                        if 'Signal CPS' in res:
                             scale_val = float(res['Signal CPS'])
                        elif 'Resultant SNR' in res:
                             scale_val = float(res['Resultant SNR'])
                        else:
                             scale_val = 1.0
                    except:
                        scale_val = 1.0
                    
                    ratio = 1.0
                    if hasattr(self, 'last_sync') and self.last_sync:
                         orig_at = self.last_sync.get('at_actual_s', 0.001) * 1000.0
                         if orig_at > 0: ratio = at_ms / orig_at
                    
                    dwell_s = (dwell_ms * ratio) / 1000.0
                    
                    channel_specs.append({
                        'name': ch_name,
                        'dwell': dwell_s,
                        'offset': 0.0 if is_simultaneous else current_offset_s,
                        'scale': scale_val
                    })
                    
                    if not is_simultaneous:
                        current_offset_s += dwell_s
                    
            else:
                # Fallback: Simultaneous or Equal Sequential
                if hasattr(self, 'spr_df') and self.spr_df is not None:
                    cols = [c for c in self.spr_df.columns if c != 'Time']
                    for c in cols:
                        channel_specs.append({
                            'name': c,
                            'dwell': at_s,
                            'offset': 0.0
                        })
            
            if not channel_specs:
                channel_specs.append({
                    'name': iso,
                    'dwell': at_s,
                    'offset': 0.0
                })
            
            IoLog.information(f"iolite Optimiser: Channel Specs (First 3): {channel_specs[:3]}")

            # Generate peak based on preference
            target_washout = self.spin_pulse_washout.value
            target_washout_s = target_washout / 1000.0

            if shape_pref == "Sawtooth Approximation":
                t_scaled = np.linspace(0.0, target_washout_s, 1000)
                y_vals = 1.0 - (t_scaled / target_washout_s)
                import pandas as pd
                composite_df = pd.DataFrame({
                    'Relative Time (s)': t_scaled,
                    'Normalised Intensity': y_vals
                })
            elif shape_pref == "Model Washout Peak (Lognormal)":
                t_vals = np.linspace(1e-6, 4.0, 1000)
                sigma = 0.5
                mu = 0.0
                s2pi = np.sqrt(2.0 * np.pi)
                y_vals = 1.0 / (t_vals * sigma * s2pi) * np.exp(-0.5 * ((np.log(t_vals) - mu) / sigma) ** 2)
                y_vals = y_vals / y_vals.max()
                
                # Scale the entire time axis exactly to [0.0, target_washout_s]
                t_scaled = (t_vals - t_vals[0]) / (t_vals[-1] - t_vals[0]) * target_washout_s
                import pandas as pd
                composite_df = pd.DataFrame({
                    'Relative Time (s)': t_scaled,
                    'Normalised Intensity': y_vals
                })
            else: # Real Composite Peak
                excluded = self.spr_excluded_peaks.get(iso, set())
                filtered_df = self.spr_raw_results_df[~self.spr_raw_results_df['Peak Index'].isin(excluded)]
                composite_df = Logic.generate_composite_peak(self.spr_df, iso, filtered_df)
                
                if composite_df is not None:
                    lvl = getattr(self, 'active_washout_level', 0.01)
                    comp_t = composite_df['Relative Time (s)'].values
                    comp_y = composite_df['Normalised Intensity'].values
                    above = np.where(comp_y >= lvl)[0]
                    if len(above) >= 2:
                        orig_width_ms = (comp_t[above[-1]] - comp_t[above[0]]) * 1000.0
                    else:
                        orig_width_ms = (comp_t.max() - comp_t.min()) * 1000.0
                    
                    if orig_width_ms > 0:
                        scale_factor = target_washout / orig_width_ms
                        if scale_factor > 0:
                            composite_df = composite_df.copy()
                            composite_df['Relative Time (s)'] = composite_df['Relative Time (s)'] * scale_factor
                        IoLog.information(f"iolite Optimiser: Scaling composite peak by {scale_factor:.4f} (Original {lvl*100}% width: {orig_width_ms:.2f} ms, Target: {target_washout:.2f} ms)")
            
            if composite_df is not None:
                try:
                    IoLog.information(f"Calling Logic.generate_pulse_train_v3 for {iso}...")
                    
                    # Calculate dt from the composite dataframe (scaled SPR resolution)
                    comp_t = composite_df['Relative Time (s)'].values
                    dt_scaled = np.median(np.diff(comp_t)) if len(comp_t) > 1 else 1e-3
                    
                    # Ensure dt_s is fine enough to:
                    # 1. Resolve the shape of the peak itself without aliasing (at least 80 steps per peak duration)
                    # 2. Resolve the minimum active dwell (at least 10 steps per dwell)
                    pulse_duration = comp_t.max() - comp_t.min() if len(comp_t) > 0 else 0.1
                    active_dwells = [spec['dwell'] for spec in channel_specs if spec['dwell'] > 0]
                    min_dwell = min(active_dwells) if active_dwells else 0.01
                    # Lock to a fine resolution of 1e-5 s (10 µs) for maximum integration accuracy
                    dt_s_calc = 1e-5
                    
                    IoLog.information(f"iolite Optimiser: Using high-resolution simulator step dt_s = {dt_s_calc:.6f} s")
                    
                    # Use v3
                    theo_df, channel_results = Logic.generate_pulse_train_v3(
                        composite_df, rr, pulses, at_s, channel_specs, background_s=bg_s, signal_s=sig_s, dt_s=dt_s_calc
                    )
                    
                    if theo_df is not None and channel_results:
                        # Store in compatible format for plotting
                        # pulse_results expects {iso: (theo_df, meas_df)}
                        # v3 returns ONE theo_df and a dict of meas_dfs
                        
                        for ch_name, meas_df in channel_results.items():
                             self.pulse_results[ch_name] = (theo_df, meas_df)
                            
                        IoLog.information(f"Simulation valid. Channels: {len(channel_results)}")
                        
                except Exception as e:
                    IoLog.error(f"Logic.generate_pulse_train FAILED for {iso}: {e}")
            
            if not self.pulse_results:
                 IoLog.warning("iolite Optimiser: Simulation returned NO results")
            
            # 3. Plot
            self.update_pulse_plot()
            
            # Store parameters of successful run
            is_override_run = False
            if hasattr(self, 'grp_pulse_override') and self.grp_pulse_override.isChecked():
                is_override_run = True

            self._last_simulated_params = {
                'rr': rr,
                'at_ms': at_ms,
                'pulses': pulses,
                'bg_s': bg_s,
                'sig_s': sig_s,
                'iso': iso,
                'washout': target_washout,
                'is_override': is_override_run
            }
            
        except Exception as e:
            IoLog.error(f"iolite Optimiser: run_pulse_simulation CRASHED: {e}")
            IoLog.error(traceback.format_exc())


    def update_pulse_plot(self, preserve_zoom=False):
        IoLog.information("iolite Optimiser: update_pulse_plot called")
        if not hasattr(self, 'pulse_results') or not self.pulse_results: 
             IoLog.warning("iolite Optimiser: No pulse_results to plot")
             return
        
        try:
            # Theme Handling
            is_dark = self.system_is_dark
            try:
                theme_sel = self._get(self.combo_theme_pulse, 'currentText')
                if theme_sel == "Dark": is_dark = True
                elif theme_sel == "Light": is_dark = False
            except: pass
            
            bg_color = '#2b2b2b' if is_dark else 'white'
            fg_color = 'white' if is_dark else 'black'
            
            self.pulse_figure.clear()
            self.pulse_figure.patch.set_facecolor(bg_color)
            
            ax = self.pulse_figure.add_subplot(111)
            ax.set_facecolor(bg_color)
            
            norm = self.chk_norm_pulse.isChecked()
            
            # Setup Color Cycler
            # Use a qualitative map like tab10
            prop_cycle = plt.rcParams['axes.prop_cycle']
            colors = prop_cycle.by_key()['color']
            
            # Plot loop
            # 1. Plot Theoretical (Using first available, as they are identical)
            # Just extract one theo_df
            first_iso = next(iter(self.pulse_results))
            theo_df, _ = self.pulse_results[first_iso]
            
            t_theo = theo_df['Time'].values
            y_theo = theo_df['Intensity'].values
            
            # 2. Calc Global Max first (needed for scaling Theo)
            all_intensities = []
            for _, (_, df) in self.pulse_results.items():
                all_intensities.extend(df['Intensity'].values)
            
            global_max_meas = 1.0
            if all_intensities:
                global_max_meas = np.max(all_intensities)
            if global_max_meas <= 0: global_max_meas = 1.0
            
            own_theo_max = np.max(y_theo) if len(y_theo) > 0 else 1.0
            if own_theo_max <= 0: own_theo_max = 1.0

            # --- CROP LOGIC ---
            show_bg = self.chk_show_bg.isChecked()
            bg_s = self.spin_pulse_bg.value
            sig_s = self.spin_pulse_sig.value
            t_min_crop = bg_s + 1.0
            t_max_crop = bg_s + sig_s - 1.0

            if not show_bg:
                # Normalise to max within signal window for overlay
                mask_theo = (t_theo >= t_min_crop) & (t_theo <= t_max_crop)
                theo_sig_max = np.max(y_theo[mask_theo]) if any(mask_theo) else own_theo_max
                y_theo = y_theo / (theo_sig_max if theo_sig_max > 0 else 1)
            elif norm:
                y_theo = y_theo / own_theo_max
            else:
                # If raw, normalize Theo to 1, then scale to match data max so it is visible
                y_theo = (y_theo / own_theo_max) * global_max_meas

            # Plot Theoretical once (Faint background)
            line_theo = ax.plot(t_theo, y_theo, label="Theoretical (Pulse Stream)", color='gray', alpha=0.3, lw=1)[0]
            if hasattr(self, 'pulse_hidden_channels') and "Theoretical (Pulse Stream)" in self.pulse_hidden_channels:
                line_theo.set_visible(False)
            
            # 3. Plot Measured for each channel

            for idx, (iso, (_, meas_df)) in enumerate(self.pulse_results.items()):

                t_meas = meas_df['Time'].values
                y_meas = meas_df['Intensity'].values
                
                if not show_bg:
                    # LOCAL NORM in signal window for shape comparison (overlay)
                    mask_meas = (t_meas >= t_min_crop) & (t_meas <= t_max_crop)
                    meas_sig_max = np.nanmax(y_meas[mask_meas]) if any(mask_meas) else (np.max(y_meas) if len(y_meas)>0 else 1)
                    y_meas = y_meas / (meas_sig_max if meas_sig_max > 0 else 1)
                elif norm:
                    # Normalize against OWN max so shapes can be compared visually
                    own_max = np.max(y_meas) if len(y_meas) > 0 else 1
                    if own_max > 0:
                        y_meas = y_meas / own_max
                
                color = colors[idx % len(colors)]
                
                # Plot Measured (Connected Lines with Markers)
                # Now that m_times are aligned to cycle, we can just plot them directly.
                
                line_meas = ax.plot(t_meas, y_meas, label=f"{iso}", color=color, marker='o', markersize=4, lw=2)[0]
                if hasattr(self, 'pulse_hidden_channels') and iso in self.pulse_hidden_channels:
                    line_meas.set_visible(False)
            
            ax.set_xlabel("Time (s)", color=fg_color)
            if not show_bg:
                t_low = t_min_crop - 2.5
                t_high = t_max_crop + 2.5
                ax.set_xlim(t_low, t_high)
                
                # --- DYNAMIC Y-SCALING to fill 90% ---
                # Determine min/max based ONLY on the stable measured signal plateau [t_min_crop, t_max_crop]
                # We EXCLUDE the theoretical trace here because it pulses down to 0, which would 
                # prevent the zoom from focusing on the measurement stability/shape.
                all_signal_plateau_y = []
                
                # Measured (normalized to 1.0 locally)
                for i_iso, (iso, (_, m_df)) in enumerate(self.pulse_results.items()):
                    tm = m_df['Time'].values
                    ym = m_df['Intensity'].values
                    
                    # Local max for this isotope (in stable region)
                    mask_m_stable = (tm >= t_min_crop) & (tm <= t_max_crop)
                    if not any(mask_m_stable): continue
                    
                    m_max = np.nanmax(ym[mask_m_stable])
                    if m_max <= 0: m_max = 1.0
                    
                    all_signal_plateau_y.extend(ym[mask_m_stable] / m_max)
                
                if all_signal_plateau_y:
                    y_min_data = np.nanmin(all_signal_plateau_y)
                    y_max_data = np.nanmax(all_signal_plateau_y) 
                    
                    data_range = y_max_data - y_min_data
                    if data_range <= 0: data_range = 0.001
                    
                    H = data_range / 0.9
                    
                    y_top = y_max_data + 0.05 * H
                    y_bottom = y_top - H
                    
                    ax.autoscale(False)
                    ax.set_ylim(y_bottom, y_top)

                ax.set_ylabel("Normalized Intensity (Signal Overlay)", color=fg_color)
            else:
                ax.set_ylabel("Normalized Intensity" if norm else "Intensity", color=fg_color)
            
            # Legend with interactivity
            self.pulse_legend_map = {}
            handles, labels = ax.get_legend_handles_labels()
            
            import matplotlib.transforms as mtransforms
            leg = ax.legend(handles, labels, loc='lower center', 
                            bbox_to_anchor=(0.5, 1.0), 
                            bbox_transform=mtransforms.blended_transform_factory(ax.figure.transFigure, ax.transAxes),
                            borderaxespad=0.5, ncol=min(10, len(handles)), frameon=True, fontsize='small',
                            handlelength=1.5, handletextpad=0.7, columnspacing=1.5)
            
            leg.get_frame().set_alpha(0.0) # Transparent frame
            
            for legline, handle in zip(leg.get_lines(), handles):
                legline.set_picker(5) # 5 pts tolerance
                self.pulse_legend_map[legline] = handle
                if not handle.get_visible():
                    legline.set_alpha(0.2)
            
            for text in leg.get_texts():
                text.set_color(fg_color)

            
            ax.tick_params(colors=fg_color)
            for spine in ax.spines.values():
                spine.set_edgecolor(fg_color)
            
            # Use EngFormatter for non-normalized
            if not norm:
                 ax.yaxis.set_major_formatter(EngFormatter(sep=" "))
            
            if self.chk_rescale_pulse.isChecked() and not preserve_zoom and show_bg:
                 self.rescale_to_visible(ax=ax, canvas=None)
            elif preserve_zoom:
                 pass 
                 
            self._update_smart_margins(self.pulse_canvas)
            self.pulse_canvas.draw()
            IoLog.information("iolite Optimiser: pulse_canvas.draw() called successfully")
            
            # --- Update RSD Stats Table ---
            self.pulse_table.setRowCount(0)
            self.pulse_table.setRowCount(len(self.pulse_results))
            
            try:
                bg_s = self.spin_pulse_bg.value
                sig_s = self.spin_pulse_sig.value
                t_min = bg_s + 1.0
                t_max = bg_s + sig_s - 1.0
                
                # Helper for SI Formatting
                def _si_pulse(val):
                    if val >= 1e9: return f"{val/1e9:.2f} G"
                    if val >= 1e6: return f"{val/1e6:.2f} M"
                    if val >= 1e3: return f"{val/1e3:.2f} k"
                    return f"{val:.2f}"
                
                # Calculate Theoretical Mean in signal window for Ideal Match
                t_theo = theo_df['Time'].values
                y_theo = theo_df['Intensity'].values
                # We need to re-apply the same scaling/normalization to y_theo that was used for plotting
                # but it's simpler to just use raw intensities for the ratio.
                # Logic: Ideal Match compares the 'area' or 'average' of measured vs theory.
                theo_mask = (t_theo >= t_min) & (t_theo <= t_max)
                y_theo_crop = y_theo[theo_mask]
                theo_mean = np.nanmean(y_theo_crop) if len(y_theo_crop) > 0 else 1.0
                if theo_mean <= 0: theo_mean = 1.0
                
                # Total Cycle Time (Acq Time)
                cycle_time_ms = self.spin_pulse_at.value
                if cycle_time_ms <= 0: cycle_time_ms = 1.0

                for i, (iso, (_, meas_df)) in enumerate(self.pulse_results.items()):
                    t_meas = meas_df['Time'].values
                    y_meas = meas_df['Intensity'].values
                    
                    mask = (t_meas >= t_min) & (t_meas <= t_max)
                    y_crop = y_meas[mask]
                    
                    if len(y_crop) > 1:
                        m = np.nanmean(y_crop)
                        s = np.nanstd(y_crop)
                        rsd = (s / m * 100.0) if m > 0 else 0.0
                        
                        # Duty Cycle: Dwell / AT
                        # We need to find the dwell and scale for this isotope
                        dwell_ms = 1.0
                        ch_scale = 1.0
                        if hasattr(self, 'optimised_results') and self.optimised_results:
                             for res in self.optimised_results:
                                 if res.get('Isotope') == iso:
                                     dwell_ms = res.get('Final Dwell (ms)', 1.0)
                                     # Get scaling factor (same logic as run_pulse_simulation)
                                     try:
                                         if 'Signal CPS' in res: ch_scale = float(res['Signal CPS'])
                                         elif 'Resultant SNR' in res: ch_scale = float(res['Resultant SNR'])
                                     except: pass
                                     break
                        
                        # Ideal Match: How close is the measured mean to the theoretical mean (scaled)?
                        # Theoretical mean is based on 0-1 pulse intensity. 
                        # Measured is based on CPS or SNR (the ch_scale).
                        scaled_theo_mean = theo_mean * ch_scale
                        ideal_match = (m / scaled_theo_mean * 100.0) if scaled_theo_mean > 0 else 0.0
                        
                        # Max Error: Largest deviation from scaled theoretical mean
                        max_err_abs = np.nanmax(np.abs(y_crop - scaled_theo_mean))
                        max_err_pct = (max_err_abs / scaled_theo_mean * 100.0) if scaled_theo_mean > 0 else 0.0
                        
                        duty_cycle = (dwell_ms / cycle_time_ms * 100.0)

                        def _it(txt):
                            it = QTableWidgetItem(txt)
                            it.setTextAlignment(Qt.AlignCenter)
                            return it
                            
                        self.pulse_table.setItem(i, 0, _it(iso))
                        self.pulse_table.setItem(i, 1, _it(_si_pulse(m)))
                        self.pulse_table.setItem(i, 2, _it(f"{rsd:.2f} %"))
                        self.pulse_table.setItem(i, 3, _it(f"{ideal_match:.1f} %"))
                        self.pulse_table.setItem(i, 4, _it(f"{max_err_pct:.1f} %"))
                        self.pulse_table.setItem(i, 5, _it(f"{duty_cycle:.1f} %"))
                    else:
                        for col in range(6):
                            txt = iso if col == 0 else ("N/A (Window too small)" if col == 1 else "-")
                            it = QTableWidgetItem(txt)
                            it.setTextAlignment(Qt.AlignCenter)
                            self.pulse_table.setItem(i, col, it)
            except Exception as e:
                IoLog.error(f"Pulse Table Update Error: {e}")
            
        except Exception as e:
            IoLog.error(f"iolite Optimiser: update_pulse_plot failed: {e}")
            IoLog.error(traceback.format_exc())



    def on_pulse_pick(self, event):
        """
        Handles click events on the pulse plot legend to toggle trace visibility.
        """
        legline = event.artist
        origline = self.pulse_legend_map.get(legline)
        
        if origline:
            vis = not origline.get_visible()
            origline.set_visible(vis)
            
            # Dim the legend line to indicate hidden state
            legline.set_alpha(1.0 if vis else 0.2)
            
            # Track hidden state by label
            label = origline.get_label()
            if not hasattr(self, 'pulse_hidden_channels'):
                self.pulse_hidden_channels = set()
            if vis:
                self.pulse_hidden_channels.discard(label)
            else:
                self.pulse_hidden_channels.add(label)
                
            self.pulse_canvas.draw()

    def _on_pulse_rescale_toggled(self, checked):
        if checked:
            self.chk_y_zoom_pulse.setChecked(False)
            if hasattr(self, 'pulse_figure') and self.pulse_figure.axes:
                 self.rescale_to_visible(ax=self.pulse_figure.axes[0], canvas=self.pulse_canvas)

    def _on_pulse_y_zoom_toggled(self, checked):
        if checked:
            self.chk_rescale_pulse.setChecked(False)

# --- UI SETUP ---
widget = None

def createUIElements():
    # 'ui' is a global object provided by iolite for UI plugins
    global widget
    IoLog.information("iolite Optimiser: createUIElements called")
    
    try:
        widget = ioliteOptimiser()
        widget.setWindowTitle(f"iolite Optimiser - v{VERSION}")
    except Exception as e:
        IoLog.error(f"iolite Optimiser: Error creating widget: {e}")
        IoLog.error(traceback.format_exc())
        return
        
    action = QAction("iolite\nOptimiser", None)
    
    try:
        from iolite.ui import CommonUIPyInterface as CUI
        cui = CUI()
        icon = None
        # Try a few different icon names
        for name in ['sliders', 'chart', 'chart-line', 'analytics', 'settings', 'funnel', 'trophy']:
            try:
                ic = cui.icon(name)
                if not ic.isNull():
                    icon = ic
                    break
            except: pass
            
        if icon:
            action.setIcon(icon)
        else:
            action.setIcon(cui.icon('sliders'))
    except Exception as e:
        IoLog.warning(f"iolite Optimiser: Could not set action icon: {e}")
        
    ui.setWidget(widget) # type: ignore
    ui.setAction(action) # type: ignore
    ui.setMenuName(['Tools']) # type: ignore
