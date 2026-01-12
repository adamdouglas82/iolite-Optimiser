#/ Name: iolite Optimiser
#/ Description: Characterises your system's single-pulse washout to calculate the optimal balance of spot size, scan speed, and repetition rate that achieves your target data quality.
#/ Type: UI
#/ Authors: Adam Douglas
#/ Version: 0.9.4
#/ Contact: Adam.Douglas@icpms.com

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
try:
    from iolite.QtGui import (QAction, QSizePolicy, QWidget, QLabel, QDoubleSpinBox, QSpinBox, QCheckBox, 
                              QComboBox, QGridLayout, QHBoxLayout, QVBoxLayout, QGroupBox, QFormLayout, 
                              QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QMenu, 
                              QColorDialog, QDialog, QPushButton, QScrollArea, QSplitter, QFrame, 
                              QLineEdit, QTabWidget, QStackedWidget, QApplication, QPalette, QColor,
                              QListWidget, QListWidgetItem, QTextEdit, QInputDialog, QMessageBox)
    from iolite.QtCore import Qt, QTimer, QSize, QEvent, QUrl
    from iolite import data, IoLog
except ImportError:
    # Fallback for non-iolite environments (VSCode)
    try:
        from PyQt5.QtWidgets import (QAction, QSizePolicy, QWidget, QLabel, QDoubleSpinBox, QSpinBox, QCheckBox, 
                                     QComboBox, QGridLayout, QHBoxLayout, QVBoxLayout, QGroupBox, QFormLayout, 
                                     QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QMenu, 
                                     QColorDialog, QDialog, QPushButton, QScrollArea, QSplitter, QFrame, 
                                     QLineEdit, QTabWidget, QStackedWidget, QApplication,
                                     QListWidget, QListWidgetItem, QTextEdit, QInputDialog, QMessageBox)
        from PyQt5.QtCore import Qt, QTimer, QSize, QEvent, QUrl
        from PyQt5.QtGui import QPalette, QColor
    except ImportError:
        pass

# --- EMBEDDED CONSTANTS ---
VERSION = "0.9.4"

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
        "iCAP MRX": {"type": "Quadrupole", "min_dwell": 0.5, "prec": 0.001},
        "iCAP MTX": {"type": "Quadrupole", "min_dwell": 0.5, "prec": 0.001},
        "Element 2":    {"type": "Sector-Field", "min_dwell": 2.0, "prec": 0.001},
        "Element 2 XR": {"type": "Sector-Field", "min_dwell": 0.1, "prec": 0.001},
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
    },
    "Custom": {
        "Custom Model": {"type": "Quadrupole", "min_dwell": 0.1, "prec": 0.1}
    }
}


MFR_DEFAULT_UNITS = {
    "Agilent": "s",
    "Thermo": "ms",
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

# --- EMBEDDED LOGIC ---
class Logic:
    @staticmethod
    def analyse_washout_peaks(df, isotope, prominence_threshold=100.0, min_distance=1, is_cps=True):
        """
        Analyses washout peaks (Single Pulse Response) for a given isotope.
        Calculates FW0.1M (10% Height) and FW0.01M (1% Height).
        Derived from logic.py
        """
        if find_peaks is None:
            return pd.DataFrame(), {"Error": "Scipy (scipy.signal) is not installed."}

        if isotope not in df.columns or "Time" not in df.columns:
            return pd.DataFrame(), {"Error": f"Columns missing: {isotope} or Time"}

        y_proc = df[isotope].values
        time_proc = df['Time'].values

        # 1. Find Peaks
        try:
            # Pass min_distance (distance in samples)
            peaks, properties = find_peaks(y_proc, prominence=prominence_threshold, distance=max(1, int(min_distance)))
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
                    area_001 = np.trapz(y_proc[idx_start_001:idx_end_001], x=time_proc[idx_start_001:idx_end_001])
                    area_01 = np.trapz(y_proc[max(0, int(np.floor(left_ips_01[i]))):min(len(time_proc), int(np.ceil(right_ips_01[i])))], 
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
    def generate_composite_peak(df, isotope, peaks_df):
        """
        Generates a composite (averaged) peak shape by aligning and normalizing detected peaks.
        Uses asymmetric windows that split the baseline gap between adjacent pulses.
        Aligned such that Relative Time 0 is at the 1% Rise Marker (l001).
        """
        if peaks_df.empty or isotope not in df.columns or "Time" not in df.columns:
            return None

        y_raw = df[isotope].values
        t_raw = df['Time'].values
        
        # Calculate typical sample rate
        dt = np.median(np.diff(t_raw)) if len(t_raw) > 1 else 1.0
        
        # 1. Calculate the core pulse dimensions
        # Alignment is at l001. We want to cover from the start of the buffer to the end of the pulse + buffer.
        max_fw001 = (peaks_df['r001'] - peaks_df['l001']).max()
        
        # 2. Asymmetric baseline buffers (User Request: Show more peak, peak to the left)
        if len(peaks_df) > 1:
            pdf = peaks_df.sort_values('Peak Time (s)')
            gaps = pdf['l001'].values[1:] - pdf['r001'].values[:-1]
            median_gap = max(0, np.median(gaps))
            
            # User Request: Exactly 10% margin on both sides
            b_left = 0.1 * max_fw001
            b_right = 0.20 * max_fw001
            
            # Safety: Ensure we don't clip next peak (40% of gap)
            b_right = min(b_right, 0.4 * median_gap)
        else:
            # Single peak fallbacks
            b_left = 0.1 * max_fw001
            b_right = 0.20 * max_fw001
            
        # Global minimums (0.1ms / 0.1ms) - Fixes "needle" peaks
        b_left = max(b_left, 0.0001)
        b_right = max(b_right, 0.0001)
        
        peak_segments = []
        
        for i, row in peaks_df.iterrows():
            # Align at l001 absolute time
            t_start_abs = row['l001'] - b_left
            t_end_abs = row['r001'] + b_right
            
            # Convert to indices
            idx_start = np.searchsorted(t_raw, t_start_abs)
            idx_end = np.searchsorted(t_raw, t_end_abs, side='right')
            
            if idx_end <= idx_start: continue
            
            seg_y = y_raw[idx_start:idx_end]
            seg_t = t_raw[idx_start:idx_end] - row['l001'] # Center t=0 at l001
            
            # Normalise intensity
            mx = row['Max Intensity']
            if mx > 0:
                seg_y_norm = seg_y / mx
                peak_segments.append((seg_t, seg_y_norm))

        if not peak_segments:
            return None

        # Interpolation to a common grid
        # Grid covers from -b_left to max_pulse_width + b_right
        grid_t = np.arange(-b_left, max_fw001 + b_right + dt, dt)
        sum_y = np.zeros_like(grid_t)
        count_y = np.zeros_like(grid_t)
        
        for seg_t, seg_y in peak_segments:
            interp_y = np.interp(grid_t, seg_t, seg_y, left=0, right=0)
            sum_y += interp_y
            count_y += 1
            
        avg_y = sum_y / count_y if len(peak_segments) > 0 else sum_y
        
        return pd.DataFrame({'Relative Time (s)': grid_t, 'Normalised Intensity': avg_y})

    @staticmethod
    def calculate_constrained_at(inputs):
        """Helper to determine constrained Acquisition Time."""
        fw_s = inputs['washout_ms'] / 1000.0
        bs = inputs['spot_size_um']
        n_val = inputs['pulses_per_pixel']
        rr_max = inputs['max_rr_hz']
        ss_max = inputs['max_speed_um_s']
        is_mc = inputs['icp_technology'] == "Multi-Collector"
        valid_times_s = [d/1000.0 for d in (inputs.get('allowed_dwells') or [])]

        at_from_washout = fw_s
        rr_ideal = n_val / fw_s
        at_from_rr = n_val / rr_max if rr_ideal > rr_max else fw_s
        ss_ideal = bs / fw_s
        at_from_ss = bs / ss_max if ss_ideal > ss_max else fw_s

        min_duty = inputs.get('min_duty_cycle', 0)
        overhead_s = inputs.get('overhead_ms', 0) / 1000.0
        min_dwell_req_s = inputs.get('min_dwell_needed_ms', 0) / 1000.0
        
        at_from_duty = fw_s
        harmonic = 1
        
        # State tracking for reasons why Acq Time was increased beyond washout
        rr_limited = False
        ss_limited = False
        duty_info = None  # Store detail like "10%" or "2x"
        dwell_info = None # Store detail like "50ms" or "2x"

        if at_from_rr > (fw_s + 1e-7):
            rr_limited = True
        
        if at_from_ss > (fw_s + 1e-7):
            ss_limited = True

        # --- HARMONIC SCALING (Hardware-Aware) ---
        # 1. Base Quantum Calculation
        # We start with the washout-compatible base, then account for hardware limits.
        rr_h1 = max(1, math.floor(n_val / fw_s))
        at_h1 = n_val / rr_h1
        p_s = inputs['precision_ms'] / 1000.0
        base_washout = round(at_h1 / p_s) * p_s 
        
        # The 'base' for harmonics must be at least the physical hardware limit
        base = max(base_washout, at_from_rr, at_from_ss)
        
        # 2. Use hardware-limit base for Duty Cycle scaling
        if min_duty > 0 and min_duty < 1.0:
            overhead_s = inputs.get('overhead_ms', 0) / 1000.0
            min_dwell_for_duty = (min_duty * overhead_s) / (1.0 - min_duty)
            required_cycle_for_duty = min_dwell_for_duty + overhead_s
            
            # For MC, we don't force Harmonic multiples (N * Washout). 
            # We just need AT >= Required Duty Cycle Time.
            if is_mc:
                at_from_duty_direct = required_cycle_for_duty
                if required_cycle_for_duty > (base + 1e-7):
                     duty_info = f"{min_duty*100:.0f}%"
            else:
                calc_harmonic = math.ceil((required_cycle_for_duty / base) - 1e-7)
                if calc_harmonic > 1:
                    duty_info = f"{calc_harmonic}x"
                harmonic = max(harmonic, calc_harmonic)
                at_from_duty_direct = 0
        else:
             at_from_duty_direct = 0

        # 3. Use hardware-limit base for Dwell Budget scaling
        if min_dwell_req_s > 0:
            required_cycle = min_dwell_req_s + overhead_s
            if is_mc:
                 at_from_dwell_direct = required_cycle
                 if required_cycle > (base + 1e-7):
                     dwell_info = f"{min_dwell_req_s*1000:.0f} ms"
            else:
                 harmonic_dwell = math.ceil((required_cycle / base) - 1e-7)
                 if harmonic_dwell > 1:
                     dwell_info = f"{harmonic_dwell}x"
                 harmonic = max(harmonic, harmonic_dwell)
                 at_from_dwell_direct = 0
        else:
             at_from_dwell_direct = 0

        # Build combined note for Acq Time increases
        notes = []
        at_reasons = []
        
        if rr_limited:
            at_reasons.append("Rep Rate Limit")
            if n_val > 1.0:
                at_reasons.append(f"Pulses per Dwell Time ({n_val:.1f})")
        if ss_limited:
            at_reasons.append("Stage Speed Limit")
        if duty_info:
            at_reasons.append(f"Duty Cycle ({duty_info})")
        if dwell_info:
            at_reasons.append(f"Dwell Budget ({dwell_info})")

        if at_reasons:
            res_str = at_reasons[0]
            if len(at_reasons) > 1:
                res_str = ", ".join(at_reasons[:-1]) + " and " + at_reasons[-1]
            notes.append(f"Constraint: Acq Time increased for {res_str}")
            
        # Re-verify against allowed list for MC (if applicable)
        if is_mc and valid_times_s:
            # For MC, 'base' refers to chosen Integration Time + Overhead
            for t in sorted(valid_times_s):
                if (t + overhead_s) >= fw_s: 
                    base = t + overhead_s
                    break
        
        at_from_duty = base * harmonic
        
        # Incorporate Direct Requirements (for MC)
        at_needed = max(at_from_washout, at_from_rr, at_from_ss, at_from_duty, at_from_duty_direct, at_from_dwell_direct)
        at_final = at_needed
        if is_mc and valid_times_s:
            at_final = max(valid_times_s) + overhead_s
            for t in sorted(valid_times_s):
                if (t + overhead_s) >= at_needed: 
                    at_final = t + overhead_s
                    break
        
        # --- REFINEMENT: Acq Time must match discrete Rep Rate Steps ---
        optimised_n = n_val
        rr_actual = 0
        
        if is_mc:
            # MC STRATEGY: Time-Driven (Integration Priority)
            # The Acquisition Time is dictated strictly by the Dwell Time + Overhead.
            # We do NOT force it to be N / RR. 
            # The Laser Rep Rate will be matched "as close as possible" in the sync function,
            # and the resulting pulses-per-pixel will float.
            
            at_actual_s = at_needed
            optimised_n = n_val 
            
            # Since we are Time-Driven, we don't calculate an integer RR here for constraints.
            # We rely on calculate_laser_sync to pick the nearest RR.
            
        else:
            # Standard Systems: Fixed Pulses, Floor RR (Pulse-Driven)
            # Refined: Use avoid_gaps logic here so scaling is accurate
            avoid_gaps = inputs.get('avoid_gaps', False)
            rr_theoretical = n_val / at_final
            prec = inputs.get('rr_prec_hz', 1.0)
            if prec <= 0: prec = 1.0
            
            allowed_rr = inputs.get('allowed_rr', None)
            if allowed_rr:
                strategy = "ceil" if avoid_gaps else "floor"
                if strategy == "ceil":
                    valid = [r for r in allowed_rr if r >= rr_theoretical - 1e-9]
                    rr_actual = min(valid) if valid else max(allowed_rr)
                else:
                    valid = [r for r in allowed_rr if r <= rr_theoretical + 1e-9]
                    rr_actual = max(valid) if valid else min(allowed_rr)
            else:
                if avoid_gaps:
                    rr_actual = math.ceil(rr_theoretical / prec) * prec
                else:
                    rr_actual = math.floor(rr_theoretical / prec) * prec
            
            # GUARD: Prevent Zero Division if Rep Rate floors to 0 (e.g. theoretical < precision)
            if rr_actual <= 1e-9:
                 rr_actual = prec if prec > 0 else 1.0
            
            # Acquisition Time must match this RR
            if avoid_gaps:
                # If Avoid Gaps (Overlapping), we prioritize maintaining the Dwell Time (and thus Overlap)
                # over the exact Pulse Count. So we keep AT fixed (or slightly higher) and let N increase.
                at_actual_s = at_final 
                optimised_n = rr_actual * at_actual_s
            else:
                # If Standard (Flooring), we prioritize the exact Pulse Count, 
                # so we adjust AT to match N / RR.
                at_actual_s = n_val / rr_actual
                optimised_n = n_val
        
        # Snap to hardware precision step (e.g. 0.1 ms for Agilent)
        p_s = inputs['precision_ms'] / 1000.0
        # Standard + MC systems: Snap to nearest hardware step
        at_actual_s = round(at_actual_s / p_s) * p_s
        
        return at_actual_s, at_needed, harmonic, valid_times_s, notes, optimised_n

    @staticmethod
    def calculate_laser_sync(inputs):
        """
        Calculates the optimal Acquisition Time (AT), Laser Repetition Rate (RR),
        and Stage Speed based on hardware limits.
        """
        fw_s = inputs['washout_ms'] / 1000.0
        bs = inputs['spot_size_um']
        n_val = inputs['pulses_per_pixel']
        n_target = n_val # Preserve original user request for error checking (Capture BEFORE valid_times/constrained_at)
        
        rr_max = inputs['max_rr_hz']
        ss_max = inputs['max_speed_um_s']
        allowed_rr = inputs.get('allowed_rr', None)
        is_mc = inputs['icp_technology'] == "Multi-Collector"
        
        valid_times_s = [d/1000.0 for d in (inputs.get('allowed_dwells') or [])]
        
        # 1. Get definitive AT from helper (Single Source of Truth)
        # Refined: at_final_s = Actual (post-rounding), at_target_s = Theoretical Needs
        at_final_s, at_target_s, _, _, notes, opt_n = Logic.calculate_constrained_at(inputs)
        
        inputs['pulses_per_pixel'] = opt_n
        
        # 2. Calculate Derived Laser Params based on AT
        # Use OPTIMISED pulses (if modified by MC logic or Avoid Gaps scaling)
        n_val = opt_n 
        rr_final_raw = n_val / at_final_s
        
        # SNAP Rep Rate to Hardware Precision (to reveal actual dosage deviation)
        rr_prec = inputs.get('rr_prec_hz', 1.0)
        if rr_prec <= 0: rr_prec = 1.0
        
        # Determine Rounding Strategy
        avoid_gaps = inputs.get('avoid_gaps', False)
        
        if is_mc:
            # MC Mode: Dwell is fixed. RR = n / AT.
            # Avoid Gaps -> CEIL (Snap UP to ensure overlap)
            # Standard   -> FLOOR (Snap DOWN to ensure budget/pulse constraint)
            if avoid_gaps:
                 rr_final = math.ceil(rr_final_raw / rr_prec - 1e-9) * rr_prec
            else:
                 rr_final = math.floor(rr_final_raw / rr_prec + 1e-9) * rr_prec
        else:
            # Standard Mode: AT is tuned to match RR.
            # calculate_constrained_at already aligned AT to be compatible with a valid RR.
            # We round to nearest to recover that integer.
            rr_final = round(rr_final_raw / rr_prec) * rr_prec
        
        # Update inputs with optimised N so downstream usage (e.g. results table) is correct
        inputs['pulses_per_pixel'] = opt_n
        
        warning = "" 
        
        # 3. Recalculate Actual Pulses based on REALIZABLE Rep Rate and AT
        actual_pulses = rr_final * at_final_s
        
        
        
        
        
        # Check for Asynchronous Mapping (Integrated pixels v Laser shots mismatch)
        # Always report this to the user as per request
        # Compare against ORIGINAL target (n_target), not the optimised one (which moves with logic)
        sync_error = actual_pulses - n_target
        err_pct = (sync_error / n_target) * 100 if n_target != 0 else 0.0
        
        if is_mc:
             # MC Time-Driven Adjustment Message
             warning = f"Optimised to match Integration Time: {actual_pulses:.4f} pulses (Error: {err_pct:+.4f}%)."
        else:
             # Standard system warning
             warning = f"Calculated Pulses per Dwell Time: {actual_pulses:.4f} (Error: {err_pct:+.4f}%)."

        # Diagnostic Metrics
        budget_ms = (at_final_s * 1000) - inputs.get('overhead_ms', 0)
        
        # Mapping-Priority Speed (Sync strictly to laser shots)
        # We derive speed from Spot Size and Rep Rate to guarantee 0.0 overlap for 1-pulse pixels.
        speed = (bs * rr_final) / n_val
        overlap_um = bs - (speed / rr_final)
        overlap_pct = (overlap_um / bs) * 100 if bs > 0 else 0
        
        
        # Diagnostic
        
        prec_ms = inputs['precision_ms']
        
        # Diagnostic
        # IoLog.debug(f"Logic: bs={bs}, rr={rr_final}, n={n_val}, calc_speed={speed}")
        
        result = {
            "Target Acquisition Time (ms)": round(at_target_s * 1000.0, 3), # Target is theoretical
            "Acquisition Time (ms)": round(at_final_s * 1000.0 / prec_ms) * prec_ms,
            "Laser Rep Rate (Hz)": rr_final,
            "Stage Speed (µm s⁻¹)": float(speed),
            "Actual Pulses": actual_pulses,
            "Dwell Budget (ms)": round(budget_ms / prec_ms) * prec_ms,
            "Overhead (ms)": round(inputs.get('overhead_ms', 0) / prec_ms) * prec_ms,
            "Overlap (%)": overlap_pct,
            "Overlap (µm)": overlap_um,
            "Warning": warning,
            "Notes": notes
        }
        return result

    @staticmethod
    def calculate_sigma_for_spot(spot_size, config, isotope_data):
        target_sigma = config['lower_sigma_limit']
        precision_s = config['precision_ms'] / 1000.0
        tech = config['icp_technology']
        system = config['system_type']
        precision_ms = config['precision_ms']
        precision_s = precision_ms / 1000.0
        snr_threshold = config.get('snr_threshold', 0.0)
        
        # Calculate Logic-Derived Budget from Constraints
        c_inputs = config.copy()
        c_inputs['spot_size_um'] = spot_size
        at_constrained, _, _, _, _, _ = Logic.calculate_constrained_at(c_inputs)
        
        dwell_budget_s = at_constrained - (config['overhead_ms'] / 1000.0)
        min_dwell_s = config['min_dwell_ms'] / 1000.0
        
        decimals = 0 if precision_ms >= 1.0 else int(math.ceil(-math.log10(precision_ms)))
        scaling_factor = (spot_size / config['ref_spot_size_um'])**2
        
        rr_scaling = 1.0
        if config.get('scale_signal', False):
            initial_rr = config.get('initial_rr', 0)
            if initial_rr > 0:
                pulses = config.get('pulses_per_pixel', 1) 
                
                # Use the ACTUAL constrained acquisition time to determine the target RR
                # This ensures consistent scaling even if washout logic differs slightly from hardware constraints
                at_s = at_constrained
                if at_s > 0:
                     target_rr = pulses / at_s
                else: 
                     washout_s = config['washout_ms'] / 1000.0
                     ideal_rr = pulses / washout_s
                     max_rr = config.get('max_rr_hz', 1000)
                     target_rr = min(ideal_rr, max_rr)
                     
                rr_scaling = target_rr / initial_rr

        processed_iso = []
        for iso in isotope_data:
            bl_dt = iso['baseline_dt_s']
            sig_cps = iso['sig_cps']
            blk_cps = iso['blk_cps']
            std_counts = iso['stdev_blank_counts']
            
            # 1. Determine the Noise Floor (Count Domain)
            total_bg_counts = blk_cps * bl_dt
            
            # Theoretical Poisson Variance
            theoretical_var = total_bg_counts
            
            # Observed Variance
            observed_var = std_counts**2
            
            # FIX 1: Use MAX, not SUM. Avoid double-counting noise.
            effective_noise_var = max(theoretical_var, observed_var)
            
            # 2. Calculate Critical Level (Lc) for Detection
            d = 0.4
            z = 1.645  # 95% Confidence
            
            if total_bg_counts == 0:
                is_zero_bg = True
                L_c_counts = (z**2)/2 + z * math.sqrt(2 * d) 
            else:
                is_zero_bg = False
                L_c_counts = z * math.sqrt(2 * effective_noise_var)

            # 3. Calculate "Robust SNR" (Detection Metric)
            net_cps = sig_cps - blk_cps
            
            # SCALING: Apply spot size and rep rate scaling to the net counts
            scaled_net_counts_expected = (net_cps * scaling_factor * rr_scaling) * bl_dt
            
            # We calculate internal Detection Ratio (DR = net / Lc)
            # then scale it back to "Sigma-equivalent" units so the UI remains familiar.
            # Factor = z * sqrt(2) ~ 2.326
            sigma_factor = z * math.sqrt(2)
            
            base_snr = 0
            if L_c_counts > 0:
                detection_ratio = scaled_net_counts_expected / L_c_counts
                base_snr = detection_ratio * sigma_factor # Now back in "Sigma SNR" language

            # Scaling logic
            scaled_net_cps = net_cps * scaling_factor * rr_scaling
            sig_cps_new_bs = scaled_net_cps + blk_cps

            initial_snr_display = iso.get('initial_snr', 0.0)
            is_optimisable = (base_snr > 0 and iso['status'] not in ["Exclude", "Set to Min", "Custom"])

            processed_iso.append({
                **iso, 
                'sig_cps_new_bs': sig_cps_new_bs, 
                'std_cps_per_sqrt_sec': math.sqrt(effective_noise_var) / bl_dt if bl_dt > 0 else 0,
                'base_snr': base_snr,  # This is now "Detection Ratio" at bl_dt
                'initial_snr_display': initial_snr_display,
                'is_optimisable': is_optimisable, 
                'final_dt': 0.0, 
                'snapped_dt': 0.0,
                'is_zero_bg': is_zero_bg
            })

        if tech == "Multi-Collector":
            # For MC, we use the allowed_dwells from config (Integration Times)
            mc_times = [d/1000.0 for d in (config.get('allowed_dwells') or [])]
            if not mc_times:
                mc_times = [0.066, 0.131, 0.262, 0.524, 1.049] if system == "Neptune" else [0.050, 0.100, 0.250, 0.500, 1.000]
            
            custom_dwells = [iso.get('custom_time_s', min_dwell_s) for iso in processed_iso if iso['status'] == "Custom"]
            
            if custom_dwells:
                selected_time = max(custom_dwells)
            else:
                potential_times = [t for t in mc_times if t >= (dwell_budget_s - 1e-5)]
                if potential_times:
                    selected_time = min(potential_times)
                else:
                    selected_time = max(mc_times)
                
            for iso in processed_iso:
                if iso['status'] == "Exclude":
                    iso['snapped_dt'] = 0
                else:
                    # Resultant SNR Adjustment: 
                    # base_snr is Detection Ratio at bl_dt. New DR scales with sqrt(dt / bl_dt).
                    bl_dt = iso['baseline_dt_s']
                    resultant_snr = iso['base_snr'] * math.sqrt(selected_time / bl_dt) if bl_dt > 0 else 0
                    if resultant_snr < snr_threshold:
                        iso['snapped_dt'] = selected_time
                        iso['constraint'] = "Min SNR"
                    else:
                        iso['snapped_dt'] = selected_time
        elif tech == "TOF":
            simultaneous_time = math.floor(dwell_budget_s / precision_s) * precision_s
            for iso in processed_iso:
                if iso['status'] == "Exclude":
                    iso['snapped_dt'] = 0
                else:
                    bl_dt = iso['baseline_dt_s']
                    resultant_snr = iso['base_snr'] * math.sqrt(simultaneous_time / bl_dt) if bl_dt > 0 else 0
                    if resultant_snr < snr_threshold:
                        iso['snapped_dt'] = simultaneous_time
                        iso['constraint'] = "Min SNR"
                    else:
                        iso['snapped_dt'] = simultaneous_time
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
                        # We normalize the SNR to 1s for the share calculation.
                        # base_snr is at bl_dt. SNR_at_1s = base_snr / sqrt(bl_dt).
                        # inv_snr = 1 / (base_snr / sqrt(bl_dt)) = sqrt(bl_dt) / base_snr.
                        bl_dt = iso['baseline_dt_s']
                        inv_snr_sum += (math.sqrt(bl_dt) / iso['base_snr']) if iso['base_snr'] > 0 else 0
                
                valid_indices = [i for i, x in enumerate(processed_iso) if x['is_optimisable'] and x['final_dt'] == 0]
                num_valid = len(valid_indices)
                reserved_time = num_valid * min_dwell_s
                extra_budget = current_budget - reserved_time
                
                if extra_budget < 0:
                    for idx in valid_indices: processed_iso[idx]['final_dt'] = min_dwell_s
                elif num_valid > 0 and inv_snr_sum > 0:
                    for idx in valid_indices:
                        bl_dt = processed_iso[idx]['baseline_dt_s']
                        share_ratio = (math.sqrt(bl_dt) / processed_iso[idx]['base_snr']) / inv_snr_sum
                        processed_iso[idx]['final_dt'] = min_dwell_s + (extra_budget * share_ratio)

                failures = []
                for idx in valid_indices:
                    dt = processed_iso[idx]['final_dt']
                    bl_dt = processed_iso[idx]['baseline_dt_s']
                    res_snr = processed_iso[idx]['base_snr'] * math.sqrt(dt / bl_dt) if bl_dt > 0 else 0
                    if res_snr < snr_threshold:
                        failures.append((idx, res_snr))
                
                if not failures: break
                else:
                    failures.sort(key=lambda x: x[1])
                    processed_iso[failures[0][0]]['is_optimisable'] = False

            total_snapped = 0.0
            for iso in processed_iso:
                if iso['status'] == "Exclude": snapped = 0.0
                else:
                    snapped = max(min_dwell_s, math.floor(iso['final_dt'] / precision_s + 1e-9) * precision_s)
                    if not iso['is_optimisable'] and iso['status'] == "Auto": iso['constraint'] = "Min SNR"
                    elif abs(snapped - min_dwell_s) < 1e-9 and iso['status'] == "Auto": iso['constraint'] = "Min ICP"
                    elif iso['status'] == "Set to Min": iso['constraint'] = "Min ICP"
                    else: iso['constraint'] = ""
                iso['snapped_dt'] = snapped; total_snapped += snapped
            
            drift = dwell_budget_s - total_snapped
            if drift >= (precision_s * 0.9):
                candidates = [i for i, x in enumerate(processed_iso) if x['is_optimisable'] and x['status'] not in ["Exclude", "Custom", "Set to Min"]]
                if candidates:
                    candidates.sort(key=lambda i: processed_iso[i]['base_snr'] * math.sqrt(processed_iso[i]['snapped_dt'] / processed_iso[i]['baseline_dt_s']))
                    num_steps = int(round(drift / precision_s))
                    for idx in range(num_steps):
                        processed_iso[candidates[idx % len(candidates)]]['snapped_dt'] += precision_s

        separations = []; output_rows = []
        for iso in processed_iso:
            dt = iso['snapped_dt']
            bl_dt = iso['baseline_dt_s']
            sep = iso['base_snr'] * math.sqrt(dt / bl_dt) if bl_dt > 0 else 0
            if iso['status'] not in ["Exclude", "Set to Min", "Custom"]: separations.append(sep)
            
            val_ms = round(dt * 1000.0, decimals)
            output_rows.append({"Isotope": iso['name'], "Final Dwell (ms)": val_ms if decimals > 0 else int(val_ms),
                                "Initial SNR": round(iso['initial_snr_display'], 2), "Sigma Sep": round(sep, 2),
                                "Status": iso['status'], "Constraint": iso.get('constraint', ""), "IsZeroBG": iso.get('is_zero_bg', False)})

        return (min(separations) if separations else 999.0), output_rows

    @staticmethod
    def calculate_minimum_required_spot_size(config, isotope_data):
        target_sigma = config['lower_sigma_limit']
        low, high = 1, 300
        best_spot, best_res = high, None
        while low <= high:
            mid = (low + high) // 2 
            sigma, rows = Logic.calculate_sigma_for_spot(mid, config, isotope_data)
            if sigma >= target_sigma: best_spot = mid; best_res = rows; high = mid - 1
            else: low = mid + 1
        if best_res is None:
            _, best_res = Logic.calculate_sigma_for_spot(300, config, isotope_data)
        return best_spot, best_res

    @staticmethod
    def calculate_signal_statistics(df, isotope_cols, bg_range, sig_range, edited_dwells_df, estimated_dwell_ms, is_raw_counts, show_counts_check):
        df_blk = df[(df["Time"] >= bg_range[0]) & (df["Time"] <= bg_range[1])]
        df_sig = df[(df["Time"] >= sig_range[0]) & (df["Time"] <= sig_range[1])]
        temp_data = []
        for col in isotope_cols:
            dwell_s = edited_dwells_df[col].iloc[0] / 1000.0 if col in edited_dwells_df.columns else estimated_dwell_ms / 1000.0
            m_sig, m_blk, s_blk = df_sig[col].mean(), df_blk[col].mean(), df_blk[col].std()
            
            # Robust SNR (Detection Ratio) Calculation for Display
            # This must match the Logic.calculate_sigma_for_spot implementation
            sig_cps = m_sig if not is_raw_counts else m_sig / dwell_s
            blk_cps = m_blk if not is_raw_counts else m_blk / dwell_s
            std_counts = s_blk if is_raw_counts else s_blk * dwell_s
            
            bl_dt = dwell_s
            total_bg_counts = blk_cps * bl_dt
            theoretical_var = total_bg_counts
            observed_var = std_counts**2
            effective_noise_var = max(theoretical_var, observed_var)
            
            # Factor = z * sqrt(2) ~ 2.326 to restore familiar Sigma scaling
            sigma_factor = z * math.sqrt(2)
            
            net_counts = (sig_cps - blk_cps) * bl_dt
            detection_ratio = net_counts / L_c_counts if L_c_counts > 0 else 0
            snr = detection_ratio * sigma_factor # Back to "Sigma SNR" for UI consistency

            factor = 1.0 if (is_raw_counts == show_counts_check) else (dwell_s if show_counts_check else 1/dwell_s)
            temp_data.append({"name": col, "current_dwell_ms": dwell_s*1000, "disp_sig": m_sig*factor, "disp_std_sig": df_sig[col].std()*factor, "disp_blk": m_blk*factor, "disp_std_blk": s_blk*factor, "initial_snr": snr})
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
            row['sig_cps'] = df_sig[name].mean() if not is_raw_counts else df_sig[name].mean() / dwell_s
            row['blk_cps'] = df_blk[name].mean() if not is_raw_counts else df_blk[name].mean() / dwell_s
            row['stdev_blank_counts'] = df_blk[name].std() if is_raw_counts else df_blk[name].std() * dwell_s
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
        
        # Add a Close button at the bottom
        btns = QHBoxLayout()
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
        import os
        home = os.path.expanduser("~")
        # Define the base directory (matches main class logic)
        base_dir = os.path.join(home, "Documents", "iolite", "iolite Optimiser")
        if not os.path.exists(base_dir):
            alt_path = os.path.join(home, "OneDrive", "Documents", "iolite", "iolite Optimiser")
            if os.path.exists(alt_path):
                base_dir = alt_path
        
        # Ensure directory exists if we're trying to save
        if not os.path.exists(base_dir):
            try:
                os.makedirs(base_dir)
            except: pass
            
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
        self.spr_debounce_timer.timeout.connect(self.run_spr_analysis)

        # Initialize hardware state
        self.icp_tech = "Quadrupole"
        self.min_dwell = 0.1
        self.precision = 0.1
        self.max_rr = 1000
        self.allowed_rr = None
        self.allowed_dwells = None
        self.max_speed = 20000
        
        try:
            home = os.path.expanduser("~")
            ol_path = os.path.join(home, "Documents", "iolite", "iolite Optimiser")
            if not os.path.exists(ol_path):
                 # Try OneDrive fallback
                 alt_path = os.path.join(home, "OneDrive", "Documents", "iolite", "iolite Optimiser")
                 if os.path.exists(os.path.join(home, "OneDrive", "Documents", "iolite")):
                     if not os.path.exists(alt_path): os.makedirs(alt_path)
                     ol_path = alt_path
                 else:
                     os.makedirs(ol_path)

            self.settings_json_path = os.path.join(ol_path, "iolite Optimiser Settings.json")
            self.persistent_settings = self.load_persistent_settings()
        except Exception as e:
            self.persistent_settings = {}
            IoLog.warning(f"iolite Optimiser: Error creating settings path: {e}")

        # Create Hardware Settings Dialog early
        self.settings_dlg = SettingsDialog(self)
        
        # Pre-hydration removed as it is handled by _handle_mfr_changed in initUI
        
        self.initUI()

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
        
        # Select Optimiser by default for now (or SPR if desired)
        self.tabs.setCurrentIndex(1)
        
        main_layout.addWidget(self.tabs)
        self.setLayout(main_layout)
        
    def init_spr_tab(self):
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()

        # --- LEFT COLUMN (CONTROLS) ---
        
        # 0. Reload Button (Fixed at Top)
        h_ctrl = QHBoxLayout()
        self.btn_spr_run = QPushButton("Reload Data from iolite")
        self.btn_spr_run.setFixedHeight(30)
        self.btn_spr_run.clicked.connect(lambda: self.refresh_data(tab="spr"))
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
        
        self.cmb_spr_iso = QComboBox()
        self.cmb_spr_iso.currentTextChanged.connect(self.run_spr_analysis)
        grid_det.addWidget(QLabel("Select Channel:"), 0, 0)
        grid_det.addWidget(self.cmb_spr_iso, 0, 1)
        
        l_prom = QHBoxLayout()
        l_prom.setContentsMargins(0, 0, 0, 0)
        self.spin_spr_prom = QDoubleSpinBox()
        self.spin_spr_prom.setRange(0, 1e9)
        self.spin_spr_prom.setDecimals(0)
        self.spin_spr_prom.setValue(100.0)
        self.spin_spr_prom.setSingleStep(100.0)
        self.spin_spr_prom.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.spin_spr_prom.valueChanged.connect(self._on_spr_prom_changed)
        l_prom.addWidget(self.spin_spr_prom)
        
        self.chk_spr_auto_prom = QCheckBox("Auto")
        self.chk_spr_auto_prom.setChecked(True)
        self.chk_spr_auto_prom.toggled.connect(self._on_spr_auto_prom_toggled)
        l_prom.addWidget(self.chk_spr_auto_prom)
        # self.spin_spr_prom.setEnabled(False) # Removed: User wants always clickable
        
        grid_det.addWidget(QLabel("Peak Cutoff:"), 1, 0)
        grid_det.addLayout(l_prom, 1, 1)
        
        self.spin_spr_dist = QSpinBox()
        self.spin_spr_dist.setRange(1, 100000)
        self.spin_spr_dist.setValue(self.persistent_settings.get('spr_min_distance', 10))
        self.spin_spr_dist.setSuffix(" data pts.")
        self.spin_spr_dist.setToolTip("Minimum distance between peaks. Prevents detecting multiple points on the same peak.")
        self.spin_spr_dist.valueChanged.connect(self._on_spr_dist_changed)
        self.spin_spr_dist.valueChanged.connect(self.save_persistent_settings)
        grid_det.addWidget(QLabel("Min. Distance:"), 2, 0)
        grid_det.addWidget(self.spin_spr_dist, 2, 1)
        
        self.cmb_spr_unit = QComboBox()
        self.cmb_spr_unit.addItems(["Seconds (s)", "Milliseconds (ms)"])
        # Load selection or default to ms
        self.cmb_spr_unit.setCurrentText(self.persistent_settings.get('spr_time_unit', "Milliseconds (ms)"))
        self.cmb_spr_unit.currentTextChanged.connect(self.run_spr_analysis)
        self.cmb_spr_unit.currentTextChanged.connect(self.save_persistent_settings)
        grid_det.addWidget(QLabel("Time Unit:"), 3, 0)
        grid_det.addWidget(self.cmb_spr_unit, 3, 1)
        
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
        
        l_settings.addStretch()
        self.spr_scroll.setWidget(scroll_content)
        left_layout.addWidget(self.spr_scroll)
        
        # --- RIGHT COLUMN (PLOT & TABLE) ---
        h_ctrl_spr = QHBoxLayout()
        
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
        
        h_ctrl_spr.addWidget(QLabel("Theme:"))
        h_ctrl_spr.addWidget(self.combo_theme_spr)
        h_ctrl_spr.addSpacing(10)
        h_ctrl_spr.addWidget(self.chk_y_zoom_spr)
        h_ctrl_spr.addWidget(self.chk_rescale_spr)
        h_ctrl_spr.addStretch()
        
        # Hide Summary Stats Checkbox (Right Aligned)
        self.chk_hide_stats = QCheckBox("Hide Summary Stats")
        self.chk_hide_stats.setChecked(False)
        self.chk_hide_stats.toggled.connect(self._on_stats_toggled)
        h_ctrl_spr.addWidget(self.chk_hide_stats)
        
        right_layout.addLayout(h_ctrl_spr)

        self.spr_figure = Figure(figsize=(5, 4), dpi=100, constrained_layout=False)
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
        
        self.spr_table = CopyableTableWidget()
        self.spr_table.setColumnCount(7)
        self.spr_table.setHorizontalHeaderLabels(["Peak", "Time (s)", "FW0.1M (10 %) (s)", "FW0.1M (10 %) Area", "FW0.01M (1 %) (s)", "FW0.01M (1 %) Area", "Max Intensity"])
        self.spr_table.verticalHeader().setVisible(False)
        self.spr_table.horizontalHeader().setDefaultAlignment(Qt.AlignCenter)
        self.spr_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        right_layout.addWidget(self.spr_table, 1)
        
        main_layout.addLayout(left_layout)
        main_layout.addLayout(right_layout)
        self.tab_spr.setLayout(main_layout)

        # Exclusion tracking
        self.spr_excluded_peaks = {} # isotope -> set of peak indices
        self.spr_peak_label_map = {} # artist -> peak index

    def run_spr_analysis(self, *args, rescale=None):
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
            
        iso = self._get(self.cmb_spr_iso, 'currentText')
        prominence = self._get(self.spin_spr_prom, 'value')
        unit_text = self._get(self.cmb_spr_unit, 'currentText')
        distance = self._get(self.spin_spr_dist, 'value')
        
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

        # 1. Detect All Peaks (Cache or find new)
        # For simplicity we re-detect if prominence/distance changed, but we keep results_df for toggling
        self.spr_raw_results_df = Logic.analyse_washout_peaks(self.spr_df, iso, prominence, min_distance=distance, is_cps=is_cps)
        
        # 2. Filter Excluded Peaks for Statistics
        if iso not in self.spr_excluded_peaks:
            self.spr_excluded_peaks[iso] = set()
            
        excluded = self.spr_excluded_peaks[iso]
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
        if "Error" in stats:
            self.lbl_spr_fw10.setText("Error")
            self.lbl_spr_fw1.setText(stats["Error"])
            return
            
        if stats.get("Count", 0) == 0:
            self.lbl_spr_fw10.setText("-")
            self.lbl_spr_rsd10.setText("-")
            self.lbl_spr_fw_max10.setText("-")
            self.lbl_spr_fw1.setText("No peaks detected")
            self.lbl_spr_rsd1.setText("-")
            self.lbl_spr_fw_max1.setText("-")
            return

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
        
        # Helper for adaptive rounding (float)
        def _adap_round(v_ms):
            if v_ms >= 100: return round(v_ms, 0)
            if v_ms >= 10:  return round(v_ms, 1)
            return round(v_ms, 2)
        
        # Store for "Apply" logic (rounded to matching UI precision)
        self._last_spr_fw10_avg_ms = _adap_round(stats['FW0.1M Mean'] * 1000.0)
        self._last_spr_fw10_max_ms = _adap_round(stats['FW0.1M Max'] * 1000.0)
        self._last_spr_fw1_avg_ms = _adap_round(stats['FW0.01M Mean'] * 1000.0)
        self._last_spr_fw1_max_ms = _adap_round(stats['FW0.01M Max'] * 1000.0)

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

        # Plot individual channel (Always seconds X)
        # Uses first color from palette (now Purple)
        line, = ax.plot(t_zeroed, self.spr_df[iso], lw=1.5, alpha=0.9, label=iso)
        
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
            comp_df = Logic.generate_composite_peak(self.spr_df, iso, filtered_df)
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
        
        # Anchored Y-Axis Label (Fixed 20px from left edge)
        import matplotlib.transforms as mtransforms
        unit = getattr(self, 'cached_unit_label', "Counts")
        ax.set_ylabel(f"Intensity ({unit})", fontsize='medium', color=fg_col)
        # Use 20px to be safe from bezel
        ax.yaxis.set_label_coords(20/self.spr_figure.dpi/self.spr_figure.get_size_inches()[0], 0.5, 
                                  transform=mtransforms.blended_transform_factory(self.spr_figure.transFigure, ax.transAxes))
        
        # Connect dynamic margin update on zoom/pan
        ax.callbacks.connect('ylim_changed', lambda event: self._update_smart_margins(self.spr_canvas))
        
        # Engineering Notation (k, M, G)
        # Matches main Optimiser plot style
        ax.yaxis.set_major_formatter(EngFormatter(places=0, sep=" "))
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
        # Use legendHandles to ensure we get all markers/lines correctly
        for legline, legtext in zip(leg.legendHandles, leg.get_texts()):
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
            
            # Round to integer for UI readability (User Request)
            auto_val = round(auto_val, 0)
                
            self._block_signals(self.spin_spr_prom, True)
            self.spin_spr_prom.setValue(auto_val)
            self._block_signals(self.spin_spr_prom, False)
        except Exception as e:
            # IoLog.error(f"iolite Optimiser: Error in auto-prominence detection: {e}")
            pass

    def _on_spr_prom_changed(self):
        """
        Called when SPR peak cutoff (prominence) is manually adjusted.
        Auto-deselects the 'Auto' checkbox and applies 300ms debouncing.
        """
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

    def _on_spr_dist_changed(self):
        """
        Called when SPR Min Distance is manually adjusted.
        Applies 500ms debouncing.
        """
        # Start/Restart the 500ms debounce timer
        self.spr_debounce_timer.start(500)

    def _on_spr_auto_prom_toggled(self, checked):
        # Removed: self.spin_spr_prom.setEnabled(not checked)
        if checked:
            self.run_spr_analysis()

    def _on_stats_toggled(self, checked):
        # Toggle visibility without re-running analysis (Preserves Zoom)
        if hasattr(self, 'spr_stats_legend') and self.spr_stats_legend:
            self.spr_stats_legend.set_visible(not checked)
            # Efficient redraw if available
            if hasattr(self.spr_canvas, 'draw_idle'):
                self.spr_canvas.draw_idle()
            else:
                self.spr_canvas.draw()

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
        self.spin_wash.setValue(val_ms)
        self.tabs.setCurrentIndex(1) # Switch back to Optimiser
        IoLog.information(f"iolite Optimiser: Applied {mode} washout value: {val_ms:.2f} ms")

    def init_optimiser_tab(self):
        main_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()
        
        # 0. Data Selection & Settings (Fixed at Top)
        h_ctrl = QHBoxLayout()
        self.btn_run = QPushButton("Reload Data from iolite")
        self.btn_run.setFixedHeight(30)
        self.btn_run.clicked.connect(lambda: self.refresh_data(tab="opt"))
        
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
        l_settings.setSpacing(2)

        # 1. Hardware Configuration Summary
        self.grp_hw_summary = QGroupBox("")
        v_hw_sum = QVBoxLayout()
        v_hw_sum.setSpacing(0)
        v_hw_sum.setContentsMargins(5, 2, 5, 2)
        
        lbl_hw_title = QLabel("Hardware Configuration")
        lbl_hw_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_hw_title.setAlignment(Qt.AlignCenter)
        v_hw_sum.addWidget(lbl_hw_title)
        
        l_hw_sum = QVBoxLayout()
        l_hw_sum.setContentsMargins(5, 0, 5, 0)
        l_hw_sum.setSpacing(0)
        self.lbl_hw_sum = QLabel("Initializing hardware...")
        self.lbl_hw_sum.setWordWrap(False)
        self.lbl_hw_sum.setStyleSheet("")
        l_hw_sum.addWidget(self.lbl_hw_sum)
        v_hw_sum.addLayout(l_hw_sum)
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
        
        self.lbl_cust_prec = QLabel("Dwell Precision (ms):"); self.lbl_cust_prec.setFixedWidth(lbl_w)
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

        self.lbl_cust_rr_prec = QLabel("Rep-Rate Precision (Hz):"); self.lbl_cust_rr_prec.setFixedWidth(lbl_w)
        self.spin_cust_rr_prec = QDoubleSpinBox()
        self.spin_cust_rr_prec.setRange(0.001, 100); self.spin_cust_rr_prec.setFixedWidth(100); self.spin_cust_rr_prec.setDecimals(3)
        self.spin_cust_rr_prec.setValue(self.persistent_settings.get('cust_rr_prec', 1.0))
        self.spin_cust_rr_prec.setToolTip("Rep Rate Precision (Hz)")
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
        self.chk_limit_vis.setChecked(self.persistent_settings.get('limit_vis', True))
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
        v_input.setContentsMargins(5, 2, 5, 2)
        
        lbl_in_title = QLabel("Input Parameters")
        lbl_in_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_in_title.setAlignment(Qt.AlignCenter)
        v_input.addWidget(lbl_in_title)
        
        grid_input = QGridLayout()
        grid_input.setContentsMargins(5, 5, 5, 5)
        grid_input.setSpacing(5)

        lbl_w_in = 140
        spin_w_in = 80

        # Row 0
        lbl_spot = QLabel("Initial Spot Size (µm):"); lbl_spot.setFixedWidth(lbl_w_in)
        grid_input.addWidget(lbl_spot, 0, 0)
        self.spin_spot = QSpinBox()
        self.spin_spot.setRange(1, 300); self.spin_spot.setValue(int(self.persistent_settings.get('spot_size', 50)))
        self.spin_spot.setFixedWidth(spin_w_in)
        grid_input.addWidget(self.spin_spot, 0, 1)

        lbl_wash = QLabel("Washout (ms):"); lbl_wash.setFixedWidth(lbl_w_in)
        grid_input.addWidget(lbl_wash, 0, 2)
        self.spin_wash = QDoubleSpinBox()
        self.spin_wash.setRange(1, 50000); self.spin_wash.setValue(self.persistent_settings.get('washout', 500))
        self.spin_wash.setFixedWidth(spin_w_in)
        grid_input.addWidget(self.spin_wash, 0, 3)

        # Row 1
        lbl_init_rr = QLabel("Initial Rep-Rate (Hz):"); lbl_init_rr.setFixedWidth(lbl_w_in)
        grid_input.addWidget(lbl_init_rr, 1, 0)
        self.spin_init_rr = QDoubleSpinBox()
        self.spin_init_rr.setRange(1, 10000); self.spin_init_rr.setValue(self.persistent_settings.get('init_rr', 20))
        self.spin_init_rr.setFixedWidth(spin_w_in)
        grid_input.addWidget(self.spin_init_rr, 1, 1)

        lbl_mode = QLabel("Mode:"); lbl_mode.setFixedWidth(lbl_w_in)
        grid_input.addWidget(lbl_mode, 1, 2)
        self.cmb_mode = QComboBox()
        self.cmb_mode.addItems(["Spot", "Line", "Imaging"])
        self.cmb_mode.setCurrentText(self.persistent_settings.get('mode', 'Spot'))
        self.cmb_mode.setFixedWidth(spin_w_in)
        grid_input.addWidget(self.cmb_mode, 1, 3)

        v_input.addLayout(grid_input)
        grp_input.setLayout(v_input)
        l_settings.addWidget(grp_input)

        # Optimised Settings / Quality Goals
        grp_qual = QGroupBox("")
        v_qual = QVBoxLayout()
        v_qual.setSpacing(0)
        v_qual.setContentsMargins(5, 2, 5, 2)
        
        lbl_qual_title = QLabel("Optimisation Parameters")
        lbl_qual_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_qual_title.setAlignment(Qt.AlignCenter)
        v_qual.addWidget(lbl_qual_title)
        
        self.grid_qual = QGridLayout()
        self.grid_qual.setContentsMargins(5, 5, 5, 5)
        self.grid_qual.setSpacing(5)

        # Row 0
        lbl_pulses = QLabel("Pulses per Acq. Time:"); lbl_pulses.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(lbl_pulses, 0, 0)
        self.spin_pulses = QDoubleSpinBox()
        self.spin_pulses.setRange(1, 100); self.spin_pulses.setValue(self.persistent_settings.get('pulses', 5)); self.spin_pulses.setDecimals(1)
        self.spin_pulses.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_pulses, 0, 1)

        lbl_sigma = QLabel("Target SNR (Sigma):"); lbl_sigma.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(lbl_sigma, 0, 2)
        self.spin_sigma = QDoubleSpinBox()
        self.spin_sigma.setRange(0.1, 100000); self.spin_sigma.setValue(self.persistent_settings.get('target_sigma', 10)); self.spin_sigma.setDecimals(1)
        self.spin_sigma.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_sigma, 0, 3)

        # Row 1 (Swapped: Duty Cycle first, then Min SNR)
        lbl_duty = QLabel("Min Duty Cycle (%):"); lbl_duty.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(lbl_duty, 1, 0)
        self.spin_duty = QDoubleSpinBox()
        self.spin_duty.setRange(0, 100); self.spin_duty.setValue(self.persistent_settings.get('min_duty', 0)); self.spin_duty.setDecimals(1)
        self.spin_duty.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_duty, 1, 1)

        lbl_snr = QLabel("Minimum SNR (Sigma):"); lbl_snr.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(lbl_snr, 1, 2)
        self.spin_snr = QDoubleSpinBox()
        self.spin_snr.setRange(0, 1000); self.spin_snr.setValue(self.persistent_settings.get('min_snr', 0)); self.spin_snr.setDecimals(1)
        self.spin_snr.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_snr, 1, 3)


        # Row 2: Image Dimensions OR Shots
        self.lbl_height = QLabel("Image Height (µm):"); self.lbl_height.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(self.lbl_height, 2, 0)
        self.spin_height = QDoubleSpinBox()
        self.spin_height.setRange(1, 1000000); self.spin_height.setValue(self.persistent_settings.get('img_height', 1000.0)); self.spin_height.setDecimals(0)
        self.spin_height.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_height, 2, 1)

        self.lbl_width = QLabel("Image Width (µm):"); self.lbl_width.setFixedWidth(lbl_w_in)
        self.grid_qual.addWidget(self.lbl_width, 2, 2)
        self.spin_width = QDoubleSpinBox()
        self.spin_width.setRange(1, 1000000); self.spin_width.setValue(self.persistent_settings.get('img_width', 1000.0)); self.spin_width.setDecimals(0)
        self.spin_width.setFixedWidth(spin_w_in)
        self.grid_qual.addWidget(self.spin_width, 2, 3)
        
        # Shots (Overlapping Height in Layout, visibility toggled)
        self.lbl_shots = QLabel("Number of Shots:"); self.lbl_shots.setFixedWidth(lbl_w_in)
        # self.grid_qual.addWidget(self.lbl_shots, 2, 0) # Managed by mode update
        self.spin_shots = QSpinBox()
        self.spin_shots.setRange(1, 100000000); self.spin_shots.setValue(int(self.persistent_settings.get('shots', 100)))
        self.spin_shots.setFixedWidth(spin_w_in)
        # self.grid_qual.addWidget(self.spin_shots, 2, 1) # Managed by mode update

        # Row 3: Bottom check inside the grid
        self.chk_scale = QCheckBox("Scale Signal with Rep-Rate")
        self.chk_scale.setToolTip("Scale sensitivity based on Rep Rate ratio")
        self.chk_scale.setChecked(self.persistent_settings.get('scale_signal', False))
        self.grid_qual.addWidget(self.chk_scale, 3, 0, 1, 2)
        
        v_qual.addLayout(self.grid_qual)
        grp_qual.setLayout(v_qual)
        
        self.chk_avoid_gaps = QCheckBox("Avoid Gaps")
        self.chk_avoid_gaps.setToolTip("Checked: Always Overlap (Increase Rep-Rate). Unchecked: Allow Gaps.")
        self.chk_avoid_gaps.setChecked(self.persistent_settings.get('avoid_gaps', False)) # Default Unchecked
        self.chk_avoid_gaps.stateChanged.connect(self._on_ui_change)
        self.grid_qual.addWidget(self.chk_avoid_gaps, 3, 2, 1, 2)
        l_settings.addWidget(grp_qual)

        # Optimised Settings (Moved from Results)
        grp_sync = QGroupBox("")
        grp_sync.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        v_sync = QVBoxLayout()
        v_sync.setSpacing(0)
        v_sync.setContentsMargins(5, 2, 5, 2)
        
        lbl_sync_title = QLabel("Optimised Settings")
        lbl_sync_title.setStyleSheet("font-weight: bold; font-size: 10pt; padding: 0px; margin: 0px;")
        lbl_sync_title.setAlignment(Qt.AlignCenter)
        v_sync.addWidget(lbl_sync_title)
        
        l_sync = QHBoxLayout()
        # Column 1: Spot Size (Left)
        # Core Optimised Values (Unified Grid)
        grid_sync = QGridLayout()
        grid_sync.setContentsMargins(5,5,5,5)
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

        # Row 1: (Empty Left), Speed, Overhead
        grid_sync.addWidget(QLabel("Speed:"), 1, 2)
        self.lbl_opt_speed = QLabel("- µm s⁻¹")
        grid_sync.addWidget(self.lbl_opt_speed, 1, 3)

        grid_sync.addWidget(QLabel("Overhead:"), 1, 4)
        self.lbl_opt_overhead = QLabel("- ms")
        grid_sync.addWidget(self.lbl_opt_overhead, 1, 5)

        # Row 2: Est Time HH:MM:SS, Overlap, Budget
        grid_sync.addWidget(QLabel("Est. Time:"), 2, 0)
        self.lbl_est_time_hms = QLabel("-")
        grid_sync.addWidget(self.lbl_est_time_hms, 2, 1)

        grid_sync.addWidget(QLabel("Overlap:"), 2, 2)
        self.lbl_opt_overlap = QLabel("- µm")
        grid_sync.addWidget(self.lbl_opt_overlap, 2, 3)

        grid_sync.addWidget(QLabel("Budget:"), 2, 4)
        self.lbl_opt_budget = QLabel("- ms")
        grid_sync.addWidget(self.lbl_opt_budget, 2, 5)

        # Row 3: Est Time Seconds, Dosage, Duty Cycle
        self.lbl_est_time_sec = QLabel("-")
        grid_sync.addWidget(self.lbl_est_time_sec, 3, 1)

        grid_sync.addWidget(QLabel("Dosage:"), 3, 2)
        self.lbl_opt_pulses = QLabel("- Pulses")
        grid_sync.addWidget(self.lbl_opt_pulses, 3, 3)

        grid_sync.addWidget(QLabel("Duty Cycle:"), 3, 4)
        self.lbl_opt_duty = QLabel("- %")
        grid_sync.addWidget(self.lbl_opt_duty, 3, 5)

        # Set stretch for all value columns
        grid_sync.setColumnStretch(1, 1)
        grid_sync.setColumnStretch(3, 1)
        v_sync.addLayout(grid_sync)
        grp_sync.setLayout(v_sync)
        l_settings.addWidget(grp_sync)
        
        self.lbl_result = QLabel("Ready")
        self.lbl_result.setWordWrap(True)
        l_settings.addWidget(self.lbl_result)
        
        l_settings.addStretch()
        self.settings_scroll_area.setWidget(scroll_content)
        self.settings_scroll_area.setFixedWidth(500)
        self.settings_scroll_area.setWidgetResizable(True)
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
        
        self.figure = Figure(figsize=(5, 3), dpi=100, constrained_layout=False); self.canvas = FigureCanvas(self.figure)
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
        
        h_ctrl = QHBoxLayout()
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
        


        h_ctrl.addWidget(QLabel("Theme:")); h_ctrl.addWidget(self.combo_theme)
        h_ctrl.addSpacing(10)
        h_ctrl.addWidget(self.chk_norm); h_ctrl.addWidget(self.chk_y_zoom)
        h_ctrl.addWidget(self.chk_rescale)
        
        h_ctrl.addStretch()
        
        h_ctrl.addWidget(self.chk_auto)
        h_ctrl.addWidget(QLabel("Background (s):"))
        h_ctrl.addWidget(self.spin_bg_start); h_ctrl.addWidget(QLabel("to")); h_ctrl.addWidget(self.spin_bg_end)
        h_ctrl.addSpacing(10)
        h_ctrl.addWidget(QLabel("Signal (s):"))
        h_ctrl.addWidget(self.spin_sig_start); h_ctrl.addWidget(QLabel("to")); h_ctrl.addWidget(self.spin_sig_end)
        
        plot_layout.addLayout(h_ctrl)
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
        self.spin_height.valueChanged.connect(self._on_ui_change)
        self.spin_width.valueChanged.connect(self._on_ui_change)
        self.spin_wash.valueChanged.connect(self._on_ui_change)
        self.spin_init_rr.valueChanged.connect(self._on_ui_change)
        self.cmb_mode.currentTextChanged.connect(self._on_ui_change)
        self.spin_shots.valueChanged.connect(self._on_ui_change)
        self.chk_scale.toggled.connect(self._on_ui_change)
        self.spin_pulses.valueChanged.connect(self._on_ui_change)
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
        
        IoLog.information(f"iolite Optimiser: Overhead updated to {self.detected_overhead_ms:.3f}ms (AT={self.detected_at_ms:.3f}, Dwell={total_init_dwell_ms:.3f})")
        
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

    def _update_icp_status(self):
        if self.allowed_dwells:
            status = f"Allowed Dwells: {', '.join(map(str, self.allowed_dwells))} ms"
        else:
            status = f"Minimum Dwell Time: {self.min_dwell} ms | Precision: {self.precision} ms"
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
            self.lbl_icp_status.setText(f"Type: {self.icp_tech} | {icp_constraints}")

        if hasattr(self, 'lbl_laser_status'):
            rr_info = f"Allowed RRs: {', '.join(map(str, self.allowed_rr))} Hz" if self.allowed_rr else f"Max RR: {self.max_rr} Hz"
            self.lbl_laser_status.setText(f"{rr_info} | Max Speed: {self.max_speed} µm/s | Prec: {self.laser_rr_prec} Hz")

        # 2. Main Window Summary (HTML formatted)
        if hasattr(self, 'lbl_hw_sum'):
            icp_header = "<b style='font-size: 1.1em;'>ICP-MS</b>"
            icp_details = f"Manufacturer: <b>{icp_mfr}</b><br/>Model: <b>{icp_mod}</b><br/>Type: <b>{self.icp_tech}</b>"
            icp_constraints = f"Allowed Dwell times: {', '.join(map(str, self.allowed_dwells))} ms" if self.allowed_dwells else f"Minimum Dwell Time: {self.min_dwell} ms | Precision: {self.precision} ms"
            icp_full = f"{icp_header}<br/>{icp_details}<br/><span style='font-size: 11px;'><b><i>{icp_constraints}</i></b></span>"
            
            las_header = "<b style='font-size: 1.1em;'>Laser</b>"
            las_details = f"Platform: <b>{las_mfr}</b><br/>Laser Source: <b>{las_mod}</b><br/>Cell Type: <b>{cell}</b>"
            rr_info = f"Allowed Rep-Rates: {', '.join(map(str, self.allowed_rr))} Hz" if self.allowed_rr else f"Maximum Rep-Rate: {self.max_rr} Hz"
            las_constraints = f"{rr_info}<br/>Rep-Rate Precision: {self.laser_rr_prec:g} Hz<br/>Max Stage Speed: {self.max_speed} µm s⁻¹"
            las_full = f"{las_header}<br/>{las_details}<br/><span style='font-size: 11px;'><b><i>{las_constraints}</i></b></span>"
            
            self.lbl_hw_sum.setText(f"{icp_full}<br/><br/>{las_full}")

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
        # REMOVED: Saving on every change is inefficient. 
        # relying on closeEvent to save settings.
        # if hasattr(self, 'save_timer'):
        #    self.save_timer.start()
        # else:
        #    self.save_persistent_settings()
        
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
        if hasattr(self, 'lbl_shots'): self.lbl_shots.setVisible(False); self.grid_qual.removeWidget(self.lbl_shots)
        if hasattr(self, 'spin_shots'): self.spin_shots.setVisible(False); self.grid_qual.removeWidget(self.spin_shots)
        
        # 2. Re-populate based on Mode
        if is_spot:
            # Spot: Shots at (2,0)
            if hasattr(self, 'lbl_shots'): 
                self.grid_qual.addWidget(self.lbl_shots, 2, 0)
                self.lbl_shots.setVisible(True)
            if hasattr(self, 'spin_shots'): 
                self.grid_qual.addWidget(self.spin_shots, 2, 1)
                self.spin_shots.setVisible(True)
                
        elif is_line:
            # Line: Length (Width) at (2,0)
            if hasattr(self, 'lbl_width'):
                self.lbl_width.setText("Line Length (µm):")
                self.grid_qual.addWidget(self.lbl_width, 2, 0)
                self.lbl_width.setVisible(True)
            if hasattr(self, 'spin_width'):
                self.grid_qual.addWidget(self.spin_width, 2, 1)
                self.spin_width.setVisible(True)
                
        else: # Imaging
            # Imaging: Height at (2,0), Width at (2,2)
            if hasattr(self, 'lbl_height'):
                self.grid_qual.addWidget(self.lbl_height, 2, 0)
                self.lbl_height.setVisible(True)
            if hasattr(self, 'spin_height'):
                self.grid_qual.addWidget(self.spin_height, 2, 1)
                self.spin_height.setVisible(True)
                
            if hasattr(self, 'lbl_width'):
                self.lbl_width.setText("Image Width (µm):")
                self.grid_qual.addWidget(self.lbl_width, 2, 2)
                self.lbl_width.setVisible(True) # FIXED: Was missing visible call
            if hasattr(self, 'spin_width'):
                self.grid_qual.addWidget(self.spin_width, 2, 3)
                self.spin_width.setVisible(True)



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
                'shots': self._get(self.spin_shots, 'value'),
                'laser_mfr': self._get(self.cmb_laser_mfr, 'currentText'),
                'laser_model': self._get(self.cmb_laser_mod, 'currentText'),
                'spot_size': self._get(self.spin_spot, 'value'),
                'washout': self._get(self.spin_wash, 'value'),
                'init_rr': self._get(self.spin_init_rr, 'value'),
                'scale_signal': self._get(self.chk_scale, 'isChecked'),
                'pulses': self._get(self.spin_pulses, 'value'),
                'target_sigma': self._get(self.spin_sigma, 'value'),
                'min_snr': self._get(self.spin_snr, 'value'),
                'min_duty': self._get(self.spin_duty, 'value'),
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
                'avoid_gaps': self._get(self.chk_avoid_gaps, 'isChecked'),
                
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
                'opt_y_zoom': self._get(self.chk_y_zoom, 'isChecked'),
                # 'opt_rescale_y': Ephemeral
                'spr_y_zoom': self._get(self.chk_y_zoom_spr, 'isChecked'),
                # 'spr_rescale_y': Ephemeral
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



    def refresh_data(self, tab=None):
        # Do NOT force tab index if None. None means "Refresh All".

        try:
            new_df, new_meta = self.get_input_dataframe()
            if new_df is not None:
                IoLog.information(f"iolite Optimiser: Initial data loaded ({len(new_df)} rows)")
                # Log channel count here for correct order
                cols = [c for c in new_df.columns if c != 'Time']
                IoLog.information(f"iolite Optimiser: Found {len(cols)} channels")
                
                # Reset Auto-Rescale to ON for new data (Ephemeral)
                if hasattr(self, 'chk_rescale'): self.chk_rescale.setChecked(True)
                if hasattr(self, 'chk_rescale_spr'): self.chk_rescale_spr.setChecked(True)
            
            # --- GLOBAL PRE-PROCESSING ---
            local_resolved_dwells = {}
            if new_df is not None:
                # Calculate Acquisition Time EARLY so DwellDialog can show it
                if 'Time' in new_df.columns and len(new_df['Time']) > 1:
                    self.detected_at_ms = new_df['Time'].diff().median() * 1000.0
                    IoLog.information(f"iolite Optimiser: Detected Acquisition Time: {self.detected_at_ms:.3f} ms")
                else:
                    self.detected_at_ms = None
                
                # Resolve Dwell Times & Units Globally (Need a snapshot for later assignment)
                cols_all = [c for c in new_df.columns if c != 'Time']
                if cols_all:
                    local_resolved_dwells = self.resolve_dwell_times(cols_all)
            
            # --- SPR TAB DATA ---
            if tab == "spr" or tab is None:
                self.spr_df = new_df
                if new_df is not None:
                    self.spr_metadata = new_meta # Store specifically for SPR
                    self.spr_dwells = local_resolved_dwells.copy() # Store specifically for SPR
                    
                if self.spr_df is not None and 'Time' in self.spr_df.columns and len(self.spr_df['Time']) > 0:
                     cols = [c for c in self.spr_df.columns if c != 'Time']
                     if hasattr(self, 'cmb_spr_iso'):
                         self._block_signals(self.cmb_spr_iso, True)
                         self.cmb_spr_iso.clear()
                         self.cmb_spr_iso.addItems(cols)
                         self._block_signals(self.cmb_spr_iso, False)
                         
                         self.run_spr_analysis()

            # --- OPTIMISER TAB DATA ---
            if tab == "opt" or tab is None:
                self.opt_df = new_df
                if new_df is not None:
                    self.opt_metadata = new_meta # Store specifically for Optimiser
                    self.channel_metadata = self.opt_metadata # Backwards compat alias
                    
                    self.opt_dwells = local_resolved_dwells.copy() # Store specifically for Optimiser
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
                         
                         # Calculate Actual AT
                         self.detected_at_ms = self.opt_df['Time'].diff().median() * 1000.0
    
                         # Resolve Dwells & Calculate Overhead
                         # Dwells already populated in local block and assigned to self.opt_dwells
                         if len(self.opt_df.columns) > 1:
                             self.lbl_result.setText(f"Loaded {len(self.opt_df.columns)-1} channels from iolite.")
                             # Unit/Dwell resolution moved global
                         
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

    def resolve_dwell_times(self, channel_names):
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
                            
                    except: pass
                
                if not found:
                    missing.append(name)
            
            # 2. If missing, check for preconfigured dwells (from TOF dialog) or show dialog
            
            # Log Unit Detection Summary
            n_cps = sum(1 for v in self.channel_is_cps.values() if v)
            n_counts = len(self.channel_is_cps) - n_cps
            IoLog.information(f"iolite Optimiser: Unit Detection - {n_cps} channels detected as CPS, {n_counts} as Counts")
            
            # Update Cached Unit Label
            self.cached_unit_label = "CPS" if n_cps > 0 else "Counts"
            
            if missing:
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

        # Pixel Margins (Fixed offsets)
        l_px = 85  # Y-axis label + SI numbers (Increased for large counts safe-zone)
        r_px = 20  # Buffer
        
        # Dynamic row-aware top padding
        rows = getattr(self, 'last_opt_rows', 1) if canvas == getattr(self, 'canvas', None) else getattr(self, 'last_spr_rows', 1)
        t_px = 22 + (rows * 20) # snug 42px for 1 row, 62px for 2
        
        b_px = 50  # X-axis label + numbers

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

        l_px = 55 # Default fallback (Reduced to 55 to allow tighter fit for small numbers)
        
        if renderer:
            try:
                # FORCE TICK UPDATE: We must force matplotlib to update the tick strings based on the current formatter
                # get_tightbbox triggers this internal update if it hasn't happened yet.
                ax.yaxis.get_tightbbox(renderer)

                # 1. Measure Y-Axis Tick Labels (Max Visible Width)
                tick_width = 0
                ylim = ax.get_ylim()
                # Sort triggers correct view interval if inverted (though not usual here)
                y_min, y_max = sorted(ylim)
                
                for tick in ax.yaxis.get_major_ticks():
                    # Only measure ticks that are actually visible within limits
                    loc = tick.get_loc()
                    if y_min <= loc <= y_max:
                        label = tick.label1 # Left label
                        if label.get_visible():
                            bbox = label.get_window_extent(renderer)
                            if bbox.width > tick_width:
                                tick_width = bbox.width
                
                # 2. Measure Y-Axis Label
                # It is anchored to figure left, but let's assume we want (LabelWidth + Gap + TickWidth + Gap)
                # But wait, we anchored label to 10px from left.
                # So Left Margin should be: 10px + LabelWidth + Gap + TickWidth + Gap
                ylabel = ax.yaxis.get_label()
                label_width = ylabel.get_window_extent(renderer).width
                
                # Gap settings
                gap = 2 
                
                # Total Required Left Margin
                # 20px (Anchor) + LabelWidth + Gap + TickWidth + Gap
                calc_l_px = 20 + label_width + gap + tick_width + 2
                l_px = max(55, calc_l_px) # Floor at 55 to prevent collapse, but allow closeness
                
            except Exception as e:
                # Fallback if calculation fails
                pass
        
        # RE-ASSERT LABEL ANCHOR (Prevent drift during pan/zoom)
        if ax:
            import matplotlib.transforms as mtransforms
            # 20px from left edge
            ax.yaxis.set_label_coords(20/w, 0.5, 
                                      transform=mtransforms.blended_transform_factory(fig.transFigure, ax.transAxes))

        # Apply Margins (w, h already calculated above)
        
        # Use overrides or defaults
        # We need to preserve current top/bottom/right if not provided?
        # But subplots_adjust takes all.
        # We should store them or recalculate them.
        
        # Recalculate defaults if not passed (similar to _on_plot_resize logic)
        if not override_margins:
            rows = getattr(self, 'last_opt_rows', 1) if canvas == getattr(self, 'canvas', None) else getattr(self, 'last_spr_rows', 1)
            t_px = 22 + (rows * 20)
            b_px = 50
            r_px = 20
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

    def on_zoom(self, event):
        if event.inaxes is None or event.inaxes == getattr(self, 'spr_ax_comp', None): return
        ax = event.inaxes
        canvas = event.canvas
        
        # Scale factor
        if event.button == 'up': scale_factor = 0.8 # Zoom IN
        elif event.button == 'down': scale_factor = 1.25 # Zoom OUT
        else: scale_factor = 1.0
        
        # Always Zoom X
        cur_xlim = ax.get_xlim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0])
        xdata = event.xdata
        new_width = cur_xrange * scale_factor
        relx = (cur_xlim[1] - xdata) / cur_xrange
        new_xlim = [xdata - new_width * (1-relx), xdata + new_width * (relx)]
        ax.set_xlim(new_xlim)
        
        # Determine if Y zoom is enabled for this canvas
        is_spr = (hasattr(self, 'spr_canvas') and canvas == self.spr_canvas)
        chk_y = getattr(self, 'chk_y_zoom_spr', None) if is_spr else getattr(self, 'chk_y_zoom', None)
        y_zoom = self._get(chk_y, 'isChecked') if chk_y else False
        
        if y_zoom:
            cur_ylim = ax.get_ylim()
            cur_yrange = (cur_ylim[1] - cur_ylim[0])
            ydata = event.ydata
            new_height = cur_yrange * scale_factor
            rely = (cur_ylim[1] - ydata) / cur_yrange
            new_ylim = [ydata - new_height * (1-rely), ydata + new_height * (rely)]
            ax.set_ylim(new_ylim)
            
        canvas.draw()

    def on_press(self, event):
        if event.inaxes is None or event.inaxes == getattr(self, 'spr_ax_comp', None): return
        if event.dblclick:
            # RESET VIEW (X and Y)
            self.rescale_to_visible(rescale_x=True, rescale_y=True, ax=event.inaxes, canvas=event.canvas)
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
            
            # 1. Handle Region Edge Picking (Optimiser Plot Only)
            if is_shift and event.canvas == getattr(self, 'canvas', None) and event.inaxes:
                if event.xdata is not None:
                    xlim = event.inaxes.get_xlim()
                    # Permissive threshold: 1% of plot width
                    threshold = (xlim[1] - xlim[0]) * 0.01
                    t_shift = getattr(self, 't_start', 0)
                    
                    found_edge = None
                    if hasattr(self, 'bg_times') and self.bg_times:
                        s_bg, e_bg = self.bg_times[0] - t_shift, self.bg_times[1] - t_shift
                        if abs(event.xdata - s_bg) < threshold: found_edge = 'bg_start'
                        elif abs(event.xdata - e_bg) < threshold: found_edge = 'bg_end'
                        
                    if not found_edge and hasattr(self, 'sig_times') and self.sig_times:
                        s_sig, e_sig = self.sig_times[0] - t_shift, self.sig_times[1] - t_shift
                        if abs(event.xdata - s_sig) < threshold: found_edge = 'sig_start'
                        elif abs(event.xdata - e_sig) < threshold: found_edge = 'sig_end'
                        
                    if found_edge:
                        self.dragged_edge = found_edge
                        if hasattr(self, 'chk_auto'):
                            self._block_signals(self.chk_auto, True)
                            self.chk_auto.setChecked(False)
                            self._block_signals(self.chk_auto, False)
                        return # Edge dragging takes precedence
            
            # 2. Initiate Panning (Standard or SPR) - STORE STARTING STATE
            self.plot_panning = True
            self.press_x_pix = event.x
            self.press_y_pix = event.y
            self.start_xlim = event.inaxes.get_xlim()
            self.start_ylim = event.inaxes.get_ylim()
            self.press_ax = event.inaxes
            # Capture data width for each axis once at start
            self.start_dx_per_pix = (self.start_xlim[1] - self.start_xlim[0]) / (event.inaxes.bbox.width or 1)
            self.start_dy_per_pix = (self.start_ylim[1] - self.start_ylim[0]) / (event.inaxes.bbox.height or 1)

    def on_drag(self, event):
        if event.inaxes is None: return
        
        # 1. Handle Edge Dragging (Priority)
        if self.dragged_edge and event.xdata is not None:
            val = max(0, event.xdata)
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
                
            if sb:
                self._block_signals(sb, True)
                sb.setValue(val)
                self._block_signals(sb, False)
            
            # Capture current limits to preserve zoom during full redraw
            if event.inaxes:
                self._opt_axis_limits = (event.inaxes.get_xlim(), event.inaxes.get_ylim())
            
            # Update plot visuals ONLY (no optimization)
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
            
            # Update limits based on displacement from STARTing limits
            if hasattr(self, 'start_xlim'):
                ax.set_xlim(self.start_xlim[0] - dx_data, self.start_xlim[1] - dx_data)
            
            is_spr = (hasattr(self, 'spr_canvas') and canvas == self.spr_canvas)
            chk_y = getattr(self, 'chk_y_zoom_spr', None) if is_spr else getattr(self, 'chk_y_zoom', None)
            if self._get(chk_y, 'isChecked') and hasattr(self, 'start_ylim'):
                ax.set_ylim(self.start_ylim[0] - dy_data, self.start_ylim[1] - dy_data)

            canvas.draw() # More robust than draw_idle
            return

        # 3. Hover Cursor Feedback
        mod = 0
        try: mod = QApplication.instance().keyboardModifiers()
        except: 
            try: mod = QApplication.keyboardModifiers()
            except: pass
        is_shift = bool(mod & Qt.ShiftModifier) or str(getattr(event, 'key', '')).lower() in ('shift', 's')
                   
        if is_shift and event.canvas == getattr(self, 'canvas', None) and event.inaxes:
            xlim = event.inaxes.get_xlim()
            threshold = (xlim[1] - xlim[0]) * 0.01
            t_shift = getattr(self, 't_start', 0)
            
            near_edge = False
            for times in [self.bg_times, self.sig_times]:
                if times:
                    if abs(event.xdata - (times[0]-t_shift)) < threshold or \
                       abs(event.xdata - (times[1]-t_shift)) < threshold:
                        near_edge = True
                        break
            
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

                # Anchored Y-Axis Label (Fixed 20px from left edge)
                import matplotlib.transforms as mtransforms
                ax.set_ylabel("Norm. Intensity" if normalise else f"Intensity ({unit})")
                ax.yaxis.set_label_coords(20/self.figure.dpi/self.figure.get_size_inches()[0], 0.5, 
                                          transform=mtransforms.blended_transform_factory(self.figure.transFigure, ax.transAxes))
                
                
                # Connect dynamic margin update on zoom/pan
                # MOVED TO AXIS CREATION TO PREVENT CALLBACK ACCUMULATION
                # ax.callbacks.connect('ylim_changed', lambda event: self._update_smart_margins(self.canvas))
                
                # Initial Smart Margin Calc (Post-draw simulation)
                # We defer this slightly or call it if renderer is ready
                # self._update_smart_margins(self.canvas) 

                # Apply SI Formatting to Y-Axis (if not normalised)
                if not normalise:
                    # EngFormatter with places=0 (no decimals) and sep=" " (space before unit)
                    ax.yaxis.set_major_formatter(EngFormatter(places=0, sep=" "))
                
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

    def get_input_dataframe(self):
        try:
            channels = data.timeSeriesList(data.Input)
            # IoLog.information(f"iolite Optimiser: Found {len(channels) if channels else 0} channels") # Silent/Moved
            if not channels: return None, {}
            
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
                
                # Calculate AT early for the dialog
                ref_time = channels[0].time()
                if len(ref_time) > 1:
                    detected_at = (ref_time[1] - ref_time[0]) * 1000.0
                else:
                    detected_at = None
                
                # Unified Dialog: Channel selection + Global dwell
                dlg = DataConfigDialog(ch_names, detected_at=detected_at, is_tof=True, parent=self)
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
            
            ref_ch = channels[0]
            time_data = ref_ch.time()
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
                if len(ch.data()) == len(time_data):
                    data_dict[ch.name] = ch.data()
            
            self.show_meta = {'Element': has_el, 'Mass': has_mass}
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
                'washout_ms': self._get(self.spin_wash, 'value'),
                'spot_size_um': fixed_spot_um if fixed_spot_um else self._get(self.spin_spot, 'value'),
                'pulses_per_pixel': self._get(self.spin_pulses, 'value'),
                'mode': self._get(self.cmb_mode, 'currentText'),
                'img_height': self._get(self.spin_height, 'value') if hasattr(self, 'spin_height') else 1000.0,
                'img_width': self._get(self.spin_width, 'value') if hasattr(self, 'spin_width') else 1000.0,
                'avoid_gaps': self._get(self.chk_avoid_gaps, 'isChecked') if hasattr(self, 'chk_avoid_gaps') else False,
                'lower_sigma_limit': self._get(self.spin_sigma, 'value'),
                'min_duty_cycle': self._get(self.spin_duty, 'value') / 100.0 if self._get(self.spin_duty, 'value') else 0.0,
                'snr_threshold': self._get(self.spin_snr, 'value'),
                'overhead_ms': getattr(self, 'detected_overhead_ms', 0.0),
                'ref_spot_size_um': self._get(self.spin_spot, 'value'),
                'scale_signal': self._get(self.chk_scale, 'isChecked'),
                'initial_rr': self._get(self.spin_init_rr, 'value')
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
            is_mc_sys = TYPE_TO_TECH_MAP.get(self.icp_tech) == "Multi-Collector"
            min_dwell_total = 0.0
            if isotope_data:
                if is_mc_sys:
                    # MC: Channels are simultaneous. The budget needed is the maximum 
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

            # Default to first available channel if not specified
            c['log_prefix'] = "Initial Sync"
            sync = Logic.calculate_laser_sync(c)
            
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
            optimum_spot_raw, df_res_opt = Logic.calculate_minimum_required_spot_size(c, isotope_data)
            optimum_spot = int(round(optimum_spot_raw))
            self.optimum_spotsize = optimum_spot
            self.override_spotsize = fixed_spot_um
            
            best_spot = optimum_spot
            df_res = df_res_opt
            
            # Identify limiting channel for the optimum spot
            limiting_iso = "None"
            if df_res_opt:
                # Filter for candidates that were part of the optimization (Auto)
                opt_rows = [r for r in df_res_opt if r.get('Status') == "Auto"]
                if not opt_rows: opt_rows = df_res_opt # Fallback
                if opt_rows:
                    min_row = min(opt_rows, key=lambda x: x['Sigma Sep'])
                    limiting_iso = min_row['Isotope']
            
            final_notes = []
            if fixed_spot_um is not None:
                # Manual Override
                best_spot = fixed_spot_um
                # Re-calculate table for this specific spot
                _, df_res = Logic.calculate_sigma_for_spot(best_spot, c, isotope_data)
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
            
            # Reset 'pulses_per_pixel' to user input to ensure Error Calculation compares against Target, not Optimized result from previous run.
            c['pulses_per_pixel'] = self._get(self.spin_pulses, 'value')
            c['log_prefix'] = "Optimised Sync"
            sync = Logic.calculate_laser_sync(c)
            
            # Inject Warning into the main Spot Size note (Separate Line, White Bullet)
            warn = sync.get('Warning', "")
            if warn and final_notes:
                # Insert after the Spot Size message (index 0)
                # Just insert as a plain string to match standard formatting (Peer bullet)
                final_notes.insert(1, warn)
            
            # Aggregate notes from Logic
            logic_notes = sync.get("Notes", [])
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
                    final_notes.append(("The following channels had zero counts in the background. The detection limit ($L_c$) has been calculated using the Square Root Transform rule for variance stabilization:", ", ".join(zero_bg_isos), "gray"))
            

            
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
            
            # SNAP TO HARDWARE PRECISION at the very end
            prec_ms = self.precision
            sync['Acquisition Time (ms)'] = round(sync['Acquisition Time (ms)'] / prec_ms) * prec_ms
            sync['Dwell Budget (ms)'] = round(sync['Dwell Budget (ms)'] / prec_ms) * prec_ms
            sync['Overhead (ms)'] = round(sync['Overhead (ms)'] / prec_ms) * prec_ms
            
            
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

            self.lbl_opt_rr.setText(f"<b>{sync['Laser Rep Rate (Hz)']:.{laser_dps}f}</b> Hz")
            self.lbl_opt_speed.setText(f"<b>{sync['Stage Speed (µm s⁻¹)']:.1f}</b> µm s⁻¹")
            
            val_at = sync['Acquisition Time (ms)'] / scaler
            val_budget = sync['Dwell Budget (ms)'] / scaler
            val_overhead = sync['Overhead (ms)'] / scaler

            self.lbl_opt_at.setText(f"<b>{val_at:.{disp_dps}f}</b> {unit}")
            self.lbl_opt_budget.setText(f"<b>{val_budget:.{disp_dps}f}</b> {unit}")
            self.lbl_opt_overhead.setText(f"<b>{val_overhead:.{disp_dps}f}</b> {unit}")

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
                
            raw_at = sync.get('Acquisition Time (ms)', 0)
            
            duty_cycle = 0.0
            if raw_at > 0:
                duty_cycle = (sum_opt_dwells / raw_at) * 100.0
            
            self.lbl_opt_duty.setText(f"<b>{duty_cycle:.1f}</b> %")



            if not self._get(self.chk_avoid_gaps, 'isChecked'):
                 # "Avoid Gaps" OFF -> Floor Strategy.
                 self.lbl_opt_pulses.setText(f"{sync['Actual Pulses']:.4f} Pulses")
            else:
                 # "Avoid Gaps" ON -> Ceil Strategy (Overlap).
                 self.lbl_opt_pulses.setText(f"<b>{sync['Actual Pulses']:.4f}</b> Pulses")

            self.lbl_opt_overlap.setText(f"<b>{sync['Overlap (µm)']:.1f}</b> µm")

            # Estimated Scan Time
            # Mode Logic
            mode = self._get(self.cmb_mode, 'currentText')
            est_sec = 0.0
            
            try:
                if mode == "Spot":
                    # Spot: Time = Shots / RepRate
                    shots = float(self._get(self.spin_shots, 'value'))
                    rr = sync.get('Laser Rep Rate (Hz)', 1.0)
                    if rr > 0:
                        est_sec = shots / rr
                        
                elif mode == "Line":
                    # Line: Time = Length / Speed
                    # Overhead is already included/negligible per user request
                    img_w = float(self._get(self.spin_width, 'value')) # Acts as Length
                    speed = float(sync.get("Stage Speed (µm s⁻¹)", 0.0))
                    
                    if speed > 1e-6:
                         est_sec = img_w / speed
                         
                else:
                    # Imaging (Default)
                    # Lines * (Length / Speed)
                    # Overhead removed per user request
                    img_h = float(self._get(self.spin_height, 'value'))
                    img_w = float(self._get(self.spin_width, 'value'))
                    spot_um = float(c.get('spot_size_um', 0.0))
                    speed_ums = float(sync.get("Stage Speed (µm s⁻¹)", 0.0))
                    
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
                
            except Exception as e:
                self.lbl_est_time_hms.setText("Error")
                self.lbl_est_time_sec.setText(f"({str(e)})")
            
            self.lbl_result.setText(self.final_status_base)
            
            # 6. Display Table
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
                
                combo = self.table.cellWidget(i, self.col_map['Mode'])
                if not combo or rebuild_widgets:
                    combo = QComboBox()
                    # Centre text in combo if possible
                    # Note: Some versions of Qt require different approaches, but standard items are always centered now
                    combo.addItems(["Auto", "Set to Min", "Exclude", "Custom"])
                    combo.setCurrentText(current_status)
                    combo.setProperty("iso_name", str(iso_name))
                    combo.activated.connect(lambda idx, n=iso_name, c=combo: self._on_mode_changed(n, c.itemText(idx)))
                    self.table.setCellWidget(i, self.col_map['Mode'], combo)
                else:
                    # Update existing combo without recreating
                    if self._get(combo, 'currentText') != current_status:
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

# --- UI SETUP ---
widget = None


def create_widget():
    global widget
    IoLog.information("iolite Optimiser: create_widget called")
    
    try:
        if widget is not None:
            # Check if C++ object is deleted
            try:
                if not widget.isVisible():
                    widget.show()
                widget.raise_()
                widget.activateWindow()
                
                # Manual Theme Trigger on Re-Open
                # Uses current dropdown value (restored from settings)
                try:
                    # widget.apply_theme() call removed to prevent loop
                    # We rely on changeEvent/showEvent, or call it safely if we are sure.
                    pass
                except Exception as e:
                    IoLog.warning(f"Failed to refresh theme on show: {e}")
                    
                IoLog.information("iolite Optimiser: Existing widget shown")
                return
            except RuntimeError:
                # Object has been deleted (C++ side), but Python wrapper remains
                IoLog.information("iolite Optimiser: Dead C++ object detected")
                widget = None
    except Exception as e:
        IoLog.error(f"iolite Optimiser: Error checking widget state: {e}")
        widget = None

    # Zombie cleanup removed by user request

    # Use 'None' as the parent for the widget per official examples
    IoLog.information("iolite Optimiser: Creating new widget (Parent: None)")
    
    try:
        widget = ioliteOptimiser()
        # widget.setAttribute(Qt.WA_DeleteOnClose) # Keep widget alive to prevent PythonQt crashes
        widget.setWindowTitle(f"iolite Optimiser - v{VERSION}")
        
        # Connect destroyed signal to cleanup global reference
        # widget.destroyed.connect(cleanup_widget)
        
        # Default size if not previously set
        widget.resize(1100, 700) 
        widget.show()
        
        # Connect signals AFTER the widget is fully initialized and shown
        # Auto-Sync Disabled: Signals cause layout instability in this environment
        # widget.connect_signals()
        try:
            screen = QApplication.primaryScreen()
            if screen:
                rect = screen.availableGeometry()
                widget.resize(int(rect.width() * 0.9), int(rect.height() * 0.9))
            else:
                widget.resize(1728, 972)
        except:
             widget.resize(1728, 972)
        
        widget.show()
        widget.raise_()
        widget.activateWindow()
        IoLog.information("iolite Optimiser: Widget initialized and shown")
    except Exception as e:
        IoLog.error(f"iolite Optimiser: Error creating widget: {str(e)}")
        IoLog.error(traceback.format_exc())

def createUIElements():
    # 'ui' is a global object provided by iolite for UI plugins
    IoLog.information("iolite Optimiser: createUIElements called")
    action = QAction("iolite Optimiser", None)
    action.triggered.connect(create_widget)
    ui.setAction(action)
    ui.setMenuName(['Tools'])
