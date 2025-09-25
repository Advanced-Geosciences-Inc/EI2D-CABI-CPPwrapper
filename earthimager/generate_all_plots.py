import re
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import pandas as pd
from matplotlib.colors import LogNorm
import warnings
warnings.filterwarnings('ignore')

class ERTDataProcessor:
    def __init__(self, out_file_path, output_dir="Generated_Images", target_max_depth=None, bottom_width=0.0):
        self.out_file = Path(out_file_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        # ROI shape/depth controls
        self.target_max_depth = target_max_depth  # e.g., 67.0 (m) to match report
        self.bottom_width = float(bottom_width)   # width at bottom of inverted pyramid (m)
        
        # Color scheme and contour options
        self.available_colormaps = ['cividis', 'viridis', 'inferno', 'coolwarm', 'grayscale']
        self.current_colormap = 'jet'  # default colormap
        self.show_contours = False     # default contour setting
        
        if not self.out_file.exists():
            raise FileNotFoundError(f"Cannot find {self.out_file.resolve()}")
        self.lines = self.out_file.read_text(encoding="utf-8", errors="ignore").splitlines()

        # Containers filled during parsing
        self.x_nodes = None
        self.y_nodes = None
        self.x_cent = None
        self.depth_cent = None
        self.electrode_positions = {}
        self.electrode_spacing = None
        self.measured_data = []
        self.convergence_data = []
        self.predicted_data = []
        self.predicted_data_per_iteration = {}  # New: stores V/I_Calc for each iteration
        self.mesh_info = {}
        self.iteration_numbers = []
        self.resistivity_models_per_iter = []  # list of 2D arrays (ney, nex)

        self._parse_data()

    # ---------------------------- low-level helpers ----------------------------
    def extract_section(self, lines, start_marker, end_marker=None, start_idx=0):
        for i in range(start_idx, len(lines)):
            if start_marker in lines[i]:
                j = i + 1
                if end_marker:
                    while j < len(lines) and end_marker not in lines[j]:
                        j += 1
                else:
                    j = len(lines)
                return i, lines[i+1:j]
        return -1, []

    def parse_numbers_sci(self, raw_lines):
        pat = r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?"
        vals = []
        for ln in raw_lines:
            if ln.strip().startswith(";"):
                continue
            tokens = re.findall(pat, ln.replace(",", " "))
            vals.extend(map(float, tokens))
        return np.array(vals, float)

    def parse_node_coordinates(self, raw_lines):
        coords = []
        for ln in raw_lines[2:]:  # skip 2-row header in node blocks
            parts = ln.strip().split(',')
            if len(parts) == 2:
                try:
                    coords.append(float(parts[1]))
                except ValueError:
                    pass
        return np.array(coords, float)

    # ---------------------------- parsing ----------------------------
    def _parse_mesh_size(self):
        """Parse the MESH SIZE block to get parameter counts for ROI depth."""
        idx, block = self.extract_section(self.lines, ";------ MESH SIZE ------", ";------")
        if idx == -1:
            return
        kv = {}
        for ln in block:
            if "=" in ln:
                k, v = ln.split("=")
                kv[k.strip()] = v.strip()
        def _getint(key, default=None):
            try:
                return int(kv.get(key, default))
            except Exception:
                return default
        self.mesh_info = {
            "n_nodes_x": _getint("Number of nodes in X"),
            "n_nodes_y": _getint("Number of nodes in Y"),
            "n_elem_x": _getint("Number of elements in X"),
            "n_elem_y": _getint("Number of elements in Y"),
            "n_param_x": _getint("Number of parameters in X"),
            "n_param_y": _getint("Number of parameters in Y"),
        }

    def _parse_measurement_data(self):
        # electrode positions
        electrode_section_start = -1
        for i, line in enumerate(self.lines):
            if ";------ ELECTRODE LOCATIONS ------" in line:
                electrode_section_start = i + 2
                break
        if electrode_section_start != -1:
            for i in range(electrode_section_start, len(self.lines)):
                line = self.lines[i].strip()
                if not line or line.startswith(";"):
                    break
                parts = line.split(',')
                if len(parts) >= 3:
                    try:
                        elec_id = int(parts[0].strip())
                        x_pos = float(parts[1].strip())
                        y_pos = float(parts[2].strip())
                        self.electrode_positions[elec_id] = (x_pos, y_pos)
                    except ValueError:
                        continue
        
        # after electrode_positions is filled
        xs_sorted = np.sort([xy[0] for xy in self.electrode_positions.values()])
        a = float(np.median(np.diff(xs_sorted))) if len(xs_sorted) > 1 else 1.0
        self.electrode_spacing = a

        # find start of measured data
        data_start_idx = -1
        for i, line in enumerate(self.lines):
            if "Commands, Raw V/I, GeomFactor, AppRes" in line:
                # Look for the actual data header (DataID   A   B   M   N...)
                for j in range(i + 1, min(i + 5, len(self.lines))):
                    if "DataID" in self.lines[j] and "A" in self.lines[j] and "B" in self.lines[j]:
                        data_start_idx = j + 1
                        break
                if data_start_idx != -1:
                    break
        if data_start_idx == -1:
            return
        for i in range(data_start_idx, len(self.lines)):
            line = self.lines[i].strip()
            if not line or line.startswith(";"):
                break
            parts = line.split(',')
            if len(parts) >= 8:
                try:
                    dct = {
                        'DataID': int(parts[0].strip()),
                        'A': int(parts[1].strip()),
                        'B': int(parts[2].strip()),
                        'M': int(parts[3].strip()),
                        'N': int(parts[4].strip()),
                        'VI_measured': float(parts[5].strip()),
                        'K': float(parts[6].strip()),
                        'App_Res_measured': float(parts[7].strip())
                    }
                    if (dct['A'] in self.electrode_positions and 
                        dct['B'] in self.electrode_positions and
                        dct['M'] in self.electrode_positions and 
                        dct['N'] in self.electrode_positions):
                        # Calculate point of attribution using proper array detection
                        dct = self._detect_array_type_and_calculate_attribution(dct)
                        self.measured_data.append(dct)
                except Exception:
                    continue

    def _detect_array_type_and_calculate_attribution(self, dct):
        """Detect array type and calculate point of attribution (x_pos, pseudodepth)"""
        xa = self.electrode_positions[dct['A']][0]
        xb = self.electrode_positions[dct['B']][0]
        xm = self.electrode_positions[dct['M']][0]
        xn = self.electrode_positions[dct['N']][0]
        
        # Calculate electrode separations
        a_spacing = self.electrode_spacing
        
        # Current electrodes midpoint
        x_ab = 0.5 * (xa + xb)
        # Potential electrodes midpoint  
        x_mn = 0.5 * (xm + xn)
        
        # Array center (point of attribution in X)
        x_center = 0.5 * (x_ab + x_mn)
        dct['x_pos'] = x_center
        
        # Determine array type and calculate pseudodepth
        ab_separation = abs(xb - xa)
        mn_separation = abs(xn - xm)
        ab_mn_separation = abs(x_ab - x_mn)
        
        # Array type detection logic
        if ab_separation == 0 and mn_separation == 0:
            # Pole-pole (both AB and MN are point sources)
            array_type = "pole-pole"
            # For pole-pole, pseudodepth is typically a * separation_factor
            n = ab_mn_separation / a_spacing
            pseudodepth = 0.5 * n * a_spacing
        elif abs(ab_separation - 3 * a_spacing) < 0.1 * a_spacing and abs(mn_separation - a_spacing) < 0.1 * a_spacing:
            # Wenner array: AB = 3a, MN = a, where a is electrode spacing
            array_type = "wenner"
            pseudodepth = 0.5 * a_spacing  # Standard Wenner pseudodepth
        elif mn_separation < 0.5 * ab_separation:
            # Schlumberger-like: small MN compared to AB
            array_type = "schlumberger"
            pseudodepth = 0.125 * (ab_separation + mn_separation)
        else:
            # Dipole-dipole or general case
            array_type = "dipole-dipole"
            # For dipole-dipole: pseudodepth ≈ 0.5 * n * a, where n is the separation factor
            n = ab_mn_separation / a_spacing
            pseudodepth = 0.5 * n * a_spacing  # or 0.513 * n * a for some conventions
        
        dct['pseudodepth'] = pseudodepth
        dct['array_type'] = array_type
        
        return dct

    def _parse_convergence_data(self):
        conv_start_idx = -1
        for i, line in enumerate(self.lines):
            if "Iteration  Res-RMS(%)    L2-Norm" in line:
                conv_start_idx = i + 1
                break
        if conv_start_idx == -1:
            return
        for i in range(conv_start_idx, len(self.lines)):
            line = self.lines[i].strip()
            if not line or line.startswith(";"):
                break
            parts = line.split(',')
            if len(parts) >= 3:
                try:
                    self.convergence_data.append({
                        'iteration': int(parts[0].strip()),
                        'rms_error': float(parts[1].strip()),
                        'l2_norm': float(parts[2].strip())
                    })
                except Exception:
                    pass

    def _parse_predicted_data(self):
        """Parse all iterations of V/I_Calc data from OUTPUT OF DATA AND MODEL sections"""
        self.predicted_data_per_iteration = {}
        self.predicted_data = []  # Keep for backward compatibility (will use latest iteration)
        
        # Find the "OUTPUT OF DATA AND MODEL OF ALL ITERATIONS" section
        output_section_start = -1
        for i, line in enumerate(self.lines):
            if ";------ OUTPUT OF DATA AND MODEL OF ALL ITERATIONS ------" in line:
                output_section_start = i
                break
        
        if output_section_start == -1:
            print("Warning: No 'OUTPUT OF DATA AND MODEL OF ALL ITERATIONS' section found")
            return
        
        # Parse all iterations starting from this section
        current_iteration = None
        i = output_section_start + 1
        
        while i < len(self.lines):
            line = self.lines[i]
            
            # Check for iteration marker
            if ";------ Iteration" in line:
                match = re.search(r';------ Iteration\s+(\d+)', line)
                if match:
                    current_iteration = int(match.group(1))
                    self.predicted_data_per_iteration[current_iteration] = []
                i += 1
                continue
            
            # Check for data header
            if current_iteration is not None and ";-Index  V/I_Meas      V/I_Calc    VI_%ERR" in line:
                # Parse data for this iteration
                i += 1
                while i < len(self.lines):
                    data_line = self.lines[i].strip()
                    if not data_line or data_line.startswith(";"):
                        # End of this iteration's data
                        break
                    
                    parts = data_line.split(',')
                    if len(parts) >= 4:
                        try:
                            data_point = {
                                'index': int(parts[0].strip()),
                                'vi_measured': float(parts[1].strip()),
                                'vi_calculated': float(parts[2].strip()),
                                'vi_error': float(parts[3].strip())
                            }
                            self.predicted_data_per_iteration[current_iteration].append(data_point)
                        except Exception as e:
                            print(f"Warning: Could not parse data line: {data_line} - {e}")
                    i += 1
                continue
            
            # Check if we've moved past the output sections (next major section)
            if line.startswith(";------") and "Iteration" not in line and current_iteration is not None:
                # We've likely moved to the next major section
                break
                
            i += 1
        
        # Set predicted_data to the latest iteration for backward compatibility
        if self.predicted_data_per_iteration:
            latest_iteration = max(self.predicted_data_per_iteration.keys())
            self.predicted_data = self.predicted_data_per_iteration[latest_iteration]
            print(f"Parsed predicted data for iterations: {sorted(self.predicted_data_per_iteration.keys())}")

    def _parse_all_resistivity_models(self):
        # node-based centroids must already be available
        nex = len(self.x_cent)
        ney = len(self.depth_cent)
        # find all iteration block indices
        iter_indices = [i for i, ln in enumerate(self.lines) if ";------ Iteration" in ln]
        self.iteration_numbers = []
        self.resistivity_models_per_iter = []
        for idx in iter_indices:
            # record iteration number if present
            m = re.search(r";------ Iteration\s+(\d+)", self.lines[idx])
            itno = int(m.group(1)) if m else len(self.iteration_numbers)
            # extract resistivity block for this iteration
            _, res_block = self.extract_section(
                self.lines[idx:],
                ";-Resistivity in Ohm-m in the elemental sequential order",
                ";-Sensitivity"
            )
            if not res_block:
                continue
            vals = self.parse_numbers_sci(res_block)
            vals = vals[:nex * ney]
            Z = vals.reshape((ney, nex), order="F")
            self.iteration_numbers.append(itno)
            self.resistivity_models_per_iter.append(Z)
        # fall-back: if none found, try last iteration logic
        if not self.resistivity_models_per_iter:
            last_iter_idx = max(i for i, ln in enumerate(self.lines) if ";------ Iteration" in ln)
            _, res_block = self.extract_section(
                self.lines[last_iter_idx:],
                ";-Resistivity in Ohm-m in the elemental sequential order",
                ";-Sensitivity"
            )
            vals = self.parse_numbers_sci(res_block)
            vals = vals[:nex * ney]
            self.iteration_numbers = [0]
            self.resistivity_models_per_iter = [vals.reshape((ney, nex), order="F")]

    def _parse_data(self):
        # nodes
        _, x_block = self.extract_section(self.lines, ";------ NODE LOCATION IN X (m) ------", ";------")
        _, y_block = self.extract_section(self.lines, ";------ NODE LOCATION IN Y (m) ------", ";------")
        self.x_nodes = self.parse_node_coordinates(x_block)
        self.y_nodes = self.parse_node_coordinates(y_block)
        # centroids (store depth positive downward)
        self.x_cent = 0.5 * (self.x_nodes[:-1] + self.x_nodes[1:])
        self.depth_cent = -0.5 * (self.y_nodes[:-1] + self.y_nodes[1:])
        # other blocks
        self._parse_mesh_size()
        self._parse_measurement_data()
        self._parse_convergence_data()
        self._parse_all_resistivity_models()
        self._parse_predicted_data()

    # ---------------------------- ROI helpers ----------------------------
    def _roi_limits(self):
        """Compute ROI limits from electrode spread (X) and parameter depth (Y)."""
        if self.electrode_positions:
            xs = [xy[0] for xy in self.electrode_positions.values()]
            x_min, x_max = min(xs), max(xs)
        else:
            x_min, x_max = float(np.nanmin(self.x_cent)), float(np.nanmax(self.x_cent))
        if self.mesh_info.get("n_param_y"):
            ny = int(self.mesh_info["n_param_y"])  # layers
            ny = min(ny, len(self.depth_cent))
            d_max = float(self.depth_cent[ny-1]) if ny > 0 else float(np.nanmax(self.depth_cent))
        else:
            d_max = float(np.nanmax(self.depth_cent))
        if self.target_max_depth is not None:
            d_max = min(d_max, float(self.target_max_depth))
        d_min = 0.0
        return x_min, x_max, d_min, d_max

    def _roi_masks(self, x_min, x_max, d_min, d_max):
        xmask = (self.x_cent >= x_min) & (self.x_cent <= x_max)
        dmask = (self.depth_cent >= d_min) & (self.depth_cent <= d_max)
        return xmask, dmask

    def _pyramid_mask(self, X, Y, x_min, x_max, d_min, d_max, slope_factor=0.5):
        """Boolean mask with inclined plane boundaries for pseudosections.
        Creates a natural wedge shape that narrows with depth using a configurable slope.
        
        Parameters:
        -----------
        slope_factor : float, default=0.5
            Controls the width at maximum depth as a fraction of surface width.
            0.3 = narrow base (30% of surface width)
            0.5 = moderate base (50% of surface width) 
            0.7 = wide base (70% of surface width)
        """
        cx = 0.5 * (x_min + x_max)
        full_width = x_max - x_min
        
        # Create a simple inclined boundary that narrows naturally with depth
        # The slope factor determines how much narrower it gets with depth
        
        # Calculate the width at each depth level
        depth_ratio = np.clip((Y - d_min) / max(d_max - d_min, 1e-9), 0, 1)
        width_at_depth = full_width * (1.0 - depth_ratio * (1.0 - slope_factor))
        
        # Create symmetric boundaries around center
        half_width = width_at_depth * 0.5
        x_left = cx - half_width
        x_right = cx + half_width
        
        return (X >= x_left) & (X <= x_right)
    
    def set_colormap(self, colormap):
        """Set the colormap for sections and pseudosections.
        
        Parameters:
        -----------
        colormap : str
            One of: 'cividis', 'viridis', 'inferno', 'coolwarm', 'grayscale', 'jet'
        """
        if colormap == 'grayscale':
            self.current_colormap = 'gray'
        elif colormap in self.available_colormaps or colormap == 'jet':
            self.current_colormap = colormap
        else:
            print(f"Warning: Unknown colormap '{colormap}'. Using default 'jet'.")
            print(f"Available colormaps: {self.available_colormaps + ['jet']}")
            self.current_colormap = 'jet'
    
    def set_contours(self, show_contours):
        """Set whether to show contour lines on sections and pseudosections.
        
        Parameters:
        -----------
        show_contours : bool
            If True, shows contour lines. If False, hides them.
        """
        self.show_contours = bool(show_contours)

    # ---------------------------- plotting / exports ----------------------------
    def _export_centroid_table(self, iter_no, Xc_roi, Dc_roi, Z_roi):
        """Save CSV and JSON with columns (x, depth, resistivity)."""
        outstem = f"iteration_{iter_no:02d}"
        # Flatten to columns
        df = pd.DataFrame({
            'x': Xc_roi.ravel(order='C'),
            'depth': Dc_roi.ravel(order='C'),
            'resistivity': Z_roi.ravel(order='C')
        })
        # CSV
        csv_path = self.output_dir / f"{outstem}_roi_centroids.csv"
        df.to_csv(csv_path, index=False)
        # JSON
        json_path = self.output_dir / f"{outstem}_roi_centroids.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(df.to_dict(orient='records'), f, ensure_ascii=False)
        print(f"Saved: {csv_path}")
        print(f"Saved: {json_path}")

    def plot_and_export_inverted_sections_all_iterations(self, slope_factor=0.5):
        """ROI clamp -> interpolate -> plot per iteration, plus CSV/JSON centroid export."""
        if not self.resistivity_models_per_iter:
            print("No resistivity models parsed")
            return
        x_min, x_max, d_min, d_max = self._roi_limits()
        xmask, dmask = self._roi_masks(x_min, x_max, d_min, d_max)
        # grids for interpolation (smooth display)
        xi = np.linspace(x_min, x_max, 400)
        di = np.linspace(d_min, d_max, 200)
        Xi, Di = np.meshgrid(xi, di)
        # mesh of ROI centroids (for export)
        x_roi = self.x_cent[xmask]
        d_roi = self.depth_cent[dmask]
        Xc_roi, Dc_roi = np.meshgrid(x_roi, d_roi)

        for itno, Z in zip(self.iteration_numbers, self.resistivity_models_per_iter):
            Z_roi = Z[dmask, :][:, xmask]
            # interpolate to regular grid
            pts = np.column_stack([Xc_roi.ravel(), Dc_roi.ravel()])
            vals = Z_roi.ravel()
            Zi = griddata(pts, vals, (Xi, Di), method="cubic")
            if np.isnan(Zi).any():
                Zi2 = griddata(pts, vals, (Xi, Di), method="linear")
                Zi[np.isnan(Zi)] = Zi2[np.isnan(Zi)]
            if np.isnan(Zi).any():
                Zi3 = griddata(pts, vals, (Xi, Di), method="nearest")
                Zi[np.isnan(Zi)] = Zi3[np.isnan(Zi)]
            Zi = gaussian_filter(Zi, sigma=0.8, mode="nearest")
            # apply inverted pyramid mask to the interpolated grid
            pyr_mask = self._pyramid_mask(Xi, Di, x_min, x_max, d_min, d_max, slope_factor)
            Zi = np.ma.array(Zi, mask=~pyr_mask)

            # plot
            fig, ax = plt.subplots(figsize=(12, 3.8))
            vmin, vmax = np.nanpercentile(Zi.compressed(), [2, 98])
            cf = ax.contourf(Xi, Di, Zi, levels=50, cmap=self.current_colormap, vmin=vmin, vmax=vmax, antialiased=True)
            
            # Add contour lines if enabled
            if self.show_contours:
                cs = ax.contour(Xi, Di, Zi, levels=10, colors='black', alpha=0.4, linewidths=0.5)
            
            ax.set_aspect('equal', adjustable='box')
            ax.invert_yaxis()
            
            # Set x-axis at top and bottom
            ax.set_xlabel('Distance (m)')
            ax.xaxis.set_label_position('top')
            ax.tick_params(top=True, labeltop=True, bottom=True, labelbottom=False)
            
            ax.set_ylabel('Depth (m)')
            # Title & stats styled like the reference figure
            rms = next((d['rms_error'] for d in self.convergence_data if d['iteration']==itno), None)
            l2 = next((d['l2_norm'] for d in self.convergence_data if d['iteration']==itno), None)
            elec_spacing_txt = ''
            if self.electrode_positions:
                xs_sorted = np.sort([xy[0] for xy in self.electrode_positions.values()])
                if len(xs_sorted) >= 2:
                    dx_med = float(np.median(np.diff(xs_sorted)))
                    elec_spacing_txt = f"    Electrode Spacing = {dx_med:.0f} m"
            stats = []
            if rms is not None: stats.append(f"RMS = {rms:.2f}%")
            if l2 is not None: 
                # Display normalized L2 (average squared residual per data point)
                if len(self.measured_data) > 0:
                    l2_normalized = l2 / len(self.measured_data)
                    stats.append(f"L2 = {l2_normalized:.2f}")
                else:
                    stats.append(f"L2 = {l2:.2f}")
            subtitle = ("    ".join([f"Iteration = {itno}"] + stats) + elec_spacing_txt)
            ax.set_title(f"Inverted Resistivity Section {subtitle}", pad=12)

            # Electrode markers along top
            if self.electrode_positions:
                xs = [xy[0] for xy in self.electrode_positions.values()]
                ax.scatter(xs, np.full_like(xs, d_min), s=16, c='k', marker='s', zorder=3)

            plt.colorbar(cf, ax=ax, label='Resistivity (Ohm·m)')
            plt.tight_layout()
            img_path = self.output_dir / f"Inverted resistivity section — iteration {itno:02d}.png"
            plt.savefig(img_path, dpi=200, bbox_inches='tight')
            print(f"Saved: {img_path}")
            plt.close()

            # exports (centroid table within ROI)
            self._export_centroid_table(itno, Xc_roi, Dc_roi, Z_roi)

    # ---------------------------- existing pseudosections & plots ----------------------------
    def create_pseudosection_grid(self, data_points, value_key, title, filename, slope_factor=0.5):
        if not data_points:
            print(f"No data available for {title}")
            return
        x_coords = np.array([d['x_pos'] for d in data_points])
        y_coords = np.array([d['pseudodepth'] for d in data_points])
        values = np.array([d[value_key] for d in data_points])
        valid_mask = ~np.isnan(values) & (values > 0)
        x_coords, y_coords, values = x_coords[valid_mask], y_coords[valid_mask], values[valid_mask]
        if len(values) == 0:
            print(f"No valid data for {title}")
            return
            
        # Use same ROI limits as inversion sections for consistent depth cutting
        x_min, x_max, d_min, d_max = self._roi_limits()
        
        # Limit pseudodepth data to the same depth as inversion sections
        depth_mask = y_coords <= d_max
        x_coords, y_coords, values = x_coords[depth_mask], y_coords[depth_mask], values[depth_mask]
        
        if len(values) == 0:
            print(f"No data within ROI depth limits for {title}")
            return
        
        # Create grid using ROI limits (same as inversion sections)
        xi = np.linspace(x_min, x_max, 200)
        yi = np.linspace(d_min, d_max, 100)
        Xi, Yi = np.meshgrid(xi, yi)
        try:
            Zi = griddata((x_coords, y_coords), values, (Xi, Yi), method='cubic')
            if np.isnan(Zi).any():
                Zi_linear = griddata((x_coords, y_coords), values, (Xi, Yi), method='linear')
                Zi[np.isnan(Zi)] = Zi_linear[np.isnan(Zi)]
            if np.isnan(Zi).any():
                Zi_nearest = griddata((x_coords, y_coords), values, (Xi, Yi), method='nearest')
                Zi[np.isnan(Zi)] = Zi_nearest[np.isnan(Zi)]
        except Exception as e:
            print(f"Interpolation failed for {title}: {e}")
            return
        
        # Apply gaussian smoothing
        Zi = gaussian_filter(Zi, sigma=0.6, mode="nearest")

        # Apply wedge ROI mask (using ROI limits already obtained above)
        mask = self._pyramid_mask(Xi, Yi, x_min, x_max, d_min, d_max, slope_factor)
        Zi = np.ma.array(Zi, mask=~mask)

        # choose figure size from data extent (no hard-coded height)
        xrange = x_max - x_min
        drange = d_max - d_min
        width = 12.0
        # scale factor keeps the section pleasantly flat but proportional across files
        height = max(2.2, min(4.0, width * (drange / max(xrange, 1e-9)) * 0.55))

        fig, ax = plt.subplots(figsize=(width, height))
        
        # Linear color limits (using percentiles)
        vmin, vmax = np.percentile(values, [2, 98])
        cf = ax.contourf(Xi, Yi, Zi, levels=50, cmap=self.current_colormap, vmin=vmin, vmax=vmax)
        
        # Add contour lines if enabled
        if self.show_contours:
            cs = ax.contour(Xi, Yi, Zi, levels=10, colors='black', alpha=0.4, linewidths=0.5)

        # Electrode markers at the top
        if self.electrode_positions:
            xs = [xy[0] for xy in self.electrode_positions.values()]
            ax.scatter(xs, np.zeros_like(xs), s=16, c="k", marker="s", zorder=3)
        
        ax.invert_yaxis()
        
        # Set x-axis at top 
        ax.set_xlabel('Distance (m)')
        ax.xaxis.set_label_position('top')
        ax.tick_params(top=True, labeltop=True, bottom=True, labelbottom=False)
        
        ax.set_ylabel('Pseudodepth (m)')
        ax.set_title(title)
        
        # don't force 1:1 aspect
        ax.set_aspect('auto')
        plt.colorbar(cf, ax=ax, label='Resistivity (Ohm·m)')
        plt.tight_layout()
        save_path = self.output_dir / filename
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close()

    def create_misfit_pseudosection(self, data_points, value_key, title, filename, slope_factor=0.5):
        if not data_points:
            print(f"No data available for {title}")
            return
        x_coords = np.array([d['x_pos'] for d in data_points])
        y_coords = np.array([d['pseudodepth'] for d in data_points])
        values = np.array([d[value_key] for d in data_points])
        valid_mask = ~np.isnan(values) & (values >= 0)
        x_coords, y_coords, values = x_coords[valid_mask], y_coords[valid_mask], values[valid_mask]
        if len(values) == 0:
            print(f"No valid data for {title}")
            return
            
        # Use same ROI limits as inversion sections for consistent depth cutting
        x_min, x_max, d_min, d_max = self._roi_limits()
        
        # Limit pseudodepth data to the same depth as inversion sections
        depth_mask = y_coords <= d_max
        x_coords, y_coords, values = x_coords[depth_mask], y_coords[depth_mask], values[depth_mask]
        
        if len(values) == 0:
            print(f"No data within ROI depth limits for {title}")
            return
        
        # Create grid using ROI limits (same as inversion sections)
        xi = np.linspace(x_min, x_max, 200)
        yi = np.linspace(d_min, d_max, 100)
        Xi, Yi = np.meshgrid(xi, yi)
        try:
            Zi = griddata((x_coords, y_coords), values, (Xi, Yi), method='cubic')
            if np.isnan(Zi).any():
                Zi_linear = griddata((x_coords, y_coords), values, (Xi, Yi), method='linear')
                Zi[np.isnan(Zi)] = Zi_linear[np.isnan(Zi)]
            if np.isnan(Zi).any():
                Zi_nearest = griddata((x_coords, y_coords), values, (Xi, Yi), method='nearest')
                Zi[np.isnan(Zi)] = Zi_nearest[np.isnan(Zi)]
        except Exception as e:
            print(f"Interpolation failed for {title}: {e}")
            return
        fig, ax = plt.subplots(figsize=(12, 6))
        vmax = np.percentile(values, 95)
        cf = ax.contourf(Xi, Yi, Zi, levels=50, cmap='Reds', vmin=0, vmax=vmax)
        ax.scatter(x_coords, y_coords, c=values, s=20, cmap='Reds', vmin=0, vmax=vmax, edgecolors='black', linewidth=0.5)
        
        ax.invert_yaxis()
        
        # Set x-axis at top
        ax.set_xlabel('Distance (m)')
        ax.xaxis.set_label_position('top')
        ax.tick_params(top=True, labeltop=True, bottom=True, labelbottom=False)
        
        ax.set_ylabel('Pseudodepth (m)')
        ax.set_title(title)
        ax.set_aspect('equal', adjustable='box')
        plt.colorbar(cf, ax=ax, label='Relative Misfit (%)')
        plt.tight_layout()
        save_path = self.output_dir / filename
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close()

    # ---------------------------- plot collection ----------------------------
    def plot_measured_apparent_resistivity(self, slope_factor=0.3):
        """Plot measured apparent resistivity pseudosection (single plot since measured data doesn't change)"""
        self.create_pseudosection_grid(self.measured_data, 'App_Res_measured', 'Measured Apparent Resistivity Pseudosection', 'Measured apparent resistivity pseudosection.jpg', slope_factor)

    def plot_measured_apparent_resistivity_all_iterations(self, slope_factor=0.3):
        """Plot measured apparent resistivity pseudosection for all iterations (for consistency with calculated iterations)"""
        if not self.predicted_data_per_iteration:
            print("No iteration information available")
            return
        
        iterations = sorted(self.predicted_data_per_iteration.keys())
        print(f"Generating measured apparent resistivity pseudosections for iterations: {iterations}")
        
        for iteration in iterations:
            filename = f'Measured apparent resistivity pseudosection — iteration {iteration:02d}.jpg'
            title = f'Measured Apparent Resistivity Pseudosection — Iteration {iteration:02d}'
            self.create_pseudosection_grid(self.measured_data, 'App_Res_measured', title, filename, slope_factor)

    def plot_calculated_apparent_resistivity(self, iteration=None, slope_factor=0.3):
        """Plot calculated apparent resistivity pseudosection.
        If iteration is None, uses the latest iteration.
        """
        if iteration is None:
            # Use latest iteration
            if not self.predicted_data_per_iteration:
                print("No predicted data available for any iteration")
                return
            iteration = max(self.predicted_data_per_iteration.keys())
        
        if iteration not in self.predicted_data_per_iteration:
            print(f"No predicted data available for iteration {iteration}")
            return
        
        predicted_data = self.predicted_data_per_iteration[iteration]
        
        calc_data = []
        # Create a lookup for predicted data by index for faster matching
        pred_lookup = {pred['index']: pred for pred in predicted_data}
        
        for measured in self.measured_data:
            data_id = measured['DataID']
            if data_id in pred_lookup:
                pred = pred_lookup[data_id]
                calc_point = measured.copy()
                # Calculate Rho_app,Calc = K * V/I_Calc
                calc_point['App_Res_calculated'] = measured['K'] * pred['vi_calculated']
                calc_data.append(calc_point)
            else:
                print(f"Warning: No predicted data found for DataID {data_id}")
        
        if not calc_data:
            print("No matching calculated data available")
            return
            
        filename = 'Calculated apparent resistivity pseudosection.jpg'
        title = 'Calculated Apparent Resistivity Pseudosection'
        
        self.create_pseudosection_grid(calc_data, 'App_Res_calculated', title, filename, slope_factor)

    def plot_convergence_curve(self):
        if not self.convergence_data:
            print("No convergence data available")
            return
        iterations = [d['iteration'] for d in self.convergence_data]
        rms = [d['rms_error'] for d in self.convergence_data]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(iterations, rms, 'bo-', linewidth=2, markersize=8)
        ax.set_xlabel('Iteration'); ax.set_ylabel('RMS Error (%)')
        ax.set_title('Convergence Curve of Resistivity Inversion'); ax.grid(True, alpha=0.3); ax.set_ylim(bottom=0)
        for it, r in zip(iterations, rms):
            ax.annotate(f'{r:.2f}%', (it, r), textcoords="offset points", xytext=(0,10), ha='center')
        plt.tight_layout()
        save_path = self.output_dir / 'Convergence curve of resistivity inversion.jpg'
        plt.savefig(save_path, dpi=200, bbox_inches='tight'); print(f"Saved: {save_path}")
        plt.close()

    def plot_all_calculated_apparent_resistivity_iterations(self, slope_factor=0.3):
        """Plot calculated apparent resistivity pseudosection for all iterations"""
        if not self.predicted_data_per_iteration:
            print("No predicted data available for any iteration")
            return
        
        # Generate pseudosections for all iterations
        iterations = sorted(self.predicted_data_per_iteration.keys())
        print(f"Generating calculated apparent resistivity pseudosections for iterations: {iterations}")
        
        for iteration in iterations:
            predicted_data = self.predicted_data_per_iteration[iteration]
            
            calc_data = []
            # Create a lookup for predicted data by index for faster matching
            pred_lookup = {pred['index']: pred for pred in predicted_data}
            
            for measured in self.measured_data:
                data_id = measured['DataID']
                if data_id in pred_lookup:
                    pred = pred_lookup[data_id]
                    calc_point = measured.copy()
                    # Calculate Rho_app,Calc = K * V/I_Calc
                    calc_point['App_Res_calculated'] = measured['K'] * pred['vi_calculated']
                    calc_data.append(calc_point)
            
            if not calc_data:
                print(f"No matching calculated data available for iteration {iteration}")
                continue
                
            filename = f'Calculated apparent resistivity pseudosection — iteration {iteration:02d}.jpg'
            title = f'Calculated Apparent Resistivity Pseudosection — Iteration {iteration:02d}'
            
            self.create_pseudosection_grid(calc_data, 'App_Res_calculated', title, filename, slope_factor)

    def plot_crossplot(self, iteration=None):
        """Plot crossplot of measured vs predicted resistivities.
        If iteration is None, uses the latest iteration.
        """
        if iteration is None:
            # Use latest iteration
            if not self.predicted_data_per_iteration:
                print("No predicted data available for crossplot")
                return
            iteration = max(self.predicted_data_per_iteration.keys())
        
        if iteration not in self.predicted_data_per_iteration:
            print(f"No predicted data available for iteration {iteration}")
            return
            
        predicted_data = self.predicted_data_per_iteration[iteration]
        
        measured_values, predicted_values = [], []
        # Create a lookup for predicted data by index
        pred_lookup = {pred['index']: pred for pred in predicted_data}
        
        for measured in self.measured_data:
            data_id = measured['DataID']
            if data_id in pred_lookup:
                pred = pred_lookup[data_id]
                # Calculate predicted apparent resistivity: Rho_app,Calc = K * V/I_Calc
                pr = measured['K'] * pred['vi_calculated']
                mr = measured['App_Res_measured']
                if mr > 0 and pr > 0:
                    measured_values.append(mr)
                    predicted_values.append(pr)
        
        if not measured_values:
            print("No valid data for crossplot")
            return
            
        measured_values = np.array(measured_values)
        predicted_values = np.array(predicted_values)
        
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.loglog(measured_values, predicted_values, 'bo', alpha=0.6, markersize=4)
        mn, mx = min(measured_values.min(), predicted_values.min()), max(measured_values.max(), predicted_values.max())
        ax.loglog([mn, mx], [mn, mx], 'r-', linewidth=2, label='1:1 Line')
        ax.loglog([mn, mx], [mn*0.9, mx*0.9], 'r--', alpha=0.5, label='±10%')
        ax.loglog([mn, mx], [mn*1.1, mx*1.1], 'r--', alpha=0.5)
        ax.loglog([mn, mx], [mn*0.8, mx*0.8], 'g--', alpha=0.5, label='±20%')
        ax.loglog([mn, mx], [mn*1.2, mx*1.2], 'g--', alpha=0.5)
        ax.set_xlabel('Measured Resistivity (Ohm·m)')
        ax.set_ylabel('Predicted Resistivity (Ohm·m)')
        title = 'Crossplot of Measured vs Predicted Resistivities'
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_aspect('equal')
        r2 = np.corrcoef(measured_values, predicted_values)[0, 1]**2
        ax.text(0.05, 0.95, f'R² = {r2:.3f}', transform=ax.transAxes, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        plt.tight_layout()
        
        filename = 'Crossplot of measured vs calculated resistivities.jpg'
        save_path = self.output_dir / filename
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close()

    def plot_data_misfit_histogram(self, iteration=None):
        """Plot data misfit histogram. If iteration is None, uses the latest iteration."""
        if iteration is None:
            # Use latest iteration
            if not self.predicted_data_per_iteration:
                print("No predicted data available for misfit histogram")
                return
            iteration = max(self.predicted_data_per_iteration.keys())
        
        if iteration not in self.predicted_data_per_iteration:
            print(f"No predicted data available for iteration {iteration}")
            return
            
        predicted_data = self.predicted_data_per_iteration[iteration]
        errors = [abs(pred['vi_error']) for pred in predicted_data]  # Use absolute values
        #errors = [(pred['vi_error']) for pred in predicted_data]  # Use absolute values        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(errors, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
        ax.set_xlabel('Data Misfit (%)')
        ax.set_ylabel('Frequency')
        title = 'Data Misfit Histogram'
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        mean_error = np.mean(errors)
        ax.axvline(mean_error, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_error:.2f}%')
        ax.legend()
        plt.tight_layout()
        
        filename = 'Data misfit histogram.png'
        save_path = self.output_dir / filename
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close()

    def plot_data_misfit_scatter(self, iteration, slope_factor=0.5):
        """Plot data misfit scatter as pseudosection-style plot with error points for a specific iteration."""
        if iteration not in self.predicted_data_per_iteration:
            print(f"No predicted data available for iteration {iteration}")
            return
            
        predicted_data = self.predicted_data_per_iteration[iteration]
        
        misfit_data = []
        # Create a lookup for predicted data by index
        pred_lookup = {pred['index']: pred for pred in predicted_data}
        
        for measured in self.measured_data:
            data_id = measured['DataID']
            if data_id in pred_lookup:
                pred = pred_lookup[data_id]
                misfit_point = measured.copy()
                # Use actual error value (not absolute) to show positive/negative misfits
                misfit_point['misfit_error'] = pred['vi_error']
                misfit_data.append(misfit_point)
        
        if not misfit_data:
            print(f"No valid misfit data for iteration {iteration}")
            return
            
        # Extract coordinates and error values
        x_coords = np.array([d['x_pos'] for d in misfit_data])
        y_coords = np.array([d['pseudodepth'] for d in misfit_data])
        errors = np.array([d['misfit_error'] for d in misfit_data])
        
        # Filter valid data points (allow negative errors now)
        valid_mask = ~np.isnan(errors)
        x_coords, y_coords, errors = x_coords[valid_mask], y_coords[valid_mask], errors[valid_mask]
        
        if len(errors) == 0:
            print(f"No valid misfit data points for iteration {iteration}")
            return
            
        # Use same ROI limits as other pseudosections for consistency
        x_min, x_max, d_min, d_max = self._roi_limits()
        
        # Limit pseudodepth data to the same depth as other plots
        depth_mask = y_coords <= d_max
        x_coords, y_coords, errors = x_coords[depth_mask], y_coords[depth_mask], errors[depth_mask]
        
        if len(errors) == 0:
            print(f"No misfit data within ROI depth limits for iteration {iteration}")
            return
        
        # Apply wedge ROI mask to filter points outside the pyramid shape
        pyramid_mask = self._pyramid_mask(
            x_coords.reshape(-1, 1), 
            y_coords.reshape(-1, 1), 
            x_min, x_max, d_min, d_max, slope_factor
        ).flatten()
        
        x_coords = x_coords[pyramid_mask]
        y_coords = y_coords[pyramid_mask]
        errors = errors[pyramid_mask]
        
        if len(errors) == 0:
            print(f"No misfit data within pyramid ROI for iteration {iteration}")
            return
        
        # Choose figure size from data extent (same as pseudosections)
        xrange = x_max - x_min
        drange = d_max - d_min
        width = 12.0
        height = max(2.2, min(4.0, width * (drange / max(xrange, 1e-9)) * 0.55))

        fig, ax = plt.subplots(figsize=(width, height))
        
        # Color scale for misfit (symmetric around zero to show positive/negative errors)
        vmax = max(abs(np.percentile(errors, 2)), abs(np.percentile(errors, 98)))
        vmin = -vmax  # Symmetric color scale
        
        # Plot scatter points only (no contour lines for error plots)
        scatter = ax.scatter(x_coords, y_coords, c=errors, s=30, cmap=self.current_colormap, 
                           vmin=vmin, vmax=vmax, edgecolors='black', linewidth=0.3, alpha=0.8)
        
        # Electrode markers at the top
        if self.electrode_positions:
            xs = [xy[0] for xy in self.electrode_positions.values()]
            ax.scatter(xs, np.zeros_like(xs), s=16, c="k", marker="s", zorder=3)
        
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(d_min, d_max)
        ax.invert_yaxis()
        
        # Set x-axis at top (same as pseudosections)
        ax.set_xlabel('Distance (m)')
        ax.xaxis.set_label_position('top')
        ax.tick_params(top=True, labeltop=True, bottom=True, labelbottom=False)
        
        ax.set_ylabel('Pseudodepth (m)')
        
        # Title with iteration information
        title = f'Data Misfit Scatter Plot — Iteration {iteration:02d}'
        ax.set_title(title)
        
        # Use auto aspect (same as pseudosections)
        ax.set_aspect('auto')
        
        # Add colorbar - label shows it includes positive and negative errors
        plt.colorbar(scatter, ax=ax, label='Data Misfit (%)')
        plt.tight_layout()
        
        # Save with iteration number in filename
        filename = f'Data misfit scatter — iteration {iteration:02d}.jpg'
        save_path = self.output_dir / filename
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close()

    def plot_all_data_misfit_scatter_iterations(self, slope_factor=0.5):
        """Plot data misfit scatter for all iterations (pseudosection-style with error points)"""
        if not self.predicted_data_per_iteration:
            print("No predicted data available for any iteration")
            return
        
        # Generate misfit scatter plots for all iterations
        iterations = sorted(self.predicted_data_per_iteration.keys())
        print(f"Generating data misfit scatter plots for iterations: {iterations}")
        
        for iteration in iterations:
            self.plot_data_misfit_scatter(iteration, slope_factor)

    # ---------------------------- master runner ----------------------------
    def generate_all_plots(self, inversion_slope_factor=0.5, pseudosection_slope_factor=0.3, 
                          all_iterations_pseudosections=False, colormap='jet', show_contours=False):
        """
        Generate all plots with configurable slope factors for different plot types.
        
        Parameters:
        -----------
        inversion_slope_factor : float, default=0.5
            Controls bottom width for inversion resistivity sections (0.3=narrow, 0.7=wide)
        pseudosection_slope_factor : float, default=0.3  
            Controls bottom width for pseudosection plots (0.3=narrow, 0.7=wide)
        all_iterations_pseudosections : bool, default=False
            If True, generates pseudosections for all iterations. If False, only final iteration.
        colormap : str, default='jet'
            Color scheme for sections and pseudosections. Options: 'cividis', 'viridis', 
            'inferno', 'coolwarm', 'grayscale', 'jet'
        show_contours : bool, default=False
            If True, adds contour lines to sections and pseudosections (not applied to error plots)
        """
        # Set colormap and contour options
        self.set_colormap(colormap)
        self.set_contours(show_contours)
        
        print(f"\nGenerating all plots in directory: {self.output_dir}")
        print(f"Inversion slope factor: {inversion_slope_factor} (bottom width = {inversion_slope_factor*100:.0f}% of surface)")
        print(f"Pseudosection slope factor: {pseudosection_slope_factor} (bottom width = {pseudosection_slope_factor*100:.0f}% of surface)")
        print(f"All iterations pseudosections: {all_iterations_pseudosections}")
        print(f"Colormap: {colormap} ({'gray' if colormap == 'grayscale' else colormap})")
        print(f"Show contours: {show_contours}")
        print("="*60)
        
        # Print detected array types for verification
        if self.measured_data:
            array_types = set(d.get('array_type', 'unknown') for d in self.measured_data)
            print(f"Detected array types: {', '.join(array_types)}")
            print(f"Total measurement points: {len(self.measured_data)}")
            
        if self.predicted_data_per_iteration:
            final_iteration = max(self.predicted_data_per_iteration.keys())
            print(f"Available iterations: {sorted(self.predicted_data_per_iteration.keys())}")
            if all_iterations_pseudosections:
                print(f"Generating pseudosections for all iterations")
            else:
                print(f"Using final iteration ({final_iteration}) for pseudosections")
        
        # Generate pseudosections
        if all_iterations_pseudosections:
            # Generate pseudosections for all iterations
            self.plot_measured_apparent_resistivity_all_iterations(pseudosection_slope_factor)
            self.plot_all_calculated_apparent_resistivity_iterations(pseudosection_slope_factor)
        else:
            # Generate pseudosections only for final iteration
            self.plot_measured_apparent_resistivity(pseudosection_slope_factor)
            # For calculated pseudosections, we need to call the single iteration method for just the final iteration
            if self.predicted_data_per_iteration:
                final_iteration = max(self.predicted_data_per_iteration.keys())
                self.plot_calculated_apparent_resistivity(final_iteration, pseudosection_slope_factor)
        
        # Plot inverted resistivity sections for all iterations (existing method)
        self.plot_and_export_inverted_sections_all_iterations(inversion_slope_factor)
        
        # Plot convergence curve (existing method)
        self.plot_convergence_curve()
        
        # Plot crossplot and histogram for final iteration only
        self.plot_crossplot()
        self.plot_data_misfit_histogram()
        
        # Plot data misfit scatter for all iterations (pseudosection-style with error points)
        self.plot_all_data_misfit_scatter_iterations(pseudosection_slope_factor)
        
        print("="*60)
        print(f"All plots generated successfully in: {self.output_dir.resolve()}")


def main():
    out_file_path = "Dip56_trial1_trial1.out"
    try:
        processor = ERTDataProcessor(out_file_path, "Generated_Images", target_max_depth=None, bottom_width=150)
        
        # Configure slope factors for different plot types:
        # inversion_slope_factor: controls bottom width for resistivity inversion sections
        # pseudosection_slope_factor: controls bottom width for apparent resistivity pseudosections
        # all_iterations_pseudosections: set to True to generate pseudosections for all iterations
        # colormap: choose from 'cividis', 'viridis', 'inferno', 'coolwarm', 'grayscale', 'jet'
        # show_contours: set to True to add contour lines to sections and pseudosections
        
        processor.generate_all_plots(
            inversion_slope_factor=0.5,        # 50% bottom width for inversions
            pseudosection_slope_factor=0.3,    # 30% bottom width for pseudosections (narrower)
            all_iterations_pseudosections=True, # Generate pseudosections for all iterations
            colormap='viridris',                 # Color scheme: 'cividis', 'viridis', 'inferno', 'coolwarm', 'grayscale', 'jet'
            show_contours=True                 # Show contour lines on sections and pseudosections
        )
    except Exception as e:
        print(f"Error processing data: {e}")
        import traceback; traceback.print_exc()

if __name__ == "__main__":
    main()
