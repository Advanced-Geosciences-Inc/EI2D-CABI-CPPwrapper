# ROI-aware EI2D OUT file parser with proper geophysical filtering
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

# Configurable input file
OUT = Path("Dip56_trial1_trial1.out")  # <-- adjust if needed

# -------------------------
# Enhanced parsing helpers
# -------------------------
def extract_section(lines, start_marker, end_marker=None, start_idx=0):
    """Extract section with improved error handling"""
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

def parse_electrode_locations(lines):
    """Parse electrode locations to determine survey extent"""
    electrodes = []
    _, elec_block = extract_section(lines, ";------ ELECTRODE LOCATIONS ------", ";------")
    
    if elec_block:
        for line in elec_block[2:]:  # Skip header
            parts = line.strip().split(',')
            if len(parts) >= 3:
                try:
                    idx = int(parts[0].strip())
                    x = float(parts[1].strip())
                    y = float(parts[2].strip())
                    electrodes.append({'idx': idx, 'x': x, 'y': y})
                except ValueError:
                    continue
    
    return electrodes

def calculate_geophysical_roi(electrodes, mesh_info=None):
    """Calculate geophysically meaningful ROI matching EarthImager 2D criteria"""
    if not electrodes:
        return None
        
    x_coords = [e['x'] for e in electrodes]
    
    # Survey extent (electrode array bounds)
    x_min = min(x_coords)
    x_max = max(x_coords)
    survey_length = x_max - x_min
    
    # EarthImager 2D ROI criteria (based on analysis):
    # - Lateral: exactly electrode array bounds (no padding)
    # - Depth: ~15% of survey length for dipole-dipole arrays
    roi_x_min = x_min
    roi_x_max = x_max
    roi_depth_max = survey_length * 0.15  # EarthImager 2D investigation depth
    
    roi = {
        'x_min': roi_x_min,
        'x_max': roi_x_max, 
        'depth_min': 0.0,
        'depth_max': roi_depth_max,
        'survey_length': survey_length,
        'electrode_spacing': survey_length / (len(electrodes) - 1) if len(electrodes) > 1 else 1.0
    }
    
    print(f"EarthImager 2D ROI: X=[{roi_x_min:.1f}, {roi_x_max:.1f}], Depth=[0.0, {roi_depth_max:.1f}]")
    print(f"Survey length: {survey_length:.1f} m, Investigation depth: {roi_depth_max:.1f} m")
    
    return roi

def parse_node_coordinates(raw_lines):
    """Parse 'Node, value' lines (skip the 2-row header)."""
    coords = []
    for ln in raw_lines[2:]:  # Skip header
        parts = ln.strip().split(',')
        if len(parts) >= 2:
            try:
                coords.append(float(parts[1].strip()))
            except ValueError:
                pass
    return np.array(coords, float)

def parse_numbers_sci(raw_lines):
    """Parse numbers including scientific notation from a lines block."""
    pat = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    vals = []
    for ln in raw_lines:
        if ln.strip().startswith(";"):
            continue
        # Replace commas with spaces for better parsing
        clean_line = ln.replace(",", " ")
        tokens = re.findall(pat, clean_line)
        vals.extend(map(float, tokens))
    return np.array(vals, float)

def infer_grid_from_resistivity_data(res_lines):
    """Try to infer grid dimensions from resistivity data structure"""
    if not res_lines:
        return None, None
    
    # Count total values
    vals = parse_numbers_sci(res_lines)
    total_vals = len(vals)
    
    # Try common aspect ratios and factorizations
    possible_factors = []
    for n_y in range(5, int(np.sqrt(total_vals)) + 20):
        if total_vals % n_y == 0:
            n_x = total_vals // n_y
            # Prefer reasonable aspect ratios (more width than depth)
            aspect_ratio = n_x / n_y
            if 2 <= aspect_ratio <= 12:  # Reasonable for geophysical surveys
                possible_factors.append((n_y, n_x, aspect_ratio))
    
    if possible_factors:
        # Choose the factor with aspect ratio closest to 5 (common in ERT)
        best_factor = min(possible_factors, key=lambda x: abs(x[2] - 5.0))
        return best_factor[0], best_factor[1]
    
    return None, None

# -------------------------
# ROI-aware mask helper
# -------------------------
def apply_geophysical_masks(Xg, Dg, Zi, roi, angle_deg=45):
    """Apply geophysically meaningful masks"""
    
    # 1. Remove data above surface (z < 0)
    surface_mask = Dg < 0
    Zi[surface_mask] = np.nan
    
    # 2. Remove data outside lateral survey extent
    lateral_mask = (Xg < roi['x_min']) | (Xg > roi['x_max'])
    Zi[lateral_mask] = np.nan
    
    # 3. Remove data below depth of investigation
    depth_mask = Dg > roi['depth_max']
    Zi[depth_mask] = np.nan
    
    # 4. Apply corner masks (standard ERT presentation)
    theta = np.deg2rad(angle_deg)
    L_cut = roi['depth_max'] / np.tan(theta) if np.tan(theta) != 0 else roi['depth_max']

    # Left curtain
    xl = Xg - roi['x_min']
    left_line = (roi['depth_max'] / L_cut) * xl
    left_region = (xl >= 0) & (xl <= L_cut) & (Dg > left_line)

    # Right curtain  
    xr = (roi['x_max'] - Xg)
    right_line = (roi['depth_max'] / L_cut) * xr
    right_region = (xr >= 0) & (xr <= L_cut) & (Dg > right_line)

    Zi[left_region] = np.nan
    Zi[right_region] = np.nan
    
    return Zi

# -------------------------
# ROI-aware renderer
# -------------------------
def render_resistivity_roi_geophysical(
    x_centroids, depth_centroids, Z_roi, roi,
    grid_dx=0.5, grid_dz=0.2,
    cmap="RdBu_r", levels=25,
    gaussian_sigma=0.6,
    mask_corners=True, mask_angle_deg=45,
    vmin=None, vmax=None,
    title="Resistivity Model (ROI)",
    save_path=None, show=True
):
    """ROI-aware renderer matching EarthImager 2D style"""
    
    # Use ROI limits for plotting grid
    x_limits = (roi['x_min'], roi['x_max'])
    d_limits = (roi['depth_min'], roi['depth_max'])
    
    print(f"Rendering ROI: X=[{x_limits[0]:.1f}, {x_limits[1]:.1f}], Depth=[{d_limits[0]:.1f}, {d_limits[1]:.1f}]")
    
    # Dense plotting grid within ROI
    xr = np.arange(x_limits[0], x_limits[1] + grid_dx*0.5, grid_dx)
    dr = np.arange(d_limits[0], d_limits[1] + grid_dz*0.5, grid_dz)
    Xg, Dg = np.meshgrid(xr, dr)

    # Interpolate from centroids within ROI
    Xc, Dc = np.meshgrid(x_centroids, depth_centroids)
    pts = np.column_stack([Xc.ravel(), Dc.ravel()])
    vals = Z_roi.ravel()

    # Multi-method interpolation
    Zi = griddata(pts, vals, (Xg, Dg), method="cubic")
    if np.isnan(Zi).any():
        Zi2 = griddata(pts, vals, (Xg, Dg), method="linear")
        Zi[np.isnan(Zi)] = Zi2[np.isnan(Zi)]
    if np.isnan(Zi).any():
        Zi3 = griddata(pts, vals, (Xg, Dg), method="nearest")
        Zi[np.isnan(Zi)] = Zi3[np.isnan(Zi)]

    # Optional Gaussian blur
    if gaussian_sigma and gaussian_sigma > 0.0:
        Zi = gaussian_filter(Zi, sigma=gaussian_sigma, mode="nearest")

    # Apply geophysical masking
    if mask_corners:
        Zi = apply_geophysical_masks(Xg, Dg, Zi, roi, angle_deg=mask_angle_deg)

    # Plot with proper aspect ratio
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Use resistivity-appropriate color limits if not specified
    if vmin is None or vmax is None:
        finite_vals = Zi[np.isfinite(Zi)]
        if len(finite_vals) > 0:
            if vmin is None:
                vmin = np.percentile(finite_vals, 5)   # Clip extreme low values
            if vmax is None:
                vmax = np.percentile(finite_vals, 95)  # Clip extreme high values
    
    cf = ax.contourf(Xg, Dg, Zi, levels=levels, cmap=cmap,
                     vmin=vmin, vmax=vmax, antialiased=True, extend='both')
    
    ax.set_aspect('equal', adjustable='box')
    ax.invert_yaxis()
    ax.set_xlabel("Distance (m)", fontsize=12)
    ax.set_ylabel("Depth (m)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add colorbar with proper label
    cbar = plt.colorbar(cf, ax=ax, label="Resistivity (Ω·m)", shrink=0.8)
    cbar.ax.tick_params(labelsize=10)
    
    # Set proper axis limits to ROI
    ax.set_xlim(x_limits)
    ax.set_ylim(d_limits[1], d_limits[0])  # Inverted for depth
    
    plt.tight_layout()

    if save_path:
        sp = Path(save_path)
        sp.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(sp, dpi=300, bbox_inches='tight')
        print(f"Saved: {sp.resolve()}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return Xg, Dg, Zi

# -------------------------
# Enhanced main function with ROI filtering
# -------------------------
def main():
    if not OUT.exists():
        raise SystemExit(f"Cannot find {OUT.resolve()}")

    lines = OUT.read_text(encoding="utf-8", errors="ignore").splitlines()
    print(f"Processing: {OUT.resolve()}")
    print(f"File contains {len(lines)} lines")

    # Step 1: Parse electrode locations to determine ROI
    print("\n=== STEP 1: Determine EarthImager 2D ROI ===")
    electrodes = parse_electrode_locations(lines)
    if electrodes:
        roi = calculate_geophysical_roi(electrodes)
        if not roi:
            raise SystemExit("❌ Could not calculate meaningful ROI")
    else:
        print("⚠️  No electrode locations found, using default ROI")
        roi = {'x_min': 0, 'x_max': 50, 'depth_min': 0, 'depth_max': 7.5, 'survey_length': 50}

    # Step 2: Parse grid and resistivity data
    print("\n=== STEP 2: Parse resistivity model ===")
    method_success = False
    
    # Method 1: Try NODE LOCATION sections (preferred)
    _, x_block = extract_section(lines, ";------ NODE LOCATION IN X (m) ------", ";------")
    _, y_block = extract_section(lines, ";------ NODE LOCATION IN Y (m) ------", ";------")

    if x_block and y_block:
        x_nodes = parse_node_coordinates(x_block)
        y_nodes = parse_node_coordinates(y_block)
        
        if len(x_nodes) > 1 and len(y_nodes) > 1:
            print(f"Found NODE sections: {len(x_nodes)} X-nodes, {len(y_nodes)} Y-nodes")
            
            # Calculate centroids
            x_cent = 0.5 * (x_nodes[:-1] + x_nodes[1:])
            depth_cent = -0.5 * (y_nodes[:-1] + y_nodes[1:])  # Negative Y becomes positive depth
            
            n_ex = len(x_cent)
            n_ey = len(depth_cent)
            
            # Find resistivity data
            last_iter = max((i for i, ln in enumerate(lines) if ";------ Iteration" in ln), default=-1)
            if last_iter >= 0:
                _, res_block = extract_section(
                    lines[last_iter:], ";-Resistivity in Ohm-m in the elemental sequential order", ";-Sensitivity"
                )
                
                if res_block:
                    vals = parse_numbers_sci(res_block)
                    expected_elements = n_ex * n_ey
                    
                    if len(vals) >= expected_elements:
                        vals = vals[:expected_elements]
                        Z = vals.reshape((n_ey, n_ex), order="F")
                        method_success = True
                        print(f"✅ Parsed {n_ex}×{n_ey} grid with NODE sections")
    
    # Method 2: Infer from resistivity data (fallback)
    if not method_success:
        print("Falling back to grid inference...")
        
        last_iter = max((i for i, ln in enumerate(lines) if ";------ Iteration" in ln), default=-1)
        if last_iter >= 0:
            _, res_block = extract_section(
                lines[last_iter:], ";-Resistivity in Ohm-m in the elemental sequential order", ";-Sensitivity"
            )
        else:
            _, res_block = extract_section(lines, "Resistivity in Ohm-m", "Sensitivity")
        
        if res_block:
            n_ey_inferred, n_ex_inferred = infer_grid_from_resistivity_data(res_block)
            
            if n_ey_inferred and n_ex_inferred:
                vals = parse_numbers_sci(res_block)
                expected_elements = n_ex_inferred * n_ey_inferred
                
                if len(vals) >= expected_elements:
                    vals = vals[:expected_elements]
                    Z = vals.reshape((n_ey_inferred, n_ex_inferred), order="F")
                    
                    # Create coordinates based on ROI
                    x_cent = np.linspace(roi['x_min'], roi['x_max'], n_ex_inferred)
                    depth_cent = np.linspace(roi['depth_min'], roi['depth_max'], n_ey_inferred)
                    
                    n_ex, n_ey = n_ex_inferred, n_ey_inferred
                    method_success = True
                    print(f"✅ Inferred {n_ex}×{n_ey} grid from resistivity data")

    if not method_success:
        raise SystemExit("❌ Failed to parse resistivity data")

    # Step 3: Apply ROI filtering to model
    print(f"\n=== STEP 3: Apply ROI filtering ===")
    print(f"Full model: {n_ex} × {n_ey} elements")
    print(f"X range: {np.min(x_cent):.1f} to {np.max(x_cent):.1f}")
    print(f"Depth range: {np.min(depth_cent):.1f} to {np.max(depth_cent):.1f}")
    print(f"Resistivity range: {np.min(Z):.1f} to {np.max(Z):.1f} Ω·m")

    # Filter model to ROI
    xmask = (x_cent >= roi['x_min']) & (x_cent <= roi['x_max'])
    dmask = (depth_cent >= roi['depth_min']) & (depth_cent <= roi['depth_max'])
    
    Z_roi = Z[dmask, :][:, xmask]
    x_roi = x_cent[xmask]
    d_roi = depth_cent[dmask]
    
    print(f"ROI model: {len(x_roi)} × {len(d_roi)} elements")
    print(f"ROI resistivity range: {np.min(Z_roi):.1f} to {np.max(Z_roi):.1f} Ω·m")
    
    # Step 4: Render with proper ROI
    print(f"\n=== STEP 4: Render geophysical model ===")
    render_resistivity_roi_geophysical(
        x_centroids=x_roi,
        depth_centroids=d_roi,
        Z_roi=Z_roi,
        roi=roi,
        cmap="RdBu_r",
        levels=30,
        gaussian_sigma=0.6,
        mask_corners=True,
        mask_angle_deg=45,
        save_path=f"ei2d_matched_{OUT.stem}.png",
        title=f"ERT Model (EI2D ROI) - {OUT.name}",
        show=True
    )

if __name__ == "__main__":
    main()