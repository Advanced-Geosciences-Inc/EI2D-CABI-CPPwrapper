# Fixed EI2D OUT file parser - handles variable grid dimensions
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

# Configurable input file
OUT = Path("earthimager_results.out")  # <-- adjust if needed

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
            if 2 <= aspect_ratio <= 10:  # Reasonable for geophysical surveys
                possible_factors.append((n_y, n_x, aspect_ratio))
    
    if possible_factors:
        # Choose the factor with aspect ratio closest to 5 (common in ERT)
        best_factor = min(possible_factors, key=lambda x: abs(x[2] - 5.0))
        return best_factor[0], best_factor[1]
    
    return None, None

# -------------------------
# Mask helper (correct orientation)  
# -------------------------
def apply_corner_masks(Xg, Dg, Zi, x_min, x_max, d_max, angle_deg=45):
    """Mask bottom corners that slope DOWN from the top corners."""
    theta = np.deg2rad(angle_deg)
    L_cut = d_max / np.tan(theta) if np.tan(theta) != 0 else d_max

    # Left curtain
    xl = Xg - x_min
    left_line = (d_max / L_cut) * xl
    left_region = (xl >= 0) & (xl <= L_cut) & (Dg > left_line)

    # Right curtain  
    xr = (x_max - Xg)
    right_line = (d_max / L_cut) * xr
    right_region = (xr >= 0) & (xr <= L_cut) & (Dg > right_line)

    Zm = Zi.copy()
    Zm[left_region] = np.nan
    Zm[right_region] = np.nan
    return Zm

# -------------------------
# Enhanced renderer
# -------------------------
def render_resistivity_roi_smart(
    x_centroids, depth_centroids, Z_roi,
    x_limits=None, d_limits=None,
    grid_dx=0.5, grid_dz=0.5,
    cmap="jet", levels=25,
    gaussian_sigma=0.8,
    mask_corners=True, mask_angle_deg=45,
    vmin=None, vmax=None,
    title="Resistivity Model",
    save_path=None, show=True
):
    """Enhanced renderer with automatic limit detection"""
    
    # Auto-detect limits if not provided
    if x_limits is None:
        x_pad = (np.max(x_centroids) - np.min(x_centroids)) * 0.1
        x_limits = (np.min(x_centroids) - x_pad, np.max(x_centroids) + x_pad)
    
    if d_limits is None:
        d_pad = (np.max(depth_centroids) - np.min(depth_centroids)) * 0.1  
        d_limits = (np.min(depth_centroids) - d_pad, np.max(depth_centroids) + d_pad)
    
    # Dense plotting grid
    xr = np.arange(x_limits[0], x_limits[1] + grid_dx*0.5, grid_dx)
    dr = np.arange(d_limits[0], d_limits[1] + grid_dz*0.5, grid_dz)
    Xg, Dg = np.meshgrid(xr, dr)

    # Interpolate from centroids  
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

    # Optional corner masking
    if mask_corners:
        Zi = apply_corner_masks(Xg, Dg, Zi, x_limits[0], x_limits[1], 
                               d_limits[1], angle_deg=mask_angle_deg)

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    cf = ax.contourf(Xg, Dg, Zi, levels=levels, cmap=cmap,
                     vmin=vmin, vmax=vmax, antialiased=True)
    ax.set_aspect('equal', adjustable='box')
    ax.invert_yaxis()
    ax.set_xlabel("Offset (m)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(title)
    plt.colorbar(cf, ax=ax, label="Resistivity (Ohm·m)")
    plt.tight_layout()

    if save_path:
        sp = Path(save_path)
        sp.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(sp, dpi=200)
        print(f"Saved: {sp.resolve()}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return Xg, Dg, Zi

# -------------------------
# Enhanced main function
# -------------------------
def main():
    if not OUT.exists():
        raise SystemExit(f"Cannot find {OUT.resolve()}")

    lines = OUT.read_text(encoding="utf-8", errors="ignore").splitlines()
    print(f"Processing: {OUT.resolve()}")
    print(f"File contains {len(lines)} lines")

    # Method 1: Try to parse NODE LOCATION sections (preferred)
    print("\n=== METHOD 1: Parsing NODE LOCATION sections ===")
    _, x_block = extract_section(lines, ";------ NODE LOCATION IN X (m) ------", ";------")
    _, y_block = extract_section(lines, ";------ NODE LOCATION IN Y (m) ------", ";------")

    if x_block and y_block:
        x_nodes = parse_node_coordinates(x_block)
        y_nodes = parse_node_coordinates(y_block)
        
        if len(x_nodes) > 1 and len(y_nodes) > 1:
            print(f"Found NODE sections: {len(x_nodes)} X-nodes, {len(y_nodes)} Y-nodes")
            
            # Calculate centroids
            x_cent = 0.5 * (x_nodes[:-1] + x_nodes[1:])
            depth_cent = -0.5 * (y_nodes[:-1] + y_nodes[1:])
            
            n_ex = len(x_cent)
            n_ey = len(depth_cent)
            expected_elements = n_ex * n_ey
            
            print(f"Grid dimensions: {n_ex} × {n_ey} = {expected_elements} elements")
            
            # Find resistivity data
            last_iter = max((i for i, ln in enumerate(lines) if ";------ Iteration" in ln), default=-1)
            if last_iter >= 0:
                _, res_block = extract_section(
                    lines[last_iter:], ";-Resistivity in Ohm-m in the elemental sequential order", ";-Sensitivity"
                )
                
                if res_block:
                    vals = parse_numbers_sci(res_block)
                    print(f"Found {len(vals)} resistivity values")
                    
                    if len(vals) >= expected_elements:
                        # Truncate to expected size
                        vals = vals[:expected_elements]
                        Z = vals.reshape((n_ey, n_ex), order="F")
                        
                        print("✅ Successfully parsed with NODE LOCATION method")
                        method_success = True
                    else:
                        print(f"❌ Not enough values: need {expected_elements}, got {len(vals)}")
                        method_success = False
                else:
                    print("❌ No resistivity block found")
                    method_success = False
            else:
                print("❌ No iteration sections found")
                method_success = False
        else:
            print("❌ Invalid node coordinates")
            method_success = False
    else:
        print("❌ NODE LOCATION sections not found")
        method_success = False

    # Method 2: Infer grid from resistivity data (fallback)
    if not method_success:
        print("\n=== METHOD 2: Inferring grid from resistivity data ===")
        
        # Find any resistivity block
        last_iter = max((i for i, ln in enumerate(lines) if ";------ Iteration" in ln), default=-1)
        if last_iter >= 0:
            _, res_block = extract_section(
                lines[last_iter:], ";-Resistivity in Ohm-m in the elemental sequential order", ";-Sensitivity"
            )
        else:
            # Try simpler pattern 
            _, res_block = extract_section(lines, "Resistivity in Ohm-m", "Sensitivity")
        
        if res_block:
            n_ey_inferred, n_ex_inferred = infer_grid_from_resistivity_data(res_block)
            
            if n_ey_inferred and n_ex_inferred:
                print(f"Inferred grid: {n_ex_inferred} × {n_ey_inferred}")
                
                vals = parse_numbers_sci(res_block)
                expected_elements = n_ex_inferred * n_ey_inferred
                
                if len(vals) >= expected_elements:
                    vals = vals[:expected_elements]
                    Z = vals.reshape((n_ey_inferred, n_ex_inferred), order="F")
                    
                    # Create synthetic coordinates
                    x_cent = np.linspace(0, 100, n_ex_inferred)
                    depth_cent = np.linspace(0, 50, n_ey_inferred)
                    
                    n_ex, n_ey = n_ex_inferred, n_ey_inferred
                    method_success = True
                    print("✅ Successfully inferred grid dimensions")
                else:
                    print(f"❌ Still not enough values: need {expected_elements}, got {len(vals)}")
                    method_success = False
            else:
                print("❌ Could not infer reasonable grid dimensions")
                method_success = False
        else:
            print("❌ No resistivity data found")
            method_success = False

    if not method_success:
        raise SystemExit("❌ Failed to parse resistivity data with any method")

    # Apply ROI filtering (adjust based on your data)
    print(f"\n=== RENDERING ===")
    print(f"Final grid: {n_ex} × {n_ey}")
    print(f"X range: {np.min(x_cent):.1f} to {np.max(x_cent):.1f}")
    print(f"Depth range: {np.min(depth_cent):.1f} to {np.max(depth_cent):.1f}")
    print(f"Resistivity range: {np.min(Z):.1f} to {np.max(Z):.1f} Ω·m")

    # Simple ROI (use all data if dimensions are reasonable)
    if n_ex <= 200 and n_ey <= 50:  # Reasonable sizes
        Z_roi = Z
        x_roi = x_cent  
        d_roi = depth_cent
    else:
        # Apply ROI filtering for large grids
        x_limits = (np.min(x_cent), np.max(x_cent))
        d_limits = (0, min(50, np.max(depth_cent)))
        
        xmask = (x_cent >= x_limits[0]) & (x_cent <= x_limits[1])
        dmask = (depth_cent >= d_limits[0]) & (depth_cent <= d_limits[1])
        
        Z_roi = Z[dmask, :][:, xmask]
        x_roi = x_cent[xmask]
        d_roi = depth_cent[dmask]
    
    # Render with smart limits
    render_resistivity_roi_smart(
        x_centroids=x_roi,
        depth_centroids=d_roi,
        Z_roi=Z_roi,
        cmap="RdBu_r",
        levels=25,
        gaussian_sigma=0.8,
        mask_corners=True,
        mask_angle_deg=45,
        save_path=f"fixed_{OUT.stem}_plot.png",
        title=f"Resistivity Model - {OUT.name}",
        show=True
    )

if __name__ == "__main__":
    main()