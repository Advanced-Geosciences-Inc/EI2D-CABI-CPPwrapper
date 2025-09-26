# ROI Analyzer for EarthImager 2D OUT files
# Extracts proper ROI information from MESH SIZE and NODE LOCATION sections

import re
import numpy as np
from pathlib import Path

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

def parse_mesh_size_section(lines):
    """Parse the MESH SIZE section to get proper dimensions"""
    mesh_info = {}
    
    _, mesh_block = extract_section(lines, ";------ MESH SIZE ------", ";------")
    
    if mesh_block:
        print("=== MESH SIZE SECTION ===")
        for line in mesh_block:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            print(f"  {line}")
            
            # Extract key mesh parameters
            if "Number of nodes" in line and "=" in line:
                mesh_info['total_nodes'] = int(line.split("=")[1].strip())
            elif "Number of nodes in X" in line and "=" in line:
                mesh_info['nodes_x'] = int(line.split("=")[1].strip())
            elif "Number of nodes in Y" in line and "=" in line:
                mesh_info['nodes_y'] = int(line.split("=")[1].strip())
            elif "Number of elements" in line and "=" in line:
                mesh_info['total_elements'] = int(line.split("=")[1].strip())
            elif "Number of elements in X" in line and "=" in line:
                mesh_info['elements_x'] = int(line.split("=")[1].strip())
            elif "Number of elements in Y" in line and "=" in line:
                mesh_info['elements_y'] = int(line.split("=")[1].strip())
            elif "Number of parameters" in line and "=" in line:
                mesh_info['total_parameters'] = int(line.split("=")[1].strip())
            elif "Number of parameters in X" in line and "=" in line:
                mesh_info['param_x'] = int(line.split("=")[1].strip())
            elif "Number of parameters in Y" in line and "=" in line:
                mesh_info['param_y'] = int(line.split("=")[1].strip())
    
    return mesh_info

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

def parse_electrode_locations(lines):
    """Parse electrode locations"""
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

def analyze_roi_information(out_file):
    """Comprehensive ROI analysis"""
    
    if not out_file.exists():
        raise FileNotFoundError(f"Cannot find {out_file}")
    
    lines = out_file.read_text(encoding="utf-8", errors="ignore").splitlines()
    print(f"Analyzing: {out_file.resolve()}")
    print(f"File contains {len(lines)} lines\n")
    
    # 1. Parse mesh size information
    mesh_info = parse_mesh_size_section(lines)
    print(f"\n=== MESH INFORMATION ===")
    for key, value in mesh_info.items():
        print(f"  {key}: {value}")
    
    # 2. Parse electrode locations
    electrodes = parse_electrode_locations(lines)
    if electrodes:
        x_coords = [e['x'] for e in electrodes]
        y_coords = [e['y'] for e in electrodes]
        print(f"\n=== ELECTRODE INFORMATION ===")
        print(f"  Number of electrodes: {len(electrodes)}")
        print(f"  X range: {min(x_coords):.1f} to {max(x_coords):.1f} m")
        print(f"  Y range: {min(y_coords):.1f} to {max(y_coords):.1f} m")
        print(f"  Survey length: {max(x_coords) - min(x_coords):.1f} m")
        
        electrode_spacing = (max(x_coords) - min(x_coords)) / (len(electrodes) - 1) if len(electrodes) > 1 else 1.0
        print(f"  Electrode spacing: {electrode_spacing:.1f} m")
    
    # 3. Parse NODE LOCATION sections
    _, x_block = extract_section(lines, ";------ NODE LOCATION IN X (m) ------", ";------")
    _, y_block = extract_section(lines, ";------ NODE LOCATION IN Y (m) ------", ";------")
    
    if x_block and y_block:
        x_nodes = parse_node_coordinates(x_block)
        y_nodes = parse_node_coordinates(y_block)
        
        print(f"\n=== NODE LOCATION ANALYSIS ===")
        print(f"  X nodes: {len(x_nodes)} nodes")
        print(f"  X range: {np.min(x_nodes):.1f} to {np.max(x_nodes):.1f} m")
        print(f"  X span: {np.max(x_nodes) - np.min(x_nodes):.1f} m")
        
        print(f"  Y nodes: {len(y_nodes)} nodes")  
        print(f"  Y range: {np.min(y_nodes):.1f} to {np.max(y_nodes):.1f} m")
        print(f"  Y span: {np.max(y_nodes) - np.min(y_nodes):.1f} m")
        
        # Calculate element centroids
        x_cent = 0.5 * (x_nodes[:-1] + x_nodes[1:])
        depth_cent = -0.5 * (y_nodes[:-1] + y_nodes[1:])  # Negative Y becomes positive depth
        
        print(f"\n=== ELEMENT CENTROID ANALYSIS ===")
        print(f"  X elements: {len(x_cent)} elements")
        print(f"  X centroid range: {np.min(x_cent):.1f} to {np.max(x_cent):.1f} m")
        
        print(f"  Y elements: {len(depth_cent)} elements")
        print(f"  Depth centroid range: {np.min(depth_cent):.1f} to {np.max(depth_cent):.1f} m")
        
        # Compare with electrode array
        if electrodes:
            elec_x_min = min(x_coords)
            elec_x_max = max(x_coords)
            
            print(f"\n=== ROI COMPARISON ===")
            print(f"  Electrode array X: [{elec_x_min:.1f}, {elec_x_max:.1f}]")
            print(f"  Full mesh X: [{np.min(x_cent):.1f}, {np.max(x_cent):.1f}]")
            print(f"  Full mesh depth: [0.0, {np.max(depth_cent):.1f}]")
            
            # Calculate EarthImager 2D style ROI (likely based on investigation depth)
            survey_length = elec_x_max - elec_x_min
            investigation_depth = survey_length * 0.15  # Conservative estimate
            
            print(f"\n=== SUGGESTED EI2D ROI ===")
            print(f"  Survey length: {survey_length:.1f} m")
            print(f"  Estimated investigation depth: {investigation_depth:.1f} m")
            print(f"  Suggested X ROI: [{elec_x_min:.1f}, {elec_x_max:.1f}]")
            print(f"  Suggested depth ROI: [0.0, {investigation_depth:.1f}]")
            
            # Find elements within this ROI
            x_roi_mask = (x_cent >= elec_x_min) & (x_cent <= elec_x_max)
            depth_roi_mask = (depth_cent >= 0.0) & (depth_cent <= investigation_depth)
            
            roi_x_elements = np.sum(x_roi_mask)
            roi_depth_elements = np.sum(depth_roi_mask)
            
            print(f"  ROI elements: {roi_x_elements} × {roi_depth_elements}")
            
            return {
                'electrodes': electrodes,
                'mesh_info': mesh_info,
                'x_nodes': x_nodes,
                'y_nodes': y_nodes,
                'x_centroids': x_cent,
                'depth_centroids': depth_cent,
                'suggested_roi': {
                    'x_min': elec_x_min,
                    'x_max': elec_x_max,
                    'depth_min': 0.0,
                    'depth_max': investigation_depth
                }
            }
    
    return None

# Test with multiple files
def main():
    # Test files (adjust paths as needed)
    test_files = [
        Path("Dip56_trial1_trial1.out"),
        Path("BeaverHoles_trial5.out"), 
        Path("earthimager_results_fixed.out")
    ]
    
    for out_file in test_files:
        if out_file.exists():
            print("=" * 80)
            try:
                analyze_roi_information(out_file)
            except Exception as e:
                print(f"Error analyzing {out_file}: {e}")
            print("=" * 80)
            print()
    
if __name__ == "__main__":
    main()