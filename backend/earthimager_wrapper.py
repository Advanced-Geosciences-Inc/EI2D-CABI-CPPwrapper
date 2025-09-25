"""
EarthImager 2D Python wrapper for the C-ABI
Direct interface to the ei2d_core shared library
"""

import ctypes
import numpy as np
from typing import List, Tuple, Optional
import os
from pathlib import Path

class EI2DWrapper:
    """Python wrapper for EarthImager 2D C-ABI"""
    
    def __init__(self, lib_path: Optional[str] = None):
        """Initialize the wrapper with the shared library"""
        if lib_path is None:
            # Default to built library location
            lib_path = "/app/earthimager/cli/build/libei2d_core.so"
        
        if not os.path.exists(lib_path):
            raise FileNotFoundError(f"EI2D library not found at {lib_path}")
            
        self.lib = ctypes.CDLL(lib_path)
        self._setup_function_signatures()
        
    def _setup_function_signatures(self):
        """Set up C function signatures for proper calling"""
        
        # ei2d_InitForwGlobals
        self.lib.ei2d_InitForwGlobals.argtypes = [
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
            ctypes.c_double, ctypes.c_double, ctypes.c_double
        ]
        self.lib.ei2d_InitForwGlobals.restype = None
        
        # ei2d_SetNumParamForward
        self.lib.ei2d_SetNumParamForward.argtypes = [ctypes.c_int32, ctypes.c_int32]
        self.lib.ei2d_SetNumParamForward.restype = None
        
        # ei2d_ForwardFD
        self.lib.ei2d_ForwardFD.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32), 
            ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int32,
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32
        ]
        self.lib.ei2d_ForwardFD.restype = None
        
    def forward_fd_simple(self, n_elec: int = 8, spacing: float = 1.0) -> dict:
        """Run a simple forward FD model (like the CLI example)"""
        
        # Setup basic parameters
        nElec = n_elec
        nInf = 1
        nNx = nElec  # Top row matches electrodes
        nNy = 6      # Depth layers
        
        nNodes = nNx * nNy
        nElem = (nNx - 1) * (nNy - 1)
        
        # Create node coordinates
        nodeX = np.zeros(nNodes, dtype=np.float64)
        nodeY = np.zeros(nNodes, dtype=np.float64)
        
        for j in range(nNy):
            for i in range(nNx):
                idx = j * nNx + i
                nodeX[idx] = i * spacing
                nodeY[idx] = j * spacing
        
        # Homogeneous half-space: rho=100 ohm-m => cond = 0.01 S/m
        cond = np.full(nElem, 0.01, dtype=np.float64)
        
        # Electrode mapping (1-based indices)
        elecNodeID = np.arange(1, nElec + nInf + 1, dtype=np.int32)
        elecNodeID[-1] = 1  # Dummy for infinity electrode
        
        # Infinite electrode indices
        inf = np.array([nElec + nInf], dtype=np.int32)
        
        # Build dipole-dipole survey for n=1,2
        stingCMD = []
        for n in range(1, 3):  # n=1,2
            for i in range(1, nElec - (2 + n) + 1):
                A = i
                B = i + 1
                M = i + 1 + n
                N = i + 2 + n
                stingCMD.extend([A, B, M, N])
        
        stingCMD = np.array(stingCMD, dtype=np.int32)
        nData = len(stingCMD) // 4
        
        # Parameter windows (full grid)
        nParamX, nParamY = 1, 1
        p1 = np.array([1], dtype=np.int32)
        p2 = np.array([nNx], dtype=np.int32)
        q1 = np.array([1], dtype=np.int32)
        q2 = np.array([nNy], dtype=np.int32)
        
        # Initialize engine
        self.lib.ei2d_InitForwGlobals(
            nData,          # NumData
            nElec + nInf,   # NumElectrodes  
            nInf,           # NumInfElectrodes
            nNx,            # NumNodeX
            nNy,            # NumNodeY
            0,              # ForwModMeth (FD)
            0,              # ForwSolver (Cholesky)
            0,              # InvMethod
            0,              # ForwAccuracy
            50,             # ForwCGIter
            0,              # BCType
            1e-6,           # ForwCGResid
            0.0,            # MinTxRxSep
            0.0             # MaxTxRxSep
        )
        
        self.lib.ei2d_SetNumParamForward(nParamX, nParamY)
        
        # Prepare output arrays
        VI = np.zeros(nData, dtype=np.float64)
        jacobian = np.zeros(nData * nParamX * nParamY, dtype=np.float64)
        
        # Call forward modeling
        try:
            self.lib.ei2d_ForwardFD(
                nodeX.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                nodeY.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                cond.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                VI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                jacobian.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                elecNodeID.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                stingCMD.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                p1.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                p2.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                q1.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                q2.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                inf.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                0,              # GetJacobian (no)
                nNodes,         # nNodes
                nElem,          # nElem  
                nData           # nData
            )
            
            return {
                "success": True,
                "nData": nData,
                "VI": [float(x) if np.isfinite(x) else 0.0 for x in VI],
                "stingCMD": stingCMD.reshape(-1, 4).tolist(),
                "nodeX": [float(x) for x in nodeX],
                "nodeY": [float(y) for y in nodeY],
                "conductivity": [float(c) for c in cond],
                "message": f"Forward modeling completed successfully. {nData} data points computed."
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "message": "Forward modeling failed"
            }


class INIParser:
    """Parse and generate EarthImager 2D INI files"""
    
    @staticmethod
    def parse_ini(ini_content: str) -> dict:
        """Parse INI file content into dictionary"""
        result = {}
        current_section = None
        
        for line in ini_content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#') or line.startswith(';'):
                continue
                
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1]
                result[current_section] = {}
            elif '=' in line and current_section:
                key, value = line.split('=', 1)
                result[current_section][key.strip()] = value.strip()
        
        return result
    
    @staticmethod
    def generate_ini(params: dict) -> str:
        """Generate INI file content from parameters"""
        ini_lines = []
        
        for section, values in params.items():
            ini_lines.append(f'[{section}]')
            for key, value in values.items():
                ini_lines.append(f'{key}={value}')
            ini_lines.append('')  # Empty line between sections
            
        return '\n'.join(ini_lines)


class STGParser:
    """Parse EarthImager 2D STG (survey) files - Real AGI format"""
    
    @staticmethod  
    def parse_stg(stg_content: str) -> dict:
        """Parse STG file content - Real AGI format"""
        lines = stg_content.split('\n')
        electrodes = []
        measurements = []
        header_info = {}
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue
                
            # Parse header information
            if i == 0 and "Advanced Geosciences Inc." in line:
                header_info["format"] = "AGI STG"
                header_info["type"] = line.split("Type: ")[-1] if "Type: " in line else "Unknown"
            elif i == 1 and "trimmed data set" in line:
                # Extract record count
                parts = line.split("Records: ")
                if len(parts) > 1:
                    header_info["records"] = int(parts[1].split()[0])
            elif i == 2 and line.startswith("Unit:"):
                header_info["unit"] = line.split(": ")[1]
            
            # Parse measurement data (lines starting with numbers)
            if line and line[0].isdigit():
                try:
                    parts = line.split(',')
                    if len(parts) >= 20:  # AGI STG format has many columns
                        measurement = {
                            "record_id": int(parts[0].strip()),
                            "voltage": float(parts[4]),
                            "current": float(parts[6]),
                            "apparent_resistivity": float(parts[7]),
                            "electrode_a": {'x': float(parts[9]), 'y': float(parts[10]), 'z': float(parts[11])},
                            "electrode_b": {'x': float(parts[12]), 'y': float(parts[13]), 'z': float(parts[14])},
                            "electrode_m": {'x': float(parts[15]), 'y': float(parts[16]), 'z': float(parts[17])},
                            "electrode_n": {'x': float(parts[18]), 'y': float(parts[19]), 'z': float(parts[20])},
                        }
                        measurements.append(measurement)
                        
                        # Extract electrode positions (A, B, M, N as indices)
                        for elec_name, elec_data in [('A', measurement['electrode_a']), 
                                                    ('B', measurement['electrode_b']),
                                                    ('M', measurement['electrode_m']), 
                                                    ('N', measurement['electrode_n'])]:
                            elec_key = f"{elec_data['x']:.1f}_{elec_data['y']:.1f}"
                            if not any(e.get('key') == elec_key for e in electrodes):
                                electrodes.append({
                                    'key': elec_key,
                                    'x': elec_data['x'],
                                    'y': elec_data['y'], 
                                    'z': elec_data['z']
                                })
                except (ValueError, IndexError) as e:
                    continue
        
        # Infer electrode count and spacing
        num_electrodes = len(electrodes) if electrodes else 0
        electrode_spacing = 1.0  # Default
        if len(electrodes) >= 2:
            # Calculate spacing from electrode positions
            x_coords = sorted([e['x'] for e in electrodes])
            if len(x_coords) >= 2:
                electrode_spacing = x_coords[1] - x_coords[0]
        
        return {
            "format": "agi_stg",
            "header_info": header_info,
            "num_electrodes": num_electrodes,
            "num_measurements": len(measurements),
            "electrode_spacing": electrode_spacing,
            "electrodes": electrodes,
            "measurements": measurements[:10],  # First 10 for preview
            "full_measurements": measurements,
            "voltage_range": {
                "min": min([m['voltage'] for m in measurements]) if measurements else 0,
                "max": max([m['voltage'] for m in measurements]) if measurements else 0
            },
            "resistivity_range": {
                "min": min([m['apparent_resistivity'] for m in measurements]) if measurements else 0,
                "max": max([m['apparent_resistivity'] for m in measurements]) if measurements else 0
            }
        }


class MDLParser:
    """Parse EarthImager 2D MDL (model) files"""
    
    @staticmethod
    def parse_mdl(mdl_content: str) -> dict:
        """Parse MDL file content"""
        lines = mdl_content.split('\n')
        sections = {}
        current_section = None
        geometry = []
        commands = []
        model_info = {}
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
                
            if line.startswith(':'):
                current_section = line[1:]
                sections[current_section] = []
                continue
            
            sections.setdefault(current_section, []).append(line)
            
            # Parse specific sections
            if current_section == "Geometry":
                try:
                    parts = line.split(',')
                    if len(parts) >= 3:
                        geometry.append({
                            "electrode": int(parts[0]),
                            "x": float(parts[1]),
                            "y": float(parts[2])
                        })
                except (ValueError, IndexError):
                    continue
            
            elif current_section == "Commands":
                if 'B' in line and 'A' in line and 'M' in line and 'N' in line:
                    # Parse ABMN command format like "2B,1A,3M,4N,R"
                    try:
                        parts = line.replace('B', '').replace('A', '').replace('M', '').replace('N', '').replace('R', '').split(',')
                        if len(parts) >= 4:
                            commands.append({
                                "A": int(parts[1]),
                                "B": int(parts[0]),
                                "M": int(parts[2]),
                                "N": int(parts[3])
                            })
                    except (ValueError, IndexError):
                        continue
            
            elif current_section == "Model":
                if "=" in line:
                    key, value = line.split("=", 1)
                    model_info[key.strip()] = value.strip()
        
        return {
            "format": "agi_mdl",
            "sections": list(sections.keys()),
            "geometry": geometry,
            "commands": commands[:10],  # First 10 for preview
            "full_commands": commands,
            "model_info": model_info,
            "electrode_count": len(geometry),
            "measurement_count": len(commands)
        }


class MODParser:
    """Parse EarthImager 2D MOD (resistivity model) files"""
    
    @staticmethod
    def parse_mod(mod_content: str) -> dict:
        """Parse MOD file content - 3-layer resistivity model"""
        lines = mod_content.split('\n')
        resistivity_blocks = []
        background_resistivity = None
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            
            # Parse background resistivity
            if line.startswith("Background Resistivity"):
                background_resistivity = float(line.split("=")[1].strip())
                continue
            
            # Parse resistivity blocks
            try:
                if ',' in line:
                    parts = line.split(',')
                    if len(parts) >= 4:
                        # Format: Xmin, Ymax, dx, dy with resistivity on previous line
                        x_min = float(parts[0])
                        y_max = float(parts[1])
                        dx = float(parts[2])
                        dy = float(parts[3])
                        
                        # Get resistivity from previous non-empty line
                        continue
                elif line.replace('.', '').isdigit():
                    # This is a resistivity value
                    current_resistivity = float(line)
                    continue
                else:
                    # Try to parse as resistivity block data
                    parts = line.split()
                    if len(parts) >= 4:
                        resistivity_blocks.append({
                            "resistivity": current_resistivity if 'current_resistivity' in locals() else background_resistivity,
                            "x_min": float(parts[0]),
                            "y_max": float(parts[1]), 
                            "dx": float(parts[2]),
                            "dy": float(parts[3])
                        })
            except (ValueError, IndexError, UnboundLocalError):
                continue
        
        # Analyze the 3-layer model structure
        layers = {}
        if resistivity_blocks:
            y_positions = sorted(list(set([block['y_max'] for block in resistivity_blocks])))
            resistivities = list(set([block['resistivity'] for block in resistivity_blocks]))
            
            layers = {
                "num_layers": len(y_positions),
                "y_positions": y_positions,
                "resistivities": resistivities,
                "layer_info": []
            }
            
            # Group blocks by layer
            for i, y_pos in enumerate(y_positions):
                layer_blocks = [b for b in resistivity_blocks if abs(b['y_max'] - y_pos) < 0.001]
                if layer_blocks:
                    layer_resistivities = list(set([b['resistivity'] for b in layer_blocks]))
                    layers["layer_info"].append({
                        "layer": i + 1,
                        "y_position": y_pos,
                        "resistivities": layer_resistivities,
                        "block_count": len(layer_blocks)
                    })
        
        return {
            "format": "agi_mod",
            "background_resistivity": background_resistivity,
            "total_blocks": len(resistivity_blocks),
            "resistivity_blocks": resistivity_blocks[:20],  # First 20 for preview
            "full_blocks": resistivity_blocks,
            "layers": layers,
            "model_summary": {
                "x_range": {
                    "min": min([b['x_min'] for b in resistivity_blocks]) if resistivity_blocks else 0,
                    "max": max([b['x_min'] + b['dx'] for b in resistivity_blocks]) if resistivity_blocks else 0
                },
                "y_range": {
                    "min": min([b['y_max'] - b['dy'] for b in resistivity_blocks]) if resistivity_blocks else 0,
                    "max": max([b['y_max'] for b in resistivity_blocks]) if resistivity_blocks else 0
                }
            }
        }