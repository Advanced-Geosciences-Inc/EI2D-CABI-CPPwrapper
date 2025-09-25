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
    """Parse EarthImager 2D STG (survey) files"""
    
    @staticmethod  
    def parse_stg(stg_content: str) -> dict:
        """Parse STG file content"""
        lines = stg_content.split('\n')
        electrodes = []
        abmn_data = []
        num_electrodes = None
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            if line.startswith('NumElectrodes='):
                num_electrodes = int(line.split('=')[1])
                continue
            
            # Look for lines with 4 or more integers (ABMN data)
            parts = line.split()
            if len(parts) >= 4:
                try:
                    # Take first 4 integers as A, B, M, N
                    abmn = [int(parts[i]) for i in range(4)]
                    abmn_data.append(abmn)
                except ValueError:
                    continue
        
        # Infer number of electrodes if not specified
        if num_electrodes is None and abmn_data:
            max_elec = max(max(row) for row in abmn_data)
            num_electrodes = max_elec
        
        return {
            "num_electrodes": num_electrodes,
            "abmn_data": abmn_data,
            "num_data": len(abmn_data)
        }