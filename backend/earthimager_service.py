"""
EarthImager 2D Web Service
FastAPI endpoints for geophysical modeling
"""

from fastapi import HTTPException, UploadFile, File
from typing import Optional, Dict, Any
import tempfile
import os
from pathlib import Path
import json
import traceback

from .earthimager_wrapper import EI2DWrapper, INIParser, STGParser

class EarthImagerService:
    """Service class for EarthImager 2D operations"""
    
    def __init__(self):
        self.wrapper = None
        self._initialize_wrapper()
        
    def _initialize_wrapper(self):
        """Initialize the EI2D wrapper"""
        try:
            self.wrapper = EI2DWrapper()
        except Exception as e:
            print(f"Warning: Could not initialize EI2D wrapper: {e}")
            self.wrapper = None
    
    async def run_forward_modeling(self, 
                                 n_electrodes: int = 8, 
                                 electrode_spacing: float = 1.0,
                                 resistivity: float = 100.0) -> Dict[str, Any]:
        """Run forward modeling with simple parameters"""
        
        if not self.wrapper:
            raise HTTPException(status_code=500, detail="EI2D engine not available")
        
        try:
            # Convert resistivity to conductivity
            conductivity = 1.0 / resistivity if resistivity > 0 else 0.01
            
            result = self.wrapper.forward_fd_simple(
                n_elec=n_electrodes, 
                spacing=electrode_spacing
            )
            
            if not result.get("success"):
                raise HTTPException(status_code=500, detail=result.get("error", "Forward modeling failed"))
            
            # Format results for web display
            formatted_result = {
                "success": True,
                "parameters": {
                    "n_electrodes": n_electrodes,
                    "electrode_spacing": electrode_spacing,
                    "resistivity": resistivity,
                    "conductivity": conductivity
                },
                "results": {
                    "num_data_points": result["nData"],
                    "vi_data": result["VI"][:10],  # First 10 for preview
                    "total_vi_count": len(result["VI"]),
                    "survey_config": result["stingCMD"][:5],  # First 5 ABMN configs
                    "mesh_info": {
                        "nodes_x": len(set(result["nodeX"])),
                        "nodes_y": len(set(result["nodeY"])),
                        "total_nodes": len(result["nodeX"])
                    }
                },
                "full_data": result  # Include all data for download
            }
            
            return formatted_result
            
        except Exception as e:
            error_msg = f"Forward modeling error: {str(e)}"
            print(f"Error: {error_msg}")
            print(f"Traceback: {traceback.format_exc()}")
            raise HTTPException(status_code=500, detail=error_msg)
    
    async def process_ini_file(self, ini_content: str) -> Dict[str, Any]:
        """Process uploaded INI file"""
        try:
            parsed_ini = INIParser.parse_ini(ini_content)
            
            # Extract key parameters
            forward_params = parsed_ini.get("Forward", {})
            initial_params = parsed_ini.get("Initial", {})
            resinv_params = parsed_ini.get("ResInv", {})
            
            return {
                "success": True,
                "sections": list(parsed_ini.keys()),
                "forward_method": forward_params.get("ForwModMeth", "0"),
                "boundary_condition": forward_params.get("BCType", "0"),
                "max_iterations": resinv_params.get("MaxNumInvIter", "20"),
                "lagrange_multiplier": resinv_params.get("Lagrange", "10"),
                "parsed_data": parsed_ini
            }
            
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"INI parsing error: {str(e)}")
    
    async def process_stg_file(self, stg_content: str) -> Dict[str, Any]:
        """Process uploaded STG file"""
        try:
            parsed_stg = STGParser.parse_stg(stg_content)
            
            return {
                "success": True,
                "num_electrodes": parsed_stg["num_electrodes"],
                "num_measurements": parsed_stg["num_data"],
                "survey_preview": parsed_stg["abmn_data"][:10],  # First 10 measurements
                "electrode_range": {
                    "min": min(min(row) for row in parsed_stg["abmn_data"]) if parsed_stg["abmn_data"] else 1,
                    "max": max(max(row) for row in parsed_stg["abmn_data"]) if parsed_stg["abmn_data"] else 1
                },
                "parsed_data": parsed_stg
            }
            
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"STG parsing error: {str(e)}")
    
    async def generate_ini_config(self, params: Dict[str, Any]) -> str:
        """Generate INI configuration from web parameters"""
        try:
            # Default INI structure based on the example
            ini_config = {
                "Forward": {
                    "ForwModMeth": str(params.get("forward_method", 0)),
                    "ForwSolver": str(params.get("forward_solver", 0)),
                    "BCType": str(params.get("bc_type", 0)),
                    "ForwAccuracy": str(params.get("forward_accuracy", 1)),
                    "ForwCGIter": str(params.get("cg_iterations", 100)),
                    "ForwCGResid": str(params.get("cg_residual", "1E-6"))
                },
                "ResInv": {
                    "MaxNumInvIter": str(params.get("max_iterations", 20)),
                    "MaxRMSRes": str(params.get("max_rms", 2)),
                    "Lagrange": str(params.get("lagrange", 10)),
                    "StartRes": str(params.get("start_resistivity", 1)),
                    "MinResis": str(params.get("min_resistivity", 1)),
                    "MaxResis": str(params.get("max_resistivity", 100000))
                },
                "Initial": {
                    "MinVoltage": str(params.get("min_voltage", 0.2)),
                    "MinVoverI": str(params.get("min_v_over_i", 0.0005)),
                    "MinAppRes": str(params.get("min_app_res", 1)),
                    "MaxAppRes": str(params.get("max_app_res", 10000))
                }
            }
            
            return INIParser.generate_ini(ini_config)
            
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"INI generation error: {str(e)}")
    
    async def health_check(self) -> Dict[str, Any]:
        """Check service health"""
        return {
            "service": "EarthImager 2D Web Interface",
            "status": "healthy",
            "ei2d_engine": "available" if self.wrapper else "unavailable",
            "version": "1.0.0"
        }