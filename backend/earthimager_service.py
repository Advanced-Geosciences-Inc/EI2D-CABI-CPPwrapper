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

from earthimager_wrapper import EI2DWrapper, INIParser, STGParser, MDLParser, MODParser, EI2DRealDataProcessor

class EarthImagerService:
    """Service class for EarthImager 2D operations"""
    
    def __init__(self):
        self.wrapper = None
        self.real_processor = None
        self._initialize_services()
        
    def _initialize_services(self):
        """Initialize EI2D services"""
        try:
            self.wrapper = EI2DWrapper()
        except Exception as e:
            print(f"Warning: Could not initialize EI2D wrapper: {e}")
            self.wrapper = None
            
        try:
            self.real_processor = EI2DRealDataProcessor()
        except Exception as e:
            print(f"Warning: Could not initialize EI2D real processor: {e}")
            self.real_processor = None
    
    async def run_forward_modeling(self, 
                                 n_electrodes: int = 8, 
                                 electrode_spacing: float = 1.0,
                                 resistivity: float = 100.0) -> Dict[str, Any]:
        """Run forward modeling with simple parameters"""
        
        try:
            # For now, create a mock result to test the interface
            # TODO: Replace with actual EI2D wrapper call once debugged
            
            # Generate mock data that mimics EI2D output
            import random
            import math
            
            n_data = (n_electrodes - 3) * 2  # Dipole-dipole for n=1,2
            
            # Mock V/I data (typical geophysical values)
            vi_data = []
            for i in range(n_data):
                # Simulate realistic V/I values based on resistivity
                base_vi = 1.0 / (resistivity * 0.01)  # Basic relationship
                noise = random.uniform(0.8, 1.2)  # Add some variation
                vi_data.append(base_vi * noise)
            
            # Mock ABMN survey configuration
            survey_config = []
            for n in range(1, 3):  # n=1,2
                for i in range(1, n_electrodes - (2 + n) + 1):
                    A = i
                    B = i + 1
                    M = i + 1 + n
                    N = i + 2 + n
                    survey_config.append([A, B, M, N])
            
            # Conductivity
            conductivity = 1.0 / resistivity if resistivity > 0 else 0.01
            
            formatted_result = {
                "success": True,
                "parameters": {
                    "n_electrodes": n_electrodes,
                    "electrode_spacing": electrode_spacing,
                    "resistivity": resistivity,
                    "conductivity": conductivity
                },
                "results": {
                    "num_data_points": n_data,
                    "vi_data": vi_data[:10],  # First 10 for preview
                    "total_vi_count": len(vi_data),
                    "survey_config": survey_config[:5],  # First 5 ABMN configs
                    "mesh_info": {
                        "nodes_x": n_electrodes,
                        "nodes_y": 6,
                        "total_nodes": n_electrodes * 6
                    }
                },
                "full_data": {
                    "VI": vi_data,
                    "survey_config": survey_config,
                    "message": f"Mock forward modeling completed successfully. {n_data} data points computed.",
                    "note": "This is mock data for interface testing. Real EI2D engine integration pending."
                }
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
    
    async def run_real_forward_modeling(self, ini_content: str, stg_content: str) -> Dict[str, Any]:
        """Run real forward modeling using INI and STG files"""
        
        if not self.real_processor:
            raise HTTPException(status_code=500, detail="EI2D real processor not available")
        
        try:
            result = self.real_processor.process_ini_stg_files(ini_content, stg_content)
            
            if not result.get("success"):
                raise HTTPException(status_code=500, detail=result.get("error", "Real forward modeling failed"))
            
            return {
                "success": True,
                "method": result.get("method", "unknown"),
                "parameters": result.get("parameters", {}),
                "results": result.get("results", {}),
                "mesh": result.get("mesh", {}),
                "message": result.get("message", "Forward modeling completed"),
                "note": result.get("note", "")
            }
            
        except Exception as e:
            error_msg = f"Real forward modeling error: {str(e)}"
            print(f"Error: {error_msg}")
            print(f"Traceback: {traceback.format_exc()}")
            raise HTTPException(status_code=500, detail=error_msg)
        """Process uploaded STG file - Real AGI format"""
        try:
            parsed_stg = STGParser.parse_stg(stg_content)
            
            return {
                "success": True,
                "format": parsed_stg.get("format"),
                "header_info": parsed_stg.get("header_info", {}),
                "num_electrodes": parsed_stg["num_electrodes"],
                "num_measurements": parsed_stg["num_measurements"],
                "electrode_spacing": parsed_stg.get("electrode_spacing", 1.0),
                "measurement_preview": parsed_stg["measurements"][:5],  # First 5 measurements
                "voltage_range": parsed_stg.get("voltage_range", {}),
                "resistivity_range": parsed_stg.get("resistivity_range", {}),
                "electrodes": parsed_stg.get("electrodes", []),
                "parsed_data": parsed_stg
            }
            
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"STG parsing error: {str(e)}")
    
    async def process_mdl_file(self, mdl_content: str) -> Dict[str, Any]:
        """Process uploaded MDL file"""
        try:
            parsed_mdl = MDLParser.parse_mdl(mdl_content)
            
            return {
                "success": True,
                "format": parsed_mdl.get("format"),
                "sections": parsed_mdl["sections"],
                "electrode_count": parsed_mdl["electrode_count"],
                "measurement_count": parsed_mdl["measurement_count"],
                "geometry_preview": parsed_mdl["geometry"][:10],  # First 10 electrodes
                "commands_preview": parsed_mdl["commands"][:5],   # First 5 ABMN commands
                "model_info": parsed_mdl.get("model_info", {}),
                "parsed_data": parsed_mdl
            }
            
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"MDL parsing error: {str(e)}")
    
    async def process_mod_file(self, mod_content: str) -> Dict[str, Any]:
        """Process uploaded MOD file - 3-layer resistivity model"""
        try:
            parsed_mod = MODParser.parse_mod(mod_content)
            
            return {
                "success": True,
                "format": parsed_mod.get("format"),
                "background_resistivity": parsed_mod["background_resistivity"],
                "total_blocks": parsed_mod["total_blocks"],
                "layers": parsed_mod.get("layers", {}),
                "model_summary": parsed_mod.get("model_summary", {}),
                "resistivity_preview": parsed_mod["resistivity_blocks"][:10],  # First 10 blocks
                "parsed_data": parsed_mod
            }
            
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"MOD parsing error: {str(e)}")
    
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