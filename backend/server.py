from fastapi import FastAPI, APIRouter, HTTPException, UploadFile, File, Form
from fastapi.responses import FileResponse
from dotenv import load_dotenv
from starlette.middleware.cors import CORSMiddleware
from motor.motor_asyncio import AsyncIOMotorClient
import os
import logging
from pathlib import Path
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
import uuid
from datetime import datetime
import tempfile
import json

# Import EarthImager service
from earthimager_service import EarthImagerService

ROOT_DIR = Path(__file__).parent
load_dotenv(ROOT_DIR / '.env')

# MongoDB connection
mongo_url = os.environ['MONGO_URL']
client = AsyncIOMotorClient(mongo_url)
db = client[os.environ['DB_NAME']]

# Initialize EarthImager service
ei_service = EarthImagerService()

# Create the main app without a prefix
app = FastAPI(title="EarthImager 2D Web Interface", version="1.0.0")

# Create a router with the /api prefix
api_router = APIRouter(prefix="/api")


# Define Models
class StatusCheck(BaseModel):
    id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    client_name: str
    timestamp: datetime = Field(default_factory=datetime.utcnow)

class StatusCheckCreate(BaseModel):
    client_name: str

class ForwardModelParams(BaseModel):
    n_electrodes: int = Field(default=8, ge=4, le=120)
    electrode_spacing: float = Field(default=1.0, gt=0)
    resistivity: float = Field(default=100.0, gt=0)

class INIConfigParams(BaseModel):
    forward_method: int = Field(default=0, ge=0, le=1)  # 0=FD, 1=FE
    forward_solver: int = Field(default=0, ge=0, le=1)  # 0=Cholesky, 1=CG
    bc_type: int = Field(default=0, ge=0, le=1)         # 0=Dirichlet, 1=Mixed
    max_iterations: int = Field(default=20, ge=1, le=100)
    lagrange: float = Field(default=10.0, gt=0)
    start_resistivity: float = Field(default=1.0, gt=0)
    min_resistivity: float = Field(default=1.0, gt=0)
    max_resistivity: float = Field(default=100000.0, gt=0)

# Basic routes
@api_router.get("/")
async def root():
    return {"message": "EarthImager 2D Web Interface", "version": "1.0.0"}

@api_router.get("/health")
async def health_check():
    return await ei_service.health_check()

# EarthImager routes
@api_router.post("/earthimager/run-inversion")
async def run_full_inversion(ini_file: UploadFile = File(...), stg_file: UploadFile = File(...)):
    """Run complete EarthImager 2D inversion workflow"""
    
    # Validate file extensions
    if not ini_file.filename.endswith('.ini'):
        raise HTTPException(status_code=400, detail="First file must be an INI file")
    if not stg_file.filename.endswith('.stg'):
        raise HTTPException(status_code=400, detail="Second file must be an STG file")
    
    try:
        # Read file contents
        ini_content = (await ini_file.read()).decode('utf-8')
        stg_content = (await stg_file.read()).decode('utf-8')
        
        # Process with EI2D inversion workflow
        result = await ei_service.run_full_inversion(ini_content, stg_content)
        result["input_files"] = {
            "ini_file": ini_file.filename,
            "stg_file": stg_file.filename
        }
        
        return result
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Inversion workflow failed: {str(e)}")

@api_router.post("/earthimager/download-out-file")
async def download_out_file(content: str = Form(...)):
    """Download generated OUT file"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.out', delete=False) as f:
        f.write(content)
        temp_path = f.name
    
    return FileResponse(
        path=temp_path,
        filename="earthimager_inversion_results.out",
        media_type="text/plain"
    )

@api_router.post("/earthimager/forward-model-real")
async def run_real_forward_model(ini_file: UploadFile = File(...), stg_file: UploadFile = File(...)):
    """Run real forward modeling using uploaded INI and STG files"""
    
    # Validate file extensions
    if not ini_file.filename.endswith('.ini'):
        raise HTTPException(status_code=400, detail="First file must be an INI file")
    if not stg_file.filename.endswith('.stg'):
        raise HTTPException(status_code=400, detail="Second file must be an STG file")
    
    try:
        # Read file contents
        ini_content = (await ini_file.read()).decode('utf-8')
        stg_content = (await stg_file.read()).decode('utf-8')
        
        # Process with real EI2D engine
        result = await ei_service.run_real_forward_modeling(ini_content, stg_content)
        result["input_files"] = {
            "ini_file": ini_file.filename,
            "stg_file": stg_file.filename
        }
        
        return result
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Real forward modeling failed: {str(e)}")

@api_router.post("/earthimager/forward-model")
async def run_forward_model(params: ForwardModelParams):
    """Run forward modeling with specified parameters"""
    return await ei_service.run_forward_modeling(
        n_electrodes=params.n_electrodes,
        electrode_spacing=params.electrode_spacing,
        resistivity=params.resistivity
    )

@api_router.post("/earthimager/validate-data")
async def validate_data_flow(ini_file: UploadFile = File(...), stg_file: UploadFile = File(...)):
    """Validate and inspect data processing without running full calculations"""
    
    try:
        # Read file contents
        ini_content = (await ini_file.read()).decode('utf-8')
        stg_content = (await stg_file.read()).decode('utf-8')
        
        print(f"\n=== DATA VALIDATION REPORT ===")
        print(f"INI File: {ini_file.filename} ({len(ini_content)} chars)")
        print(f"STG File: {stg_file.filename} ({len(stg_content)} chars)")
        
        # Parse files and show detailed breakdown
        from earthimager_wrapper import INIParser, STGParser
        
        ini_data = INIParser.parse_ini(ini_content)
        stg_data = STGParser.parse_stg(stg_content)
        
        # Extract specific values for validation
        forward_params = ini_data.get("Forward", {})
        resinv_params = ini_data.get("ResInv", {})
        
        validation_report = {
            "success": True,
            "files": {
                "ini_sections": list(ini_data.keys()),
                "stg_measurements": len(stg_data["full_measurements"]),
                "stg_electrodes": len(stg_data["electrodes"])
            },
            "parsed_parameters": {
                "forward_method": forward_params.get("ForwModMeth", "Not found"),
                "forward_solver": forward_params.get("ForwSolver", "Not found"),
                "bc_type": forward_params.get("BCType", "Not found"),
                "max_iterations": resinv_params.get("MaxNumInvIter", "Not found"),
                "lagrange": resinv_params.get("Lagrange", "Not found"),
                "target_rms": resinv_params.get("MaxRMSRes", "Not found")
            },
            "survey_validation": {
                "electrode_spacing": stg_data.get("electrode_spacing", 0),
                "voltage_range": stg_data.get("voltage_range", {}),
                "resistivity_range": stg_data.get("resistivity_range", {}),
                "measurement_sample": stg_data["measurements"][:3] if stg_data["measurements"] else []
            },
            "data_integrity": {
                "all_measurements_have_coordinates": all(
                    "electrode_a" in m and "electrode_b" in m and "electrode_m" in m and "electrode_n" in m
                    for m in stg_data["full_measurements"][:10]
                ),
                "voltage_current_present": all(
                    "voltage" in m and "current" in m and "apparent_resistivity" in m
                    for m in stg_data["full_measurements"][:10]
                ),
                "coordinate_consistency": len(set([
                    f"{m['electrode_a']['x']:.1f}" for m in stg_data["full_measurements"][:10]
                ])) > 1  # Check if we have multiple electrode positions
            }
        }
        
        print(f"Forward Method: {validation_report['parsed_parameters']['forward_method']}")
        print(f"Survey measurements: {validation_report['files']['stg_measurements']}")
        print(f"Electrode spacing: {validation_report['survey_validation']['electrode_spacing']}")
        
        return validation_report
        
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "message": "Validation failed"
        }

@api_router.post("/earthimager/debug-processing")
async def debug_processing_steps(ini_file: UploadFile = File(...), stg_file: UploadFile = File(...)):
    """Debug each step of the processing pipeline with detailed logging"""
    
    try:
        ini_content = (await ini_file.read()).decode('utf-8')
        stg_content = (await stg_file.read()).decode('utf-8')
        
        debug_info = {"steps": [], "success": True}
        
        # Step 1: File parsing
        debug_info["steps"].append({
            "step": 1,
            "name": "File Parsing",
            "status": "completed",
            "details": {
                "ini_size": len(ini_content),
                "stg_size": len(stg_content)
            }
        })
        
        # Step 2: Data extraction
        from earthimager_wrapper import EI2DRealDataProcessor
        processor = EI2DRealDataProcessor()
        
        from earthimager_wrapper import INIParser, STGParser
        ini_data = INIParser.parse_ini(ini_content)
        stg_data = STGParser.parse_stg(stg_content)
        
        debug_info["steps"].append({
            "step": 2, 
            "name": "Data Extraction",
            "status": "completed",
            "details": {
                "ini_sections": len(ini_data),
                "measurements_found": len(stg_data["full_measurements"]),
                "electrodes_found": len(stg_data["electrodes"])
            }
        })
        
        # Step 3: Mesh generation
        mesh_result = processor._generate_inversion_mesh(stg_data["electrodes"], stg_data["full_measurements"])
        
        debug_info["steps"].append({
            "step": 3,
            "name": "Mesh Generation", 
            "status": "completed",
            "details": {
                "mesh_nodes_x": mesh_result["nodes_x"],
                "mesh_nodes_y": mesh_result["nodes_y"],
                "total_nodes": mesh_result["total_nodes"],
                "mesh_elements": mesh_result["total_elements"],
                "parameters": mesh_result["num_parameters"]
            }
        })
        
        # Step 4: Compare with expected values (your toy example reference)
        expected_electrodes = 14
        expected_measurements = 74
        
        validation_checks = {
            "electrode_count_match": len(stg_data["electrodes"]) == expected_electrodes,
            "measurement_count_match": len(stg_data["full_measurements"]) == expected_measurements,
            "mesh_reasonable": 300 <= mesh_result["total_nodes"] <= 1000,
            "parameters_reasonable": 100 <= mesh_result["num_parameters"] <= 500
        }
        
        debug_info["steps"].append({
            "step": 4,
            "name": "Validation Checks",
            "status": "completed" if all(validation_checks.values()) else "warnings",
            "details": validation_checks
        })
        
        return debug_info
        
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "debug_info": debug_info if 'debug_info' in locals() else {}
        }

@api_router.post("/earthimager/upload-ini")
async def upload_ini_file(file: UploadFile = File(...)):
    """Upload and process INI configuration file"""
    if not file.filename.endswith('.ini'):
        raise HTTPException(status_code=400, detail="File must be an INI file")
    
    content = await file.read()
    ini_content = content.decode('utf-8')
    
    result = await ei_service.process_ini_file(ini_content)
    result["filename"] = file.filename
    return result

@api_router.post("/earthimager/upload-mdl")
async def upload_mdl_file(file: UploadFile = File(...)):
    """Upload and process MDL model file"""
    if not file.filename.endswith('.mdl'):
        raise HTTPException(status_code=400, detail="File must be an MDL file")
    
    content = await file.read()
    mdl_content = content.decode('utf-8')
    
    result = await ei_service.process_mdl_file(mdl_content)
    result["filename"] = file.filename
    return result

@api_router.post("/earthimager/upload-mod")
async def upload_mod_file(file: UploadFile = File(...)):
    """Upload and process MOD resistivity model file"""
    if not file.filename.endswith('.mod'):
        raise HTTPException(status_code=400, detail="File must be a MOD file")
    
    content = await file.read()
    mod_content = content.decode('utf-8')
    
    result = await ei_service.process_mod_file(mod_content)
    result["filename"] = file.filename
    return result
@api_router.post("/earthimager/upload-stg")
async def upload_stg_file(file: UploadFile = File(...)):
    """Upload and process STG survey file"""
    try:
        print(f"Received STG file: {file.filename}, size: {file.size}")
        
        if not file.filename.endswith('.stg'):
            raise HTTPException(status_code=400, detail="File must be an STG file")
        
        content = await file.read()
        print(f"Read {len(content)} bytes")
        
        stg_content = content.decode('utf-8')
        print(f"Decoded to {len(stg_content)} characters")
        print(f"First 200 chars: {stg_content[:200]}")
        
        result = await ei_service.process_stg_file(stg_content)
        result["filename"] = file.filename
        return result
        
    except UnicodeDecodeError as e:
        print(f"Unicode decode error: {e}")
        raise HTTPException(status_code=400, detail=f"File encoding error: {str(e)}")
    except Exception as e:
        print(f"STG upload error: {e}")
        raise HTTPException(status_code=500, detail=f"STG processing failed: {str(e)}")

@api_router.post("/earthimager/generate-ini")
async def generate_ini_config(params: INIConfigParams):
    """Generate INI configuration file from parameters"""
    ini_content = await ei_service.generate_ini_config(params.dict())
    
    # Save to temporary file and return download link
    with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as f:
        f.write(ini_content)
        temp_path = f.name
    
    return {
        "success": True,
        "message": "INI file generated successfully",
        "content": ini_content,
        "download_url": f"/api/earthimager/download-ini/{Path(temp_path).name}"
    }

@api_router.get("/earthimager/download-ini/{filename}")
async def download_ini_file(filename: str):
    """Download generated INI file"""
    file_path = Path(tempfile.gettempdir()) / filename
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")
    
    return FileResponse(
        path=file_path,
        filename=f"earthimager_config_{filename}",
        media_type="text/plain"
    )

# Legacy status check routes (keep for compatibility)
@api_router.post("/status", response_model=StatusCheck)
async def create_status_check(input: StatusCheckCreate):
    status_dict = input.dict()
    status_obj = StatusCheck(**status_dict)
    _ = await db.status_checks.insert_one(status_obj.dict())
    return status_obj

@api_router.get("/status", response_model=List[StatusCheck])
async def get_status_checks():
    status_checks = await db.status_checks.find().to_list(1000)
    return [StatusCheck(**status_check) for status_check in status_checks]

# Include the router in the main app
app.include_router(api_router)

app.add_middleware(
    CORSMiddleware,
    allow_credentials=True,
    allow_origins=os.environ.get('CORS_ORIGINS', '*').split(','),
    allow_methods=["*"],
    allow_headers=["*"],
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@app.on_event("shutdown")
async def shutdown_db_client():
    client.close()
