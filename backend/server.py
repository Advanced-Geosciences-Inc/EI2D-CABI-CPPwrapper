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
@api_router.post("/earthimager/forward-model")
async def run_forward_model(params: ForwardModelParams):
    """Run forward modeling with specified parameters"""
    return await ei_service.run_forward_modeling(
        n_electrodes=params.n_electrodes,
        electrode_spacing=params.electrode_spacing,
        resistivity=params.resistivity
    )

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
async def upload_stg_file(file: UploadFile = File(...)):
    """Upload and process STG survey file"""
    if not file.filename.endswith('.stg'):
        raise HTTPException(status_code=400, detail="File must be an STG file")
    
    content = await file.read()
    stg_content = content.decode('utf-8')
    
    result = await ei_service.process_stg_file(stg_content)
    result["filename"] = file.filename
    return result

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
