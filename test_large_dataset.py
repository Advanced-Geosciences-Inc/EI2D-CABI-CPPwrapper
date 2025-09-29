#!/usr/bin/env python3
"""
Test JSON serialization fix with larger datasets
Simulates the Amistad dataset mentioned in the review request (56 electrodes, 424 measurements)
"""

import requests
import json
import tempfile
import os

# Configuration
BACKEND_URL = "https://ei2d-integration.preview.emergentagent.com/api"

def create_large_stg_file(num_electrodes=56, num_measurements=424):
    """Create a larger STG file similar to Amistad dataset"""
    
    stg_content = []
    stg_content.append("Advanced Geosciences Inc. (AGI) Sting/SuperSting measured data (*.stg)   Type: XYZ")
    stg_content.append(f"A large dataset test. Version: 1.0.0. Records: {num_measurements}")
    stg_content.append("Unit: meter")
    
    # Generate measurements with realistic electrode spacing
    electrode_spacing = 1.0
    
    measurement_id = 1
    for n in range(1, 8):  # Different n-spacings (1-7)
        for a in range(1, num_electrodes - (2 + n) + 1):
            if measurement_id > num_measurements:
                break
                
            A = (a - 1) * electrode_spacing
            B = a * electrode_spacing  
            M = (a + n) * electrode_spacing
            N = (a + n + 1) * electrode_spacing
            
            # Generate realistic voltage and apparent resistivity
            voltage = -0.1 - (measurement_id * 0.001)  # Decreasing voltage
            current = 1000  # mA
            app_res = 100 + (measurement_id % 100)  # Varying resistivity
            
            line = f"{measurement_id:d}    ,USER,   20250924, 10:07:23,{voltage:.5E},   0, {current:d}, {app_res:.5E}, large_dataset_test, {A:.4f},  0.0000,  0.0000,  {B:.4f},  0.0000,  0.0000,  {M:.4f},  0.0000,  0.0000,  {N:.4f},  0.0000,  0.0000"
            stg_content.append(line)
            
            measurement_id += 1
            
        if measurement_id > num_measurements:
            break
    
    return '\n'.join(stg_content)

def create_large_ini_file():
    """Create INI file for large dataset"""
    
    ini_content = """[Startup]
Xpos=10
Ypos=9
Width=1262
Height=1018
State=0
FontName=Times New Roman
FontSize=10
LastProjDir=C:\\AGI\\EarthImager2D\\Demo
LastLineDir=C:\\AGI\\EarthImager2D\\Demo\\large-dataset
LastTrialDir=C:\\AGI\\EarthImager2D\\Demo\\large-dataset\\trial1
LastDataFile=C:\\AGI\\EarthImager2D\\Demo\\large-dataset.cmd
MetersFeet=Meters
JobCodeK=Job Code
ProjSiteK=Project Site
ApprovedByK=Approved By
SurveyDateK=Survey Date
InstrumentK=Instrument
InvSoftwareK=Software
DataFileK=Data File
EmployerName=Advanced Geosciences, Inc.
JobCode=AGI #002
ProjSite=Large Dataset Test
ApprovedBy=
Instrument=SuperSting R8
IPUnit=1

[Initial]
MinVoltage=0.2
MinVoverI=0.0005
MinAppRes=1
MaxAppRes=10000
MaxRepeatErr=2
MaxRecipErr=2
RemoveNegERT=1
KeepAllData=0
ExcludeEdges=0
InvMethod=0
VerticalAxis=0
YAxis=0
MinElecSepX=0.003
MinElecSepZ=0.003
ContourLines=0
BlnAspectRatio=1

[Forward]
ForwModMeth=1
ForwSolver=0
BCType=0
ForwAccuracy=1
Anisotropy=1
ElemDivision=2
ThickFactor=1.1
DepthFactor=1.1
ForwCGIter=100
ForwCGResid=1E-6

[ResInv]
MaxNumInvIter=15
MaxRMSRes=3
MinErrReduction=2
StopOnMaxIter=1
StopOnMaxRMS=0
StopOnMinErrDiff=0
DataReweight=0
UseRecipErr=0
StopOnL2Norm=1
Lagrange=5
RoughConditioner=0.2
ResStartID=1
StartRes=100
MinResis=10
MaxResis=10000
ParameterWidth=1
ParameterHeight=1
HVRoughRatio=0.5
ResNoisePC=2
MaxNumIterInvCG=10
DampFactorRes=10
EpsilonD=1
EpsilonM=1
QuasiNewton=20

[IPInv]
IPInvMethod=0
MaxNumInvIterIP=8
MaxRMSIP=2
MinErrReductionIP=5
StopOnMaxIterIP=1
StopOnMaxRMSIP=1
StopOnMinErrDiffIP=1
DataReweightIP=0
StopOnL2NormIP=0
IPNoisePC=2
IPNoiseMin=-0.2
IPNoiseMax=0.5
MinIPCorr=0.7
RemoveNegIP=0
Positivity=0
LagrangeIP=10
DampFactorIP=10
IPWinStart=1
IPWinEnd=6
IPStartID=1
StartIP=0.01

[Terrain]
MeshXformMethod=0

[CRP]
CRPNumElec=120
CRPOverlapPC=60
"""
    return ini_content

def test_large_dataset_json_serialization():
    """Test JSON serialization with larger dataset"""
    
    print("🔍 Testing JSON serialization with large dataset (56 electrodes, 424 measurements)...")
    
    # Create large dataset files
    stg_content = create_large_stg_file(56, 424)
    ini_content = create_large_ini_file()
    
    print(f"   Generated STG file: {len(stg_content)} characters")
    print(f"   Generated INI file: {len(ini_content)} characters")
    
    # Test forward modeling with large dataset
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.stg', delete=False) as stg_f:
            stg_f.write(stg_content)
            stg_path = stg_f.name
            
        with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as ini_f:
            ini_f.write(ini_content)
            ini_path = ini_f.name
        
        # Test forward modeling
        print("📡 Testing forward modeling with large dataset...")
        with open(ini_path, 'rb') as ini_f, open(stg_path, 'rb') as stg_f:
            files = {
                'ini_file': ('large_dataset.ini', ini_f, 'text/plain'),
                'stg_file': ('large_dataset.stg', stg_f, 'text/plain')
            }
            
            response = requests.post(
                f"{BACKEND_URL}/earthimager/forward-model-real",
                files=files,
                timeout=60
            )
        
        if response.status_code == 200:
            try:
                data = response.json()
                print("✅ Large dataset forward modeling JSON parsing successful!")
                
                # Check data sizes
                vi_data = data.get("results", {}).get("vi_data", [])
                apparent_res = data.get("results", {}).get("apparent_resistivities", [])
                
                print(f"   VI data points: {len(vi_data)}")
                print(f"   Apparent resistivity points: {len(apparent_res)}")
                
                # Check for NaN/infinity
                import math
                vi_invalid = any(math.isnan(x) or math.isinf(x) for x in vi_data if isinstance(x, (int, float)))
                ar_invalid = any(math.isnan(x) or math.isinf(x) for x in apparent_res if isinstance(x, (int, float)))
                
                print(f"   VI data has NaN/inf: {vi_invalid}")
                print(f"   Apparent resistivities have NaN/inf: {ar_invalid}")
                
                forward_ok = not (vi_invalid or ar_invalid)
                
            except json.JSONDecodeError as e:
                print(f"❌ Large dataset forward modeling JSON parsing failed: {e}")
                forward_ok = False
        else:
            print(f"❌ Large dataset forward modeling HTTP error: {response.status_code}")
            forward_ok = False
        
        # Test inversion with large dataset
        print("📡 Testing inversion workflow with large dataset...")
        with open(ini_path, 'rb') as ini_f, open(stg_path, 'rb') as stg_f:
            files = {
                'ini_file': ('large_dataset.ini', ini_f, 'text/plain'),
                'stg_file': ('large_dataset.stg', stg_f, 'text/plain')
            }
            
            response = requests.post(
                f"{BACKEND_URL}/earthimager/run-inversion",
                files=files,
                timeout=120  # Longer timeout for large dataset
            )
        
        if response.status_code == 200:
            try:
                data = response.json()
                print("✅ Large dataset inversion JSON parsing successful!")
                
                # Check data sizes
                results = data.get("results", {})
                resistivity_model = results.get("resistivity_model", [])
                calculated_data = results.get("calculated_data", [])
                data_residuals = results.get("data_residuals", [])
                
                print(f"   Resistivity model points: {len(resistivity_model)}")
                print(f"   Calculated data points: {len(calculated_data)}")
                print(f"   Data residuals points: {len(data_residuals)}")
                
                # Check for NaN/infinity
                import math
                rm_invalid = any(math.isnan(x) or math.isinf(x) for x in resistivity_model if isinstance(x, (int, float)))
                cd_invalid = any(math.isnan(x) or math.isinf(x) for x in calculated_data if isinstance(x, (int, float)))
                dr_invalid = any(math.isnan(x) or math.isinf(x) for x in data_residuals if isinstance(x, (int, float)))
                
                print(f"   Resistivity model has NaN/inf: {rm_invalid}")
                print(f"   Calculated data has NaN/inf: {cd_invalid}")
                print(f"   Data residuals have NaN/inf: {dr_invalid}")
                
                inversion_ok = not (rm_invalid or cd_invalid or dr_invalid)
                
            except json.JSONDecodeError as e:
                print(f"❌ Large dataset inversion JSON parsing failed: {e}")
                print(f"   Response content (first 500 chars): {response.text[:500]}")
                inversion_ok = False
        else:
            print(f"❌ Large dataset inversion HTTP error: {response.status_code}")
            print(f"   Response: {response.text[:500]}")
            inversion_ok = False
        
        # Cleanup
        os.unlink(stg_path)
        os.unlink(ini_path)
        
        return forward_ok and inversion_ok
        
    except Exception as e:
        print(f"❌ Large dataset test failed: {e}")
        return False

def main():
    """Main test execution"""
    print("🧪 Large Dataset JSON Serialization Test")
    print("=" * 50)
    
    success = test_large_dataset_json_serialization()
    
    print("\n" + "=" * 50)
    if success:
        print("🎉 Large dataset JSON serialization test PASSED!")
        print("   No 'Out of range float values are not JSON compliant' errors with large datasets")
    else:
        print("💥 Large dataset JSON serialization test FAILED!")
        print("   JSON serialization issues detected with large datasets")
    
    return success

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)