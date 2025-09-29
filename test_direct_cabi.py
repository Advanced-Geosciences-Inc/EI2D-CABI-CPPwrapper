#!/usr/bin/env python3
"""
Direct test of EI2D C-ABI to verify array bounds fix
"""

import sys
import os
sys.path.append('/app/backend')

from earthimager_wrapper import EI2DRealDataProcessor, INIParser, STGParser
from pathlib import Path

def test_direct_cabi():
    """Test direct C-ABI calls to verify array bounds fix"""
    
    print("🔧 DIRECT C-ABI ARRAY BOUNDS TEST")
    print("=" * 50)
    
    # Load test data
    stg_file_path = Path("/app/backend/test_toy_14_dd.stg")
    ini_file_path = Path("/app/backend/test_toy_14_dd.ini")
    
    if not stg_file_path.exists() or not ini_file_path.exists():
        print("❌ Test files not found")
        return False
    
    try:
        # Read files
        with open(ini_file_path, 'r') as f:
            ini_content = f.read()
        with open(stg_file_path, 'r') as f:
            stg_content = f.read()
        
        print("✅ Test files loaded")
        
        # Initialize processor
        processor = EI2DRealDataProcessor()
        
        # Check if library is loaded
        if processor.lib is None:
            print("❌ EI2D library not loaded")
            return False
        
        print("✅ EI2D library loaded successfully")
        
        # Parse data
        ini_data = INIParser.parse_ini(ini_content)
        stg_data = STGParser.parse_stg(stg_content)
        
        electrodes = stg_data["electrodes"]
        measurements = stg_data["full_measurements"]
        
        print(f"✅ Parsed: {len(electrodes)} electrodes, {len(measurements)} measurements")
        
        # Test the real forward modeling function directly
        print("\n🎯 Testing _run_real_forward_modeling directly...")
        
        forward_params = ini_data.get("Forward", {})
        forw_mod_meth = int(forward_params.get("ForwModMeth", "1"))  # 0=FD, 1=FE
        forw_solver = int(forward_params.get("ForwSolver", "0"))
        bc_type = int(forward_params.get("BCType", "0"))
        forw_accuracy = int(forward_params.get("ForwAccuracy", "1"))
        forw_cg_iter = int(forward_params.get("ForwCGIter", "100"))
        forw_cg_resid = float(forward_params.get("ForwCGResid", "1E-6"))
        
        print(f"Forward parameters: method={forw_mod_meth}, solver={forw_solver}, bc={bc_type}")
        
        # Call the real forward modeling function
        result = processor._run_real_forward_modeling(
            electrodes, measurements, forw_mod_meth, forw_solver, 
            bc_type, forw_accuracy, forw_cg_iter, forw_cg_resid
        )
        
        print(f"\n📊 RESULT ANALYSIS:")
        print(f"Success: {result.get('success', False)}")
        print(f"Method: {result.get('method', 'unknown')}")
        print(f"Message: {result.get('message', 'no message')}")
        
        if result.get('success'):
            vi_data = result.get('results', {}).get('vi_data', [])
            if vi_data:
                finite_count = sum(1 for v in vi_data if isinstance(v, (int, float)) and abs(v) > 1e-12)
                print(f"VI data: {len(vi_data)} values, {finite_count} finite non-zero")
                if finite_count > 0:
                    print(f"Sample values: {vi_data[:5]}")
                    print(f"✅ REAL C-ABI WORKING: Finite VI values computed")
                else:
                    print(f"⚠️  VI values are zero - may indicate array bounds issue")
            
            # Check if real C-ABI was used
            method = result.get('method', '')
            if 'real_ei2d' in method:
                print(f"✅ REAL C-ABI CONFIRMED: Method = {method}")
            else:
                print(f"❌ FALLBACK USED: Method = {method}")
        else:
            print(f"❌ Forward modeling failed")
            error = result.get('error', 'no error details')
            print(f"Error: {error}")
        
        return result.get('success', False)
        
    except Exception as e:
        print(f"❌ Direct C-ABI test failed: {e}")
        import traceback
        print(f"Traceback: {traceback.format_exc()}")
        return False

if __name__ == "__main__":
    test_direct_cabi()