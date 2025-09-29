#!/usr/bin/env python3
"""
CRITICAL DATA TRANSFER DEBUGGING TEST
Focus: Test /api/earthimager/forward-model-real endpoint to identify where non-zero V/I values are lost
"""

import requests
import json
import numpy as np
from pathlib import Path

# Configuration
BACKEND_URL = "https://ei2d-integration.preview.emergentagent.com/api"
TEST_DATA_DIR = "/app/backend"

def test_forward_model_real_data_transfer():
    """
    CRITICAL TEST: Debug data transfer from backend computation to frontend display
    Focus on the /api/earthimager/forward-model-real endpoint
    """
    print("🔍 CRITICAL DATA TRANSFER DEBUGGING TEST")
    print("=" * 60)
    
    # Load test files
    stg_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.stg"
    ini_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.ini"
    
    if not stg_file_path.exists() or not ini_file_path.exists():
        print("❌ Test files not found!")
        return False
    
    print(f"✓ Test files found:")
    print(f"  INI: {ini_file_path}")
    print(f"  STG: {stg_file_path}")
    
    try:
        # Make API call to forward-model-real endpoint
        with open(ini_file_path, 'rb') as ini_f, open(stg_file_path, 'rb') as stg_f:
            files = {
                'ini_file': ('test_toy_14_dd.ini', ini_f, 'text/plain'),
                'stg_file': ('test_toy_14_dd.stg', stg_f, 'text/plain')
            }
            
            print(f"\n🚀 Calling /api/earthimager/forward-model-real...")
            response = requests.post(
                f"{BACKEND_URL}/earthimager/forward-model-real",
                files=files,
                timeout=60
            )
        
        print(f"📡 Response Status: {response.status_code}")
        
        if response.status_code != 200:
            print(f"❌ API call failed: {response.status_code}")
            print(f"Response: {response.text}")
            return False
        
        # Parse JSON response
        try:
            data = response.json()
        except json.JSONDecodeError as e:
            print(f"❌ JSON decode error: {e}")
            print(f"Raw response: {response.text[:500]}...")
            return False
        
        print(f"✓ JSON response parsed successfully")
        
        # CRITICAL ANALYSIS: Inspect the response structure
        print(f"\n🔍 RESPONSE STRUCTURE ANALYSIS:")
        print(f"  Success: {data.get('success', 'NOT_FOUND')}")
        print(f"  Method: {data.get('method', 'NOT_FOUND')}")
        print(f"  Message: {data.get('message', 'NOT_FOUND')}")
        
        # Check if results section exists
        results = data.get('results', {})
        if not results:
            print(f"❌ CRITICAL: No 'results' section in response!")
            print(f"Available keys: {list(data.keys())}")
            return False
        
        print(f"✓ Results section found")
        print(f"  Results keys: {list(results.keys())}")
        
        # CRITICAL: Check vi_data array
        vi_data = results.get('vi_data', [])
        if not vi_data:
            print(f"❌ CRITICAL: No 'vi_data' in results!")
            print(f"Results structure: {json.dumps(results, indent=2)[:500]}...")
            return False
        
        print(f"✓ vi_data array found")
        print(f"  vi_data length: {len(vi_data)}")
        print(f"  vi_data type: {type(vi_data)}")
        
        # CRITICAL ANALYSIS: Inspect vi_data values
        print(f"\n🔍 VI_DATA ARRAY ANALYSIS:")
        
        # Convert to numpy for analysis
        vi_array = np.array(vi_data)
        
        print(f"  Array shape: {vi_array.shape}")
        print(f"  Array dtype: {vi_array.dtype}")
        print(f"  Sum: {np.sum(vi_array)}")
        print(f"  Non-zero count: {np.count_nonzero(vi_array)}")
        print(f"  Zero count: {np.sum(vi_array == 0)}")
        print(f"  NaN count: {np.sum(np.isnan(vi_array))}")
        print(f"  Infinite count: {np.sum(np.isinf(vi_array))}")
        print(f"  Min/Max: {np.min(vi_array)} / {np.max(vi_array)}")
        
        # Show first 10 values
        print(f"  First 10 values: {vi_data[:10]}")
        print(f"  Last 10 values: {vi_data[-10:]}")
        
        # CRITICAL: Check if backend computed non-zero values but they got lost
        backend_message = data.get('message', '')
        computed_realistic = "Computed realistic V/I values" in backend_message
        
        print(f"\n🔍 BACKEND COMPUTATION ANALYSIS:")
        print(f"  Backend message: {backend_message}")
        print(f"  Contains 'Computed realistic V/I values': {computed_realistic}")
        
        # Check apparent resistivities too
        app_res = results.get('apparent_resistivities', [])
        if app_res:
            app_res_array = np.array(app_res)
            print(f"\n🔍 APPARENT RESISTIVITIES ANALYSIS:")
            print(f"  Length: {len(app_res)}")
            print(f"  Sum: {np.sum(app_res_array)}")
            print(f"  Non-zero count: {np.count_nonzero(app_res_array)}")
            print(f"  Min/Max: {np.min(app_res_array)} / {np.max(app_res_array)}")
            print(f"  First 5: {app_res[:5]}")
        
        # CRITICAL DIAGNOSIS
        print(f"\n🎯 CRITICAL DIAGNOSIS:")
        
        if computed_realistic and np.count_nonzero(vi_array) == 0:
            print(f"❌ CONFIRMED BUG: Backend computed realistic V/I values but vi_data contains only zeros!")
            print(f"   This confirms the data transfer/conversion bug.")
            print(f"   The VI_realistic values are being lost during JSON serialization.")
            
            # Check if there are any non-zero values in other fields
            geometric_factors = results.get('geometric_factors', [])
            if geometric_factors:
                gf_array = np.array(geometric_factors)
                print(f"   Geometric factors non-zero count: {np.count_nonzero(gf_array)}")
                print(f"   Geometric factors sample: {geometric_factors[:5]}")
            
            return False
            
        elif np.count_nonzero(vi_array) > 0:
            print(f"✅ SUCCESS: vi_data contains {np.count_nonzero(vi_array)} non-zero values!")
            print(f"   Data transfer working correctly.")
            print(f"   Sample non-zero values: {[v for v in vi_data[:20] if v != 0][:5]}")
            return True
            
        else:
            print(f"⚠️  Backend did not compute realistic V/I values (expected for some methods)")
            print(f"   This may be normal depending on the forward modeling method used.")
            return True
    
    except Exception as e:
        print(f"❌ Test failed with exception: {e}")
        import traceback
        print(f"Traceback: {traceback.format_exc()}")
        return False

def main():
    """Run the critical data transfer debugging test"""
    print("CRITICAL DATA TRANSFER DEBUGGING TASK")
    print("Testing /api/earthimager/forward-model-real endpoint")
    print("Focus: Identify where non-zero V/I values are lost between backend and frontend")
    print()
    
    success = test_forward_model_real_data_transfer()
    
    if success:
        print(f"\n🎉 Data transfer test completed successfully!")
    else:
        print(f"\n💥 Data transfer bug confirmed - V/I values lost in conversion!")
    
    return success

if __name__ == "__main__":
    main()