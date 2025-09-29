#!/usr/bin/env python3
"""
Test JSON serialization fix for EarthImager 2D backend
Specifically tests the inversion workflow that was failing with "Out of range float values are not JSON compliant"
"""

import requests
import json
import os
from pathlib import Path

# Configuration
BACKEND_URL = "https://ei2d-integration.preview.emergentagent.com/api"
TEST_DATA_DIR = "/app/backend"

def test_json_serialization_fix():
    """Test that the JSON serialization fix prevents 'Out of range float values are not JSON compliant' errors"""
    
    print("🔍 Testing JSON serialization fix for inversion workflow...")
    
    stg_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.stg"
    ini_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.ini"
    
    if not stg_file_path.exists() or not ini_file_path.exists():
        print("❌ Test files not found")
        return False
    
    try:
        with open(ini_file_path, 'rb') as ini_f, open(stg_file_path, 'rb') as stg_f:
            files = {
                'ini_file': ('test_toy_14_dd.ini', ini_f, 'text/plain'),
                'stg_file': ('test_toy_14_dd.stg', stg_f, 'text/plain')
            }
            
            print("📡 Calling /api/earthimager/run-inversion...")
            response = requests.post(
                f"{BACKEND_URL}/earthimager/run-inversion",
                files=files,
                timeout=90
            )
        
        print(f"📊 Response status: {response.status_code}")
        
        if response.status_code == 200:
            try:
                data = response.json()
                print("✅ JSON parsing successful - no serialization errors!")
                
                # Verify key fields are present and finite
                success = data.get("success", False)
                workflow = data.get("workflow", "")
                parameters = data.get("parameters", {})
                results = data.get("results", {})
                
                print(f"   Success: {success}")
                print(f"   Workflow: {workflow}")
                print(f"   Final RMS: {parameters.get('final_rms', 'N/A')}")
                print(f"   Convergence: {parameters.get('convergence', 'N/A')}")
                
                # Check for NaN/infinity in critical arrays
                resistivity_model = results.get("resistivity_model", [])
                calculated_data = results.get("calculated_data", [])
                data_residuals = results.get("data_residuals", [])
                
                print(f"   Resistivity model length: {len(resistivity_model)}")
                print(f"   Calculated data length: {len(calculated_data)}")
                print(f"   Data residuals length: {len(data_residuals)}")
                
                # Verify no NaN/infinity values in arrays
                import math
                
                def has_invalid_values(arr):
                    return any(math.isnan(x) or math.isinf(x) for x in arr if isinstance(x, (int, float)))
                
                resistivity_invalid = has_invalid_values(resistivity_model)
                calculated_invalid = has_invalid_values(calculated_data)
                residuals_invalid = has_invalid_values(data_residuals)
                
                print(f"   Resistivity model has NaN/inf: {resistivity_invalid}")
                print(f"   Calculated data has NaN/inf: {calculated_invalid}")
                print(f"   Data residuals have NaN/inf: {residuals_invalid}")
                
                if not (resistivity_invalid or calculated_invalid or residuals_invalid):
                    print("✅ All float arrays are JSON-compliant (no NaN/infinity values)")
                    return True
                else:
                    print("❌ Some arrays still contain NaN/infinity values")
                    return False
                    
            except json.JSONDecodeError as e:
                print(f"❌ JSON parsing failed: {e}")
                print(f"   Response content (first 500 chars): {response.text[:500]}")
                return False
                
        else:
            print(f"❌ HTTP error: {response.status_code}")
            print(f"   Response: {response.text[:500]}")
            return False
            
    except Exception as e:
        print(f"❌ Request failed: {e}")
        return False

def test_forward_modeling_json():
    """Test forward modeling JSON serialization (should already be working)"""
    
    print("\n🔍 Testing forward modeling JSON serialization...")
    
    stg_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.stg"
    ini_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.ini"
    
    try:
        with open(ini_file_path, 'rb') as ini_f, open(stg_file_path, 'rb') as stg_f:
            files = {
                'ini_file': ('test_toy_14_dd.ini', ini_f, 'text/plain'),
                'stg_file': ('test_toy_14_dd.stg', stg_f, 'text/plain')
            }
            
            print("📡 Calling /api/earthimager/forward-model-real...")
            response = requests.post(
                f"{BACKEND_URL}/earthimager/forward-model-real",
                files=files,
                timeout=45
            )
        
        if response.status_code == 200:
            try:
                data = response.json()
                print("✅ Forward modeling JSON parsing successful!")
                
                # Check VI data for NaN/infinity
                vi_data = data.get("results", {}).get("vi_data", [])
                apparent_res = data.get("results", {}).get("apparent_resistivities", [])
                
                import math
                vi_invalid = any(math.isnan(x) or math.isinf(x) for x in vi_data if isinstance(x, (int, float)))
                ar_invalid = any(math.isnan(x) or math.isinf(x) for x in apparent_res if isinstance(x, (int, float)))
                
                print(f"   VI data has NaN/inf: {vi_invalid}")
                print(f"   Apparent resistivities have NaN/inf: {ar_invalid}")
                
                return not (vi_invalid or ar_invalid)
                
            except json.JSONDecodeError as e:
                print(f"❌ Forward modeling JSON parsing failed: {e}")
                return False
        else:
            print(f"❌ Forward modeling HTTP error: {response.status_code}")
            return False
            
    except Exception as e:
        print(f"❌ Forward modeling request failed: {e}")
        return False

def main():
    """Main test execution"""
    print("🧪 JSON Serialization Fix Testing Suite")
    print("=" * 50)
    
    # Test 1: Forward modeling (should already work)
    forward_ok = test_forward_modeling_json()
    
    # Test 2: Inversion workflow (the main fix target)
    inversion_ok = test_json_serialization_fix()
    
    print("\n" + "=" * 50)
    print("📋 TEST SUMMARY:")
    print(f"   Forward modeling JSON: {'✅ PASS' if forward_ok else '❌ FAIL'}")
    print(f"   Inversion workflow JSON: {'✅ PASS' if inversion_ok else '❌ FAIL'}")
    
    if forward_ok and inversion_ok:
        print("\n🎉 JSON serialization fix SUCCESSFUL!")
        print("   No 'Out of range float values are not JSON compliant' errors detected")
        return True
    else:
        print("\n💥 JSON serialization fix INCOMPLETE!")
        print("   Some endpoints still have NaN/infinity serialization issues")
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)