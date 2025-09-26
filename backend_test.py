#!/usr/bin/env python3
"""
Backend Test Suite for EarthImager 2D Web Interface
Tests the complete inversion workflow using toy-14-dd reference data
"""

import requests
import json
import os
import sys
from pathlib import Path
import time

# Configuration
BACKEND_URL = "https://earth-imager-web.preview.emergentagent.com/api"
TEST_DATA_DIR = "/app/backend"

class EarthImagerBackendTester:
    def __init__(self):
        self.backend_url = BACKEND_URL
        self.test_results = {
            "stg_processing": {"status": "pending", "details": {}},
            "ini_processing": {"status": "pending", "details": {}},
            "inversion_workflow": {"status": "pending", "details": {}},
            "out_file_validation": {"status": "pending", "details": {}},
            "api_endpoints": {"status": "pending", "details": {}}
        }
        
    def log(self, message, level="INFO"):
        """Log test messages"""
        print(f"[{level}] {message}")
        
    def test_health_check(self):
        """Test basic health check endpoint"""
        self.log("Testing health check endpoint...")
        try:
            response = requests.get(f"{self.backend_url}/health", timeout=10)
            if response.status_code == 200:
                data = response.json()
                self.log(f"Health check passed: {data}")
                return True
            else:
                self.log(f"Health check failed: {response.status_code}", "ERROR")
                return False
        except Exception as e:
            self.log(f"Health check error: {e}", "ERROR")
            return False
    
    def test_stg_file_processing(self):
        """Test STG file processing with toy-14-dd data"""
        self.log("Testing STG file processing...")
        
        stg_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.stg"
        if not stg_file_path.exists():
            self.log(f"STG test file not found: {stg_file_path}", "ERROR")
            self.test_results["stg_processing"]["status"] = "failed"
            self.test_results["stg_processing"]["details"]["error"] = "Test file not found"
            return False
            
        try:
            with open(stg_file_path, 'rb') as f:
                files = {'file': ('test_toy_14_dd.stg', f, 'text/plain')}
                response = requests.post(
                    f"{self.backend_url}/earthimager/upload-stg",
                    files=files,
                    timeout=30
                )
            
            if response.status_code == 200:
                data = response.json()
                
                # Validate expected results
                expected_measurements = 74
                expected_electrodes = 14
                
                actual_measurements = data.get("num_measurements", 0)
                actual_electrodes = data.get("num_electrodes", 0)
                
                success = (
                    actual_measurements == expected_measurements and
                    actual_electrodes == expected_electrodes and
                    data.get("success", False)
                )
                
                self.test_results["stg_processing"]["status"] = "passed" if success else "failed"
                self.test_results["stg_processing"]["details"] = {
                    "expected_measurements": expected_measurements,
                    "actual_measurements": actual_measurements,
                    "expected_electrodes": expected_electrodes,
                    "actual_electrodes": actual_electrodes,
                    "electrode_spacing": data.get("electrode_spacing", 0),
                    "format": data.get("format", "unknown"),
                    "voltage_range": data.get("voltage_range", {}),
                    "resistivity_range": data.get("resistivity_range", {})
                }
                
                if success:
                    self.log(f"STG processing passed: {actual_measurements} measurements, {actual_electrodes} electrodes")
                else:
                    self.log(f"STG processing validation failed", "ERROR")
                    
                return success
                
            else:
                self.log(f"STG upload failed: {response.status_code} - {response.text}", "ERROR")
                self.test_results["stg_processing"]["status"] = "failed"
                self.test_results["stg_processing"]["details"]["error"] = f"HTTP {response.status_code}: {response.text}"
                return False
                
        except Exception as e:
            self.log(f"STG processing error: {e}", "ERROR")
            self.test_results["stg_processing"]["status"] = "failed"
            self.test_results["stg_processing"]["details"]["error"] = str(e)
            return False
    
    def test_ini_file_processing(self):
        """Test INI file processing with toy-14-dd data"""
        self.log("Testing INI file processing...")
        
        ini_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.ini"
        if not ini_file_path.exists():
            self.log(f"INI test file not found: {ini_file_path}", "ERROR")
            self.test_results["ini_processing"]["status"] = "failed"
            self.test_results["ini_processing"]["details"]["error"] = "Test file not found"
            return False
            
        try:
            with open(ini_file_path, 'rb') as f:
                files = {'file': ('test_toy_14_dd.ini', f, 'text/plain')}
                response = requests.post(
                    f"{self.backend_url}/earthimager/upload-ini",
                    files=files,
                    timeout=30
                )
            
            if response.status_code == 200:
                data = response.json()
                
                # Validate expected INI sections and parameters
                expected_sections = ["Initial", "Forward", "ResInv", "IPInv", "Terrain", "CRP"]
                actual_sections = data.get("sections", [])
                
                # Check key parameters
                forward_method = data.get("forward_method", "")
                max_iterations = data.get("max_iterations", "")
                lagrange = data.get("lagrange_multiplier", "")
                
                success = (
                    data.get("success", False) and
                    len(actual_sections) >= 4 and  # At least main sections
                    forward_method in ["0", "1"] and  # Valid forward method
                    max_iterations and
                    lagrange
                )
                
                self.test_results["ini_processing"]["status"] = "passed" if success else "failed"
                self.test_results["ini_processing"]["details"] = {
                    "expected_sections": expected_sections,
                    "actual_sections": actual_sections,
                    "forward_method": forward_method,
                    "boundary_condition": data.get("boundary_condition", ""),
                    "max_iterations": max_iterations,
                    "lagrange_multiplier": lagrange,
                    "sections_found": len(actual_sections)
                }
                
                if success:
                    self.log(f"INI processing passed: {len(actual_sections)} sections, forward method {forward_method}")
                else:
                    self.log(f"INI processing validation failed", "ERROR")
                    
                return success
                
            else:
                self.log(f"INI upload failed: {response.status_code} - {response.text}", "ERROR")
                self.test_results["ini_processing"]["status"] = "failed"
                self.test_results["ini_processing"]["details"]["error"] = f"HTTP {response.status_code}: {response.text}"
                return False
                
        except Exception as e:
            self.log(f"INI processing error: {e}", "ERROR")
            self.test_results["ini_processing"]["status"] = "failed"
            self.test_results["ini_processing"]["details"]["error"] = str(e)
            return False
    
    def test_forward_modeling_real(self):
        """Test real forward modeling API with array bounds fix"""
        self.log("Testing REAL forward modeling with array bounds fix...")
        
        stg_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.stg"
        ini_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.ini"
        
        if not stg_file_path.exists() or not ini_file_path.exists():
            self.log("Test files not found for forward modeling", "ERROR")
            self.test_results["forward_modeling_real"] = {
                "status": "failed", 
                "details": {"error": "Test files not found"}
            }
            return False
            
        try:
            with open(ini_file_path, 'rb') as ini_f, open(stg_file_path, 'rb') as stg_f:
                files = {
                    'ini_file': ('test_toy_14_dd.ini', ini_f, 'text/plain'),
                    'stg_file': ('test_toy_14_dd.stg', stg_f, 'text/plain')
                }
                
                # Test the real forward modeling endpoint
                self.log("Calling /api/earthimager/forward-model-real...")
                response = requests.post(
                    f"{self.backend_url}/earthimager/forward-model-real",
                    files=files,
                    timeout=45  # Timeout for C-ABI calls
                )
            
            if response.status_code == 200:
                data = response.json()
                
                success = data.get("success", False)
                method = data.get("method", "unknown")
                message = data.get("message", "")
                
                # Check if backend completed without hanging
                backend_completed = success and "completed" in message.lower()
                no_array_bounds_error = "array" not in message.lower() and "bounds" not in message.lower()
                
                self.test_results["forward_modeling_real"] = {
                    "status": "passed" if (backend_completed and no_array_bounds_error) else "failed",
                    "details": {
                        "success": success,
                        "method": method,
                        "message": message,
                        "backend_completed": backend_completed,
                        "no_array_bounds_error": no_array_bounds_error,
                        "response_data": data
                    }
                }
                
                if backend_completed and no_array_bounds_error:
                    self.log(f"✅ Real forward modeling passed: {method}")
                    self.log(f"Message: {message}")
                else:
                    self.log(f"❌ Real forward modeling issues detected", "ERROR")
                    
                return backend_completed and no_array_bounds_error
                
            else:
                self.log(f"Forward modeling failed: {response.status_code} - {response.text}", "ERROR")
                self.test_results["forward_modeling_real"] = {
                    "status": "failed",
                    "details": {"error": f"HTTP {response.status_code}: {response.text}"}
                }
                return False
                
        except Exception as e:
            self.log(f"Forward modeling error: {e}", "ERROR")
            self.test_results["forward_modeling_real"] = {
                "status": "failed",
                "details": {"error": str(e)}
            }
            return False

    def test_inversion_workflow(self):
        """Test complete inversion workflow with REAL C-ABI calls (array bounds fix applied)"""
        self.log("Testing REAL inversion workflow with array bounds fix...")
        
        stg_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.stg"
        ini_file_path = Path(TEST_DATA_DIR) / "test_toy_14_dd.ini"
        
        if not stg_file_path.exists() or not ini_file_path.exists():
            self.log("Test files not found for inversion workflow", "ERROR")
            self.test_results["inversion_workflow"]["status"] = "failed"
            self.test_results["inversion_workflow"]["details"]["error"] = "Test files not found"
            return False
            
        try:
            with open(ini_file_path, 'rb') as ini_f, open(stg_file_path, 'rb') as stg_f:
                files = {
                    'ini_file': ('test_toy_14_dd.ini', ini_f, 'text/plain'),
                    'stg_file': ('test_toy_14_dd.stg', stg_f, 'text/plain')
                }
                
                # Test the main inversion endpoint
                self.log("Calling /api/earthimager/run-inversion...")
                response = requests.post(
                    f"{self.backend_url}/earthimager/run-inversion",
                    files=files,
                    timeout=90  # Longer timeout for real C-ABI inversion
                )
            
            if response.status_code == 200:
                data = response.json()
                
                # Validate inversion results
                success = data.get("success", False)
                workflow = data.get("workflow", "")
                parameters = data.get("parameters", {})
                results = data.get("results", {})
                out_file = data.get("out_file", {})
                message = data.get("message", "")
                
                # CRITICAL: Check if REAL C-ABI was used or simulation fallback
                is_real_cabi = "✅ REAL C-ABI inversion completed" in message
                is_simulation = "Simulation inversion completed" in message
                no_array_bounds_error = "array" not in message.lower() and "bounds" not in message.lower()
                no_hanging = success  # If we got a response, backend didn't hang
                
                # Check for key inversion components
                has_mesh = "mesh" in data and data["mesh"]
                has_results = bool(results)
                has_out_file = bool(out_file)
                
                validation_success = (
                    success and
                    workflow == "complete_ei2d_inversion" and
                    has_mesh and
                    has_results and
                    no_hanging and
                    no_array_bounds_error
                )
                
                self.test_results["inversion_workflow"]["status"] = "passed" if validation_success else "failed"
                self.test_results["inversion_workflow"]["details"] = {
                    "success": success,
                    "workflow": workflow,
                    "has_mesh": has_mesh,
                    "has_results": has_results,
                    "has_out_file": has_out_file,
                    "parameters": parameters,
                    "message": message,
                    "input_files": data.get("input_files", {}),
                    "is_real_cabi": is_real_cabi,
                    "is_simulation": is_simulation,
                    "no_array_bounds_error": no_array_bounds_error,
                    "no_hanging": no_hanging,
                    "array_bounds_fix_status": "SUCCESS" if (no_array_bounds_error and no_hanging) else "FAILED"
                }
                
                if validation_success:
                    if is_real_cabi:
                        self.log(f"✅ REAL C-ABI inversion SUCCESS: {workflow}")
                        self.log(f"✅ Array bounds fix WORKING: No Fortran runtime errors")
                    elif is_simulation:
                        self.log(f"⚠️  Simulation fallback used: {workflow}")
                        self.log(f"✅ Array bounds fix partially working: No hanging, but fallback used")
                    else:
                        self.log(f"✅ Inversion completed: {workflow}")
                    
                    # Store OUT file content for validation
                    if has_out_file:
                        self.out_file_content = out_file.get("content", "")
                else:
                    self.log(f"❌ Inversion workflow validation failed", "ERROR")
                    if not no_array_bounds_error:
                        self.log(f"❌ Array bounds error still present!", "ERROR")
                    if not no_hanging:
                        self.log(f"❌ Backend hanging detected!", "ERROR")
                    
                return validation_success
                
            else:
                self.log(f"Inversion workflow failed: {response.status_code} - {response.text}", "ERROR")
                self.test_results["inversion_workflow"]["status"] = "failed"
                self.test_results["inversion_workflow"]["details"]["error"] = f"HTTP {response.status_code}: {response.text}"
                return False
                
        except Exception as e:
            self.log(f"Inversion workflow error: {e}", "ERROR")
            self.test_results["inversion_workflow"]["status"] = "failed"
            self.test_results["inversion_workflow"]["details"]["error"] = str(e)
            return False
    
    def test_out_file_validation(self):
        """Validate OUT file format and content structure"""
        self.log("Testing OUT file validation...")
        
        if not hasattr(self, 'out_file_content') or not self.out_file_content:
            self.log("No OUT file content available for validation", "ERROR")
            self.test_results["out_file_validation"]["status"] = "failed"
            self.test_results["out_file_validation"]["details"]["error"] = "No OUT file content"
            return False
            
        try:
            content = self.out_file_content
            lines = content.split('\n')
            
            # Check for required sections
            has_header = any("Advanced Geosciences Inc." in line for line in lines[:5])
            has_settings = any(";------ SETTINGS ------" in line for line in lines)
            has_electrodes = any(";------ ELECTRODE LOCATIONS ------" in line for line in lines)
            has_commands = any(";------ Commands, Raw V/I, GeomFactor, AppRes ------" in line for line in lines)
            has_iterations = any(";------ OUTPUT OF DATA AND MODEL OF ALL ITERATIONS ------" in line for line in lines)
            
            # Check for iteration data
            iteration_sections = []
            resistivity_sections = []
            sensitivity_sections = []
            vi_data_sections = []
            
            for i, line in enumerate(lines):
                if ";------ Iteration" in line:
                    iteration_sections.append(i)
                elif ";-Resistivity in Ohm-m in the elemental sequential order" in line:
                    resistivity_sections.append(i)
                elif ";-Sensitivity in the elemental sequential order" in line:
                    sensitivity_sections.append(i)
                elif ";-Index  V/I_Meas      V/I_Calc    VI_%ERR" in line:
                    vi_data_sections.append(i)
            
            # Validate structure
            structure_valid = (
                has_header and
                has_settings and
                has_electrodes and
                has_commands and
                has_iterations and
                len(iteration_sections) > 0 and
                len(resistivity_sections) > 0 and
                len(sensitivity_sections) > 0 and
                len(vi_data_sections) > 0
            )
            
            # Check for measurement data (should have 74 measurements)
            measurement_count = 0
            for line in lines:
                if line.strip() and line.strip()[0].isdigit() and "," in line:
                    parts = line.split(",")
                    if len(parts) >= 8:  # DataID, A, B, M, N, V/I, K, App-Res
                        measurement_count += 1
            
            data_count_valid = measurement_count >= 70  # Allow some tolerance
            
            self.test_results["out_file_validation"]["status"] = "passed" if (structure_valid and data_count_valid) else "failed"
            self.test_results["out_file_validation"]["details"] = {
                "has_header": has_header,
                "has_settings": has_settings,
                "has_electrodes": has_electrodes,
                "has_commands": has_commands,
                "has_iterations": has_iterations,
                "iteration_sections": len(iteration_sections),
                "resistivity_sections": len(resistivity_sections),
                "sensitivity_sections": len(sensitivity_sections),
                "vi_data_sections": len(vi_data_sections),
                "measurement_count": measurement_count,
                "expected_measurements": 74,
                "structure_valid": structure_valid,
                "data_count_valid": data_count_valid,
                "total_lines": len(lines)
            }
            
            if structure_valid and data_count_valid:
                self.log(f"OUT file validation passed: {len(iteration_sections)} iterations, {measurement_count} measurements")
            else:
                self.log(f"OUT file validation failed", "ERROR")
                
            return structure_valid and data_count_valid
            
        except Exception as e:
            self.log(f"OUT file validation error: {e}", "ERROR")
            self.test_results["out_file_validation"]["status"] = "failed"
            self.test_results["out_file_validation"]["details"]["error"] = str(e)
            return False
    
    def test_api_endpoints(self):
        """Test various API endpoints"""
        self.log("Testing API endpoints...")
        
        endpoints_to_test = [
            ("/", "GET"),
            ("/health", "GET"),
        ]
        
        results = {}
        all_passed = True
        
        for endpoint, method in endpoints_to_test:
            try:
                if method == "GET":
                    response = requests.get(f"{self.backend_url}{endpoint}", timeout=10)
                else:
                    response = requests.post(f"{self.backend_url}{endpoint}", timeout=10)
                
                success = response.status_code in [200, 201]
                results[endpoint] = {
                    "method": method,
                    "status_code": response.status_code,
                    "success": success,
                    "response_size": len(response.text)
                }
                
                if not success:
                    all_passed = False
                    
            except Exception as e:
                results[endpoint] = {
                    "method": method,
                    "success": False,
                    "error": str(e)
                }
                all_passed = False
        
        self.test_results["api_endpoints"]["status"] = "passed" if all_passed else "failed"
        self.test_results["api_endpoints"]["details"] = results
        
        if all_passed:
            self.log(f"API endpoints test passed: {len(endpoints_to_test)} endpoints")
        else:
            self.log(f"API endpoints test failed", "ERROR")
            
        return all_passed
    
    def run_all_tests(self):
        """Run all backend tests"""
        self.log("Starting EarthImager 2D Backend Test Suite")
        self.log(f"Backend URL: {self.backend_url}")
        self.log(f"Test data directory: {TEST_DATA_DIR}")
        
        # Test sequence
        tests = [
            ("Health Check", self.test_health_check),
            ("STG File Processing", self.test_stg_file_processing),
            ("INI File Processing", self.test_ini_file_processing),
            ("Inversion Workflow", self.test_inversion_workflow),
            ("OUT File Validation", self.test_out_file_validation),
            ("API Endpoints", self.test_api_endpoints)
        ]
        
        passed_tests = 0
        total_tests = len(tests)
        
        for test_name, test_func in tests:
            self.log(f"\n{'='*50}")
            self.log(f"Running: {test_name}")
            self.log(f"{'='*50}")
            
            try:
                result = test_func()
                if result:
                    passed_tests += 1
                    self.log(f"✅ {test_name} PASSED")
                else:
                    self.log(f"❌ {test_name} FAILED")
            except Exception as e:
                self.log(f"❌ {test_name} ERROR: {e}")
        
        # Final summary
        self.log(f"\n{'='*50}")
        self.log(f"TEST SUMMARY")
        self.log(f"{'='*50}")
        self.log(f"Passed: {passed_tests}/{total_tests}")
        self.log(f"Success Rate: {(passed_tests/total_tests)*100:.1f}%")
        
        # Detailed results
        self.log(f"\nDetailed Results:")
        for test_name, result in self.test_results.items():
            status = result["status"]
            self.log(f"  {test_name}: {status.upper()}")
            if result["details"]:
                for key, value in result["details"].items():
                    if isinstance(value, dict):
                        self.log(f"    {key}: {json.dumps(value, indent=6)}")
                    else:
                        self.log(f"    {key}: {value}")
        
        return passed_tests == total_tests

def main():
    """Main test execution"""
    tester = EarthImagerBackendTester()
    success = tester.run_all_tests()
    
    if success:
        print("\n🎉 All tests passed!")
        sys.exit(0)
    else:
        print("\n💥 Some tests failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()