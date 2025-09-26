#!/usr/bin/env python3
"""
Debug forward modeling to identify why V/I values are zero
"""

import sys
sys.path.append('/app/backend')

import numpy as np
from earthimager_wrapper import EI2DRealDataProcessor

def debug_forward_modeling():
    """Debug the forward modeling step by step"""
    
    print("=== FORWARD MODELING DEBUG ===\n")
    
    # Read test files
    with open('/app/backend/test_toy_14_dd.ini', 'r') as f:
        ini_content = f.read()
    
    with open('/app/backend/test_toy_14_dd.stg', 'r') as f:
        stg_content = f.read()
    
    # Create processor
    processor = EI2DRealDataProcessor()
    
    print("1. TESTING FILE PROCESSING:")
    try:
        result = processor.process_ini_stg_files(ini_content, stg_content)
        
        print(f"   Success: {result.get('success')}")
        print(f"   Method: {result.get('method')}")
        
        if result.get('results'):
            vi_data = result['results'].get('vi_data', [])
            app_res = result['results'].get('apparent_resistivities', [])
            geom_factors = result['results'].get('geometric_factors', [])
            
            print(f"   V/I data: {len(vi_data)} values")
            print(f"     Sample values: {vi_data[:5] if vi_data else 'None'}")
            print(f"     Sum total: {sum(vi_data) if vi_data else 0}")
            
            print(f"   Apparent resistivities: {len(app_res)} values")
            print(f"     Sample values: {app_res[:5] if app_res else 'None'}")
            print(f"     Sum total: {sum(app_res) if app_res else 0}")
            
            print(f"   Geometric factors: {len(geom_factors)} values")
            print(f"     Sample values: {geom_factors[:5] if geom_factors else 'None'}")
            print(f"     Sum total: {sum(geom_factors) if geom_factors else 0}")
        
        if result.get('mesh'):
            mesh = result['mesh']
            print(f"   Mesh info:")
            print(f"     Conductivity elements: {len(mesh.get('conductivity', []))}")
            if mesh.get('conductivity'):
                cond_array = mesh['conductivity']
                print(f"     Conductivity range: {min(cond_array):.6f} to {max(cond_array):.6f} S/m")
                print(f"     Unique values: {len(set(cond_array))} (should be >1 for heterogeneous)")
        
        print(f"\n2. DETAILED ANALYSIS:")
        print(f"   Forward method from INI: {result.get('parameters', {}).get('forward_method')}")
        print(f"   Number of measurements: {len(result.get('stg_data', {}).get('measurements', []))}")
        print(f"   Number of electrodes: {len(result.get('stg_data', {}).get('electrodes', []))}")
        
    except Exception as e:
        print(f"   ERROR: {e}")
        import traceback
        traceback.print_exc()
    
    print(f"\n=== DIAGNOSIS ===")
    print("If V/I values are all zero:")
    print("1. Homogeneous conductivity → No voltage differences")
    print("2. Incorrect current injection setup")
    print("3. Wrong Fortran method (FD vs FE mismatch)")
    print("4. Missing boundary conditions")
    print("5. Array setup issues in Fortran code")

if __name__ == "__main__":
    debug_forward_modeling()