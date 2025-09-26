#!/usr/bin/env python3
"""
Debug the actual mesh dimensions vs parameter array values in the forward modeling
"""

import sys
sys.path.append('/app/backend')

from earthimager_wrapper import EI2DRealDataProcessor

def debug_mesh_parameters():
    """Debug the actual mesh parameters being set"""
    
    print("=== DEBUGGING ACTUAL MESH PARAMETERS ===\n")
    
    # Read the actual test files
    with open('/app/backend/test_toy_14_dd.ini', 'r') as f:
        ini_content = f.read()
    
    with open('/app/backend/test_toy_14_dd.stg', 'r') as f:
        stg_content = f.read()
    
    # Create processor and parse files
    processor = EI2DRealDataProcessor()
    
    # Parse the files to get the same data the forward modeling uses
    ini_result = processor.process_ini_stg_files(ini_content, stg_content)
    
    if not ini_result["success"]:
        print(f"Error parsing files: {ini_result.get('error')}")
        return
    
    # Extract the same parameters that forward modeling uses
    forward_params = ini_result["forward_parameters"] 
    stg_data = ini_result["stg_data"]
    
    print("=== INI FILE PARAMETERS ===")
    forw_mod_meth = int(forward_params.get("ForwModMeth", "1"))
    print(f"ForwModMeth: {forw_mod_meth} ({'FE' if forw_mod_meth == 1 else 'FD'})")
    
    print("\n=== STG DATA ===")
    electrodes = stg_data["electrodes"]
    measurements = stg_data["full_measurements"] 
    num_electrodes = len(electrodes)
    num_measurements = len(measurements)
    
    print(f"Electrodes: {num_electrodes}")
    print(f"Measurements: {num_measurements}")
    print(f"Electrode positions: {[e['position'] for e in electrodes[:5]]}... (first 5)")
    
    print("\n=== MESH DIMENSIONS CALCULATION ===")
    # This is the same calculation as in _run_real_forward_modeling
    
    # Conservative mesh setup (as used in current code)
    nElec = num_electrodes  # 14
    nNx = nElec + 2        # 16 nodes in X direction (conservative padding)
    nNy = 6                # 6 nodes in Y direction (depth layers)
    
    print(f"Node dimensions: nNx = {nNx}, nNy = {nNy}")
    print(f"Total nodes: {nNx * nNy}")
    
    # Element dimensions (as calculated by Fortran InitForwGlobals)
    gNumElemX = nNx - 1    # 15 elements in X
    gNumElemY = nNy - 1    # 5 elements in Y  
    gNumElem = gNumElemX * gNumElemY  # 75 elements
    
    print(f"Element dimensions: gNumElemX = {gNumElemX}, gNumElemY = {gNumElemY}")
    print(f"Total elements: gNumElem = {gNumElem}")
    
    print("\n=== PARAMETER ARRAY VALUES (CURRENT) ===")
    # Current parameter setup (as fixed)
    p1_current = 1
    p2_current = nNx - 1  # 15 (element bound)
    q1_current = 1  
    q2_current = nNy - 1  # 5 (element bound)
    
    print(f"Parameter bounds: p1={p1_current}, p2={p2_current}, q1={q1_current}, q2={q2_current}")
    print(f"ParamX range: [1, {p2_current}] (should be [1, {gNumElemX}])")
    print(f"ParamY range: [1, {q2_current}] (should be [1, {gNumElemY}])")
    
    if p2_current == gNumElemX and q2_current == gNumElemY:
        print("✅ Parameter bounds match element dimensions - should work")
    else:
        print("❌ Parameter bounds don't match element dimensions!")
        print(f"  p2={p2_current} vs gNumElemX={gNumElemX}")
        print(f"  q2={q2_current} vs gNumElemY={gNumElemY}")
    
    print("\n=== MAXIMUM iElem CALCULATION ===")
    # Calculate max iElem with current parameters
    max_iElem = q2_current + (p2_current - 1) * gNumElemY
    print(f"Max iElem = {q2_current} + ({p2_current}-1) * {gNumElemY} = {max_iElem}")
    print(f"Array size (gNumElem) = {gNumElem}")
    
    if max_iElem <= gNumElem:
        print("✅ Max iElem ≤ array size - bounds check passes")
    else:
        print(f"❌ Max iElem ({max_iElem}) > array size ({gNumElem}) - BOUNDS ERROR!")
    
    print("\n=== TESTING ERROR SCENARIO ===")
    # Test what would cause the "Index 61 above bound 60" error
    error_gNumElem = 60
    error_bad_index = 61
    
    print(f"To get gNumElem = {error_gNumElem}, we need different mesh dimensions.")
    print("Possible scenarios:")
    
    # Test some likely scenarios
    scenarios = [
        (11, 6),   # 11x6 nodes → 10x5 elements = 50 (too small)
        (12, 6),   # 12x6 nodes → 11x5 elements = 55 (too small) 
        (13, 6),   # 13x6 nodes → 12x5 elements = 60 (matches!)
        (14, 6),   # 14x6 nodes → 13x5 elements = 65 (too big)
        (16, 5),   # 16x5 nodes → 15x4 elements = 60 (matches!)
    ]
    
    for test_nNx, test_nNy in scenarios:
        test_gNumElemX = test_nNx - 1
        test_gNumElemY = test_nNy - 1
        test_gNumElem = test_gNumElemX * test_gNumElemY
        
        print(f"  {test_nNx}x{test_nNy} nodes → {test_gNumElemX}x{test_gNumElemY} elements = {test_gNumElem}")
        
        if test_gNumElem == error_gNumElem:
            print(f"    ✅ This gives gNumElem = {error_gNumElem}")
            
            # Test what parameter values would cause index 61
            if test_nNx - 1 + 1 > test_gNumElemX:  # If we set p2 = nNx instead of nNx-1
                print(f"    ❌ Setting p2 = {test_nNx} (instead of {test_nNx-1}) would cause bounds error")
            
            if test_nNy - 1 + 1 > test_gNumElemY:  # If we set q2 = nNy instead of nNy-1  
                print(f"    ❌ Setting q2 = {test_nNy} (instead of {test_nNy-1}) would cause bounds error")
    
    print(f"\n=== CONCLUSION ===")
    print("The issue likely occurs when:")
    print("1. Different mesh generation creates different nNx/nNy values")
    print("2. Parameter arrays are not correctly synchronized with actual mesh dimensions")
    print("3. There's a mismatch between what InitForwGlobals sets up vs what we pass")

if __name__ == "__main__":
    debug_mesh_parameters()