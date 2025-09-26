#!/usr/bin/env python3
"""
Debug the Sensitivity.f90 array bounds issue by analyzing the exact parameter values
"""

import numpy as np

def analyze_sensitivity_bounds():
    """Analyze the exact bounds issue in Sensitivity.f90"""
    
    print("=== SENSITIVITY.F90 ARRAY BOUNDS ANALYSIS ===\n")
    
    # From error message: "Index 61 above upper bound of 60"
    # This means: n1 array size = 60, trying to access n1(61)
    target_array_size = 60
    target_bad_index = 61
    
    print(f"Target: n1 array size = {target_array_size}, bad access = n1({target_bad_index})")
    
    # n1 array is allocated as: allocate(n1(gNumElem))
    gNumElem = target_array_size  # 60
    print(f"Therefore: gNumElem = {gNumElem}")
    
    # Find gNumElemX and gNumElemY such that gNumElemX * gNumElemY = 60
    print(f"\nPossible mesh dimensions for gNumElem = {gNumElem}:")
    factorizations = []
    for gNumElemY in range(1, int(gNumElem**0.5) + 1):
        if gNumElem % gNumElemY == 0:
            gNumElemX = gNumElem // gNumElemY
            factorizations.append((gNumElemX, gNumElemY))
            print(f"  gNumElemX = {gNumElemX}, gNumElemY = {gNumElemY}")
    
    print(f"\n=== ANALYZING INDEX CALCULATION ===")
    print("In Sensitivity.f90 lines 346-349:")
    print("  do iElemY = ParamY1(ipy), ParamY2(ipy)")
    print("     do iElemX = ParamX1(ipx), ParamX2(ipx)")
    print("        iElem = iElemY + (iElemX-1) * gNumElemY")
    
    # Test each factorization
    for gNumElemX, gNumElemY in factorizations:
        print(f"\n--- Testing mesh: {gNumElemX}x{gNumElemY} elements ---")
        
        # Valid element indices should be: iElemX ∈ [1, gNumElemX], iElemY ∈ [1, gNumElemY]
        print(f"Valid ranges: iElemX ∈ [1,{gNumElemX}], iElemY ∈ [1,{gNumElemY}]")
        
        # Calculate max valid iElem
        max_valid_iElem = gNumElemY + (gNumElemX - 1) * gNumElemY
        print(f"Max valid iElem = {gNumElemY} + ({gNumElemX}-1) * {gNumElemY} = {max_valid_iElem}")
        
        if max_valid_iElem == gNumElem:
            print("✅ This is correct - max iElem equals array size")
        else:
            print(f"❌ Mismatch: max iElem = {max_valid_iElem} vs array size = {gNumElem}")
        
        # Now find what ParamX2/ParamY2 values would cause iElem = 61
        print(f"\nFinding ParamX2/ParamY2 that cause iElem = {target_bad_index}:")
        
        for test_iElemY in range(1, gNumElemY + 5):  # Test beyond valid range
            for test_iElemX in range(1, gNumElemX + 5):  # Test beyond valid range
                test_iElem = test_iElemY + (test_iElemX - 1) * gNumElemY
                
                if test_iElem == target_bad_index:
                    print(f"  FOUND: iElemX={test_iElemX}, iElemY={test_iElemY} → iElem={test_iElem}")
                    
                    if test_iElemX > gNumElemX:
                        print(f"    ❌ PROBLEM: iElemX={test_iElemX} > valid range [1,{gNumElemX}]")
                        print(f"       This means ParamX2 = {test_iElemX} but should be ≤ {gNumElemX}")
                    
                    if test_iElemY > gNumElemY:
                        print(f"    ❌ PROBLEM: iElemY={test_iElemY} > valid range [1,{gNumElemY}]") 
                        print(f"       This means ParamY2 = {test_iElemY} but should be ≤ {gNumElemY}")
    
    print(f"\n=== CONCLUSION ===")
    print("The array bounds error occurs because:")
    print("1. Parameter arrays (ParamX2, ParamY2) contain values larger than mesh dimensions")
    print("2. This causes iElem calculation to exceed n1 array bounds")
    print("3. The fix is to ensure ParamX2 ≤ gNumElemX and ParamY2 ≤ gNumElemY")
    
    print(f"\nIn Python wrapper, we set:")
    print("  p2 = np.array([nNx - 1])  # Should be gNumElemX")  
    print("  q2 = np.array([nNy - 1])  # Should be gNumElemY")
    print("Where nNx-1 = gNumElemX and nNy-1 = gNumElemY")
    
    print(f"\nThe issue might be that our mesh generation produces different")
    print(f"dimensions than what we pass to the parameter arrays.")

if __name__ == "__main__":
    analyze_sensitivity_bounds()