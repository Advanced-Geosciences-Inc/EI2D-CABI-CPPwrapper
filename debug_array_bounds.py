#!/usr/bin/env python3
"""
Toy test to debug the Fortran array bounds issue in EarthImager 2D
This will help us understand and fix the array indexing problem
"""

def calculate_array_bounds():
    """Calculate and verify array bounds for EarthImager 2D mesh"""
    
    print("=== EarthImager 2D Array Bounds Debug Test ===\n")
    
    # Test case: 14 electrodes (from toy-14-dd)
    nElec = 14
    nNx = nElec    # 14 nodes in X
    nNy = 6        # 6 nodes in Y (depth layers)
    
    print(f"Input dimensions:")
    print(f"  nElec = {nElec}")
    print(f"  nNx (NumNodeX) = {nNx}")  
    print(f"  nNy (NumNodeY) = {nNy}")
    
    # Python calculation (what we currently do)
    nElem_python = (nNx - 1) * (nNy - 1)
    print(f"\nPython calculation:")
    print(f"  nElem = (nNx-1) * (nNy-1) = ({nNx}-1) * ({nNy}-1) = {nElem_python}")
    
    # Fortran calculation (what InitForwGlobals.f90 does)
    gNumNodeX = nNx
    gNumNodeY = nNy
    gNumElemX = gNumNodeX - 1  # 13
    gNumElemY = gNumNodeY - 1  # 5
    gNumElem = gNumElemX * gNumElemY  # 65
    
    print(f"\nFortran calculation (InitForwGlobals.f90):")
    print(f"  gNumNodeX = {gNumNodeX}")
    print(f"  gNumNodeY = {gNumNodeY}")
    print(f"  gNumElemX = gNumNodeX - 1 = {gNumElemX}")
    print(f"  gNumElemY = gNumNodeY - 1 = {gNumElemY}")
    print(f"  gNumElem = gNumElemX * gNumElemY = {gNumElem}")
    
    # Array allocation in Sensitivity.f90 line 288
    print(f"\nArray allocation:")
    print(f"  n1 array allocated with size: gNumElem = {gNumElem}")
    print(f"  Valid indices: 1 to {gNumElem} (Fortran 1-based)")
    
    # Index calculation in Sensitivity.f90 lines 346-349
    print(f"\nIndex calculation in loops:")
    print(f"  iElemX ranges: 1 to gNumElemX = 1 to {gNumElemX}")
    print(f"  iElemY ranges: 1 to gNumElemY = 1 to {gNumElemY}")
    print(f"  iElem = iElemY + (iElemX-1) * gNumElemY")
    
    # Calculate min and max iElem
    iElem_min = 1 + (1-1) * gNumElemY  # iElemY=1, iElemX=1
    iElem_max = gNumElemY + (gNumElemX-1) * gNumElemY  # iElemY=gNumElemY, iElemX=gNumElemX
    
    print(f"\n  Minimum iElem = 1 + (1-1) * {gNumElemY} = {iElem_min}")
    print(f"  Maximum iElem = {gNumElemY} + ({gNumElemX}-1) * {gNumElemY} = {iElem_max}")
    
    # Check bounds
    print(f"\nBounds check:")
    print(f"  Array size: {gNumElem}")
    print(f"  Max index: {iElem_max}")
    
    if iElem_max > gNumElem:
        print(f"  ❌ BOUNDS ERROR: Max index {iElem_max} > Array size {gNumElem}")
        print(f"     Overflow by: {iElem_max - gNumElem}")
    elif iElem_max == gNumElem:
        print(f"  ⚠️  EDGE CASE: Max index {iElem_max} == Array size {gNumElem}")
        print(f"     This could be valid in Fortran 1-based indexing")
    else:
        print(f"  ✅ OK: Max index {iElem_max} < Array size {gNumElem}")
    
    # Analysis
    print(f"\n=== ANALYSIS ===")
    print(f"The issue occurs because:")
    print(f"1. We have {gNumElemX} elements in X direction")
    print(f"2. We have {gNumElemY} elements in Y direction")  
    print(f"3. Total elements = {gNumElem}")
    print(f"4. But indexing formula assumes column-major storage:")
    print(f"   iElem = iElemY + (iElemX-1) * gNumElemY")
    print(f"5. This gives valid range 1 to {iElem_max}")
    print(f"6. Array n1 is allocated for indices 1 to {gNumElem}")
    
    if iElem_max > gNumElem:
        print(f"\n7. FIX NEEDED: The indexing is consistent, but we need to check")
        print(f"   if the loop bounds are correct in the calling code.")
    
    return {
        'gNumElemX': gNumElemX,
        'gNumElemY': gNumElemY, 
        'gNumElem': gNumElem,
        'iElem_max': iElem_max,
        'bounds_ok': iElem_max <= gNumElem
    }

if __name__ == "__main__":
    result = calculate_array_bounds()
    print(f"\n=== RESULT ===")
    print(f"Bounds check: {'✅ PASSED' if result['bounds_ok'] else '❌ FAILED'}")