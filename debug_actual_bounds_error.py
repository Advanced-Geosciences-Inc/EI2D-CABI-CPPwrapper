#!/usr/bin/env python3
"""
Test to reproduce the exact array bounds error: Index 540 > upper bound 539
This will help us identify the actual mesh dimensions that cause the problem
"""

def reverse_engineer_error_case():
    """Find mesh dimensions that would cause n1 array size = 539 and iElem = 540"""
    
    print("=== Reverse Engineering Array Bounds Error ===\n")
    print("Error: Index 540 of dimension 1 of array n1 above upper bound of 539")
    print("This means: n1 array size = 539, attempting to access n1(540)")
    print()
    
    # Target: gNumElem = 539, max iElem = 540
    target_gNumElem = 539
    target_max_iElem = 540
    
    print(f"Target: gNumElem = {target_gNumElem}, max iElem = {target_max_iElem}")
    print()
    
    # Find gNumElemX and gNumElemY such that:
    # gNumElem = gNumElemX * gNumElemY = 539
    # max iElem = gNumElemY + (gNumElemX-1) * gNumElemY = gNumElemX * gNumElemY = 540
    
    # From the second equation: gNumElemX * gNumElemY = 540
    # From the first equation: gNumElemX * gNumElemY = 539
    # This is inconsistent! There's a mismatch.
    
    print("Checking if indexing formula is correct...")
    print("iElem = iElemY + (iElemX-1) * gNumElemY")
    print("Max iElem occurs when iElemY = gNumElemY, iElemX = gNumElemX")
    print("Max iElem = gNumElemY + (gNumElemX-1) * gNumElemY = gNumElemY * gNumElemX")
    print()
    
    if target_max_iElem != target_gNumElem:
        print(f"⚠️  INCONSISTENCY DETECTED:")
        print(f"   Array size: {target_gNumElem}")
        print(f"   Max index:  {target_max_iElem}")
        print(f"   Difference: {target_max_iElem - target_gNumElem}")
        print()
    
    # Let's find possible factorizations
    print("Finding possible mesh dimensions for gNumElem = 539:")
    factors_539 = []
    for gNumElemY in range(1, int(target_gNumElem**0.5) + 1):
        if target_gNumElem % gNumElemY == 0:
            gNumElemX = target_gNumElem // gNumElemY
            factors_539.append((gNumElemX, gNumElemY))
    
    print("Possible (gNumElemX, gNumElemY) for gNumElem = 539:")
    for gNumElemX, gNumElemY in factors_539:
        max_iElem = gNumElemY * gNumElemX
        print(f"  gNumElemX={gNumElemX}, gNumElemY={gNumElemY} → max iElem = {max_iElem}")
    
    print()
    print("Finding possible mesh dimensions for max iElem = 540:")
    factors_540 = []
    for gNumElemY in range(1, int(target_max_iElem**0.5) + 1):
        if target_max_iElem % gNumElemY == 0:
            gNumElemX = target_max_iElem // gNumElemY
            factors_540.append((gNumElemX, gNumElemY))
    
    print("Possible (gNumElemX, gNumElemY) for max iElem = 540:")
    for gNumElemX, gNumElemY in factors_540:
        gNumElem = gNumElemX * gNumElemY
        print(f"  gNumElemX={gNumElemX}, gNumElemY={gNumElemY} → gNumElem = {gNumElem}")
    
    print()
    print("=== ANALYSIS ===")
    print("The error suggests one of these scenarios:")
    print("1. Array allocated with size 539, but indexing assumes size 540")
    print("2. There's an off-by-one error in either allocation or indexing")
    print("3. The mesh dimensions are calculated incorrectly somewhere")
    
    # Let's check what 74 measurements with 14 electrodes might produce
    print()
    print("=== CHECKING ACTUAL STG DATA ===")
    print("toy-14-dd.stg has 74 measurements with 14 electrodes")
    
    # If we have 14 electrodes, what are likely mesh dimensions?
    nElec = 14
    print(f"With {nElec} electrodes:")
    
    # Conservative: surface nodes = electrodes, + padding
    nNx_conservative = nElec + 2  # 16 nodes in X
    nNy_conservative = 6          # 6 nodes in Y (depth)
    gNumElemX_conservative = nNx_conservative - 1  # 15 elements
    gNumElemY_conservative = nNy_conservative - 1  # 5 elements
    gNumElem_conservative = gNumElemX_conservative * gNumElemY_conservative  # 75
    
    print(f"  Conservative mesh: {nNx_conservative}x{nNy_conservative} nodes → {gNumElemX_conservative}x{gNumElemY_conservative} elements = {gNumElem_conservative}")
    
    # Aggressive: more elements for better resolution
    nNx_aggressive = nElec * 2     # 28 nodes in X
    nNy_aggressive = 20            # 20 nodes in Y (deeper)
    gNumElemX_aggressive = nNx_aggressive - 1  # 27 elements
    gNumElemY_aggressive = nNy_aggressive - 1  # 19 elements  
    gNumElem_aggressive = gNumElemX_aggressive * gNumElemY_aggressive  # 513
    
    print(f"  Aggressive mesh: {nNx_aggressive}x{nNy_aggressive} nodes → {gNumElemX_aggressive}x{gNumElemY_aggressive} elements = {gNumElem_aggressive}")
    
    # Very aggressive: even finer mesh  
    nNx_very = nElec * 3           # 42 nodes in X
    nNy_very = 15                  # 15 nodes in Y
    gNumElemX_very = nNx_very - 1  # 41 elements
    gNumElemY_very = nNy_very - 1  # 14 elements
    gNumElem_very = gNumElemX_very * gNumElemY_very  # 574
    
    print(f"  Very aggressive mesh: {nNx_very}x{nNy_very} nodes → {gNumElemX_very}x{gNumElemY_very} elements = {gNumElem_very}")
    
    if gNumElem_aggressive <= target_gNumElem <= gNumElem_very:
        print(f"\n✅ Target gNumElem {target_gNumElem} falls in the aggressive range!")
        print("This suggests the error occurs with a moderately fine mesh.")
    
    # Check which factorization is closest to realistic mesh dimensions
    print(f"\n=== LIKELY ERROR SCENARIO ===")
    for gNumElemX, gNumElemY in factors_539:
        nNx = gNumElemX + 1
        nNy = gNumElemY + 1
        if 10 <= nNx <= 50 and 5 <= nNy <= 25:  # Reasonable mesh dimensions
            print(f"Likely scenario: {nNx}x{nNy} nodes → {gNumElemX}x{gNumElemY} elements")
            print(f"  This could happen with {nNx-2} electrodes + padding")
            
            # Check if max index exceeds array size
            max_iElem = gNumElemY + (gNumElemX-1) * gNumElemY
            if max_iElem > target_gNumElem:
                print(f"  ❌ PROBLEM: max iElem = {max_iElem} > array size {target_gNumElem}")
            else:
                print(f"  ✅ OK: max iElem = {max_iElem} <= array size {target_gNumElem}")

if __name__ == "__main__":
    reverse_engineer_error_case()