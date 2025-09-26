#!/usr/bin/env python3
"""
Test the array bounds fix by simulating the corrected parameter setup
"""

import numpy as np

def test_corrected_bounds():
    """Test that the corrected parameter bounds prevent array overflow"""
    
    print("=== TESTING CORRECTED ARRAY BOUNDS ===\n")
    
    # Test with problematic case from error: gNumElem = 539, trying to access 540
    # Let's use realistic mesh dimensions that could lead to this
    
    test_cases = [
        {"nNx": 25, "nNy": 22, "desc": "Medium mesh"},
        {"nNx": 16, "nNy": 6, "desc": "Conservative mesh (14 electrodes + padding)"},
        {"nNx": 28, "nNy": 20, "desc": "Aggressive mesh"},
        {"nNx": 50, "nNy": 12, "desc": "High resolution mesh (49x11 elements = 539)"},
    ]
    
    for i, case in enumerate(test_cases):
        nNx, nNy = case["nNx"], case["nNy"]
        desc = case["desc"]
        
        print(f"Test Case {i+1}: {desc}")
        print(f"  Nodes: {nNx}x{nNy} = {nNx * nNy}")
        
        # Calculate Fortran variables (as done in InitForwGlobals.f90)
        gNumNodeX = nNx
        gNumNodeY = nNy
        gNumElemX = gNumNodeX - 1
        gNumElemY = gNumNodeY - 1
        gNumElem = gNumElemX * gNumElemY
        
        print(f"  Elements: {gNumElemX}x{gNumElemY} = {gNumElem}")
        
        # OLD (buggy) parameter setup
        p2_old = nNx  # Wrong: using node dimension
        q2_old = nNy  # Wrong: using node dimension
        max_iElem_old = q2_old + (p2_old - 1) * gNumElemY
        
        print(f"  OLD bounds: iElemX ∈ [1,{p2_old}], iElemY ∈ [1,{q2_old}]")
        print(f"    Max iElem = {q2_old} + ({p2_old}-1) * {gNumElemY} = {max_iElem_old}")
        print(f"    Array size = {gNumElem}")
        
        if max_iElem_old > gNumElem:
            print(f"    ❌ BOUNDS ERROR: {max_iElem_old} > {gNumElem} (overflow by {max_iElem_old - gNumElem})")
        else:
            print(f"    ✅ OK: {max_iElem_old} <= {gNumElem}")
        
        # NEW (fixed) parameter setup  
        p2_new = gNumElemX  # Correct: using element dimension
        q2_new = gNumElemY  # Correct: using element dimension
        max_iElem_new = q2_new + (p2_new - 1) * gNumElemY
        
        print(f"  NEW bounds: iElemX ∈ [1,{p2_new}], iElemY ∈ [1,{q2_new}]")
        print(f"    Max iElem = {q2_new} + ({p2_new}-1) * {gNumElemY} = {max_iElem_new}")
        print(f"    Array size = {gNumElem}")
        
        if max_iElem_new > gNumElem:
            print(f"    ❌ STILL BROKEN: {max_iElem_new} > {gNumElem}")
        else:
            print(f"    ✅ FIXED: {max_iElem_new} <= {gNumElem}")
        
        print()
    
    # Demonstrate the specific error case
    print("=== SPECIFIC ERROR CASE ===")
    print("Error message: 'Index 540 of dimension 1 of array n1 above upper bound of 539'")
    
    # Find mesh that gives gNumElem = 539
    target_gNumElem = 539
    
    # Try different factorizations
    for gNumElemY in range(5, 25):
        if target_gNumElem % gNumElemY == 0:
            gNumElemX = target_gNumElem // gNumElemY
            nNx = gNumElemX + 1
            nNy = gNumElemY + 1
            
            print(f"\nFound: {nNx}x{nNy} nodes → {gNumElemX}x{gNumElemY} elements = {gNumElemX * gNumElemY}")
            
            # Test old vs new parameter setup
            p2_old = nNx
            q2_old = nNy 
            max_iElem_old = q2_old + (p2_old - 1) * gNumElemY
            
            p2_new = gNumElemX
            q2_new = gNumElemY
            max_iElem_new = q2_new + (p2_new - 1) * gNumElemY
            
            print(f"  OLD: max iElem = {q2_old} + ({p2_old}-1) * {gNumElemY} = {max_iElem_old}")
            print(f"  NEW: max iElem = {q2_new} + ({p2_new}-1) * {gNumElemY} = {max_iElem_new}")
            
            if max_iElem_old == 540 and target_gNumElem == 539:
                print(f"  🎯 THIS IS THE ERROR CASE!")
                print(f"     Array n1[539], trying to access n1[540]")
                print(f"     With fix: max access = n1[{max_iElem_new}] ✅")
                break

if __name__ == "__main__":
    test_corrected_bounds()