#!/usr/bin/env python3
"""
DEBUG PARAMETER BOUNDS - Find exact cause of array bounds error
"""

import ctypes
import numpy as np
import os

def analyze_parameter_setup():
    """Analyze the parameter bounds setup for minimal test case"""
    
    print("ANALYZING PARAMETER BOUNDS FOR MINIMAL TEST CASE")
    print("="*60)
    
    # Test case parameters 
    nNx = 4  # 4 nodes in X
    nNy = 3  # 3 nodes in Y
    nNodes = nNx * nNy  # 12 nodes
    
    gNumElemX = nNx - 1  # 3 elements in X
    gNumElemY = nNy - 1  # 2 elements in Y  
    nElem = gNumElemX * gNumElemY  # 6 elements total
    
    print(f"Grid setup:")
    print(f"  Nodes: {nNx} x {nNy} = {nNodes}")
    print(f"  Elements: {gNumElemX} x {gNumElemY} = {nElem}")
    print(f"  gNumNodeY = {nNy} (used in n1 calculation)")
    
    # My current parameter bounds
    nParamX = 1
    nParamY = 1
    ParamX1 = [1]
    ParamX2 = [3]  # This might be wrong!
    ParamY1 = [1] 
    ParamY2 = [2]  # This might be wrong!
    
    print(f"\nCurrent parameter bounds:")
    print(f"  nParamX={nParamX}, nParamY={nParamY}")
    print(f"  ParamX1={ParamX1}, ParamX2={ParamX2}")
    print(f"  ParamY1={ParamY1}, ParamY2={ParamY2}")
    
    print(f"\nTracing iElem calculation in Jacobian loop:")
    print(f"  Formula: iElem = iElemY + (iElemX-1) * gNumElemY")
    
    max_iElem = 0
    for ipx in range(nParamX):
        for ipy in range(nParamY):
            print(f"\n  Parameter region [{ipx}][{ipy}]:")
            for iElemX in range(ParamX1[ipx], ParamX2[ipx] + 1):
                for iElemY in range(ParamY1[ipy], ParamY2[ipy] + 1):
                    iElem = iElemY + (iElemX - 1) * gNumElemY
                    print(f"    iElemX={iElemX}, iElemY={iElemY} → iElem={iElem}")
                    max_iElem = max(max_iElem, iElem)
    
    print(f"\n  Maximum iElem generated: {max_iElem}")
    print(f"  n1 array size (gNumElem): {nElem}")
    
    if max_iElem > nElem:
        print(f"  ❌ PROBLEM: iElem={max_iElem} exceeds n1 array bounds [1,{nElem}]")
        print(f"     This will cause 'Index {max_iElem} above upper bound of {nElem}'")
    else:
        print(f"  ✓ OK: All iElem values within bounds")
    
    # Analyze the root cause
    print(f"\nROOT CAUSE ANALYSIS:")
    print(f"  The parameter bounds define ELEMENT ranges, not node ranges")
    print(f"  Valid element indices should be 1 to gNumElem = {nElem}")
    print(f"  But element indexing in 2D grid uses: iElem = iElemY + (iElemX-1) * gNumElemY")
    
    print(f"\nCORRECT PARAMETER BOUNDS should be:")
    print(f"  For elements arranged as:")
    print(f"    X-direction: 1 to {gNumElemX}")  
    print(f"    Y-direction: 1 to {gNumElemY}")
    print(f"  So ParamX2 should be ≤ {gNumElemX}, ParamY2 should be ≤ {gNumElemY}")
    
    # Test corrected bounds
    print(f"\nTESTING CORRECTED BOUNDS:")
    ParamX2_corrected = [gNumElemX]  # Should be 3
    ParamY2_corrected = [gNumElemY]  # Should be 2
    
    print(f"  Corrected: ParamX1={ParamX1}, ParamX2={ParamX2_corrected}")
    print(f"  Corrected: ParamY1={ParamY1}, ParamY2={ParamY2_corrected}")
    
    max_iElem_corrected = 0
    for ipx in range(nParamX):
        for ipy in range(nParamY):
            for iElemX in range(ParamX1[ipx], ParamX2_corrected[ipx] + 1):
                for iElemY in range(ParamY1[ipy], ParamY2_corrected[ipy] + 1):
                    iElem = iElemY + (iElemX - 1) * gNumElemY
                    max_iElem_corrected = max(max_iElem_corrected, iElem)
    
    print(f"  Maximum iElem with corrected bounds: {max_iElem_corrected}")
    if max_iElem_corrected <= nElem:
        print(f"  ✓ CORRECTED bounds should work!")
    else:
        print(f"  ❌ Still problematic...")

def main():
    analyze_parameter_setup()

if __name__ == "__main__":
    main()