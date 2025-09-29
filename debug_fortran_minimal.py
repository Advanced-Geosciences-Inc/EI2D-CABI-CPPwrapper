#!/usr/bin/env python3
"""
MINIMAL FORTRAN DEBUGGING TEST
Create the smallest possible test case to isolate Fortran issues
"""

import ctypes
import numpy as np
import os
from pathlib import Path

class MinimalEI2DTest:
    """Minimal test case for debugging Fortran issues"""
    
    def __init__(self, lib_path=None):
        if lib_path is None:
            lib_path = "/app/earthimager/cli/build/libei2d_core.so"
        
        if not os.path.exists(lib_path):
            raise FileNotFoundError(f"Library not found: {lib_path}")
            
        self.lib = ctypes.CDLL(lib_path)
        self._setup_function_signatures()
    
    def _setup_function_signatures(self):
        """Setup C function signatures exactly as defined in C-ABI"""
        
        # ei2d_InitForwGlobals
        self.lib.ei2d_InitForwGlobals.argtypes = [
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
            ctypes.c_double, ctypes.c_double, ctypes.c_double
        ]
        self.lib.ei2d_InitForwGlobals.restype = None
        
        # ei2d_SetNumParamForward
        self.lib.ei2d_SetNumParamForward.argtypes = [ctypes.c_int32, ctypes.c_int32]
        self.lib.ei2d_SetNumParamForward.restype = None
        
        # ei2d_ForwardFD
        self.lib.ei2d_ForwardFD.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_int32), ctypes.c_int32,
            ctypes.c_int32, ctypes.c_int32, ctypes.c_int32
        ]
        self.lib.ei2d_ForwardFD.restype = None

    def create_minimal_test_case(self):
        """
        Create minimal 4-electrode, 2-measurement test case
        This is the simplest possible EI2D configuration for debugging
        """
        
        # MINIMAL CONFIGURATION: 4 electrodes in a line
        nElec = 4
        nInf = 1
        nNx = 4    # Match electrodes
        nNy = 3    # Minimal depth
        
        nNodes = nNx * nNy  # 12 nodes total
        nElem = (nNx - 1) * (nNy - 1)  # 6 elements total
        
        print(f"MINIMAL TEST CASE SETUP:")
        print(f"  Electrodes: {nElec}")
        print(f"  Nodes: {nNodes} ({nNx} x {nNy})")
        print(f"  Elements: {nElem} ({nNx-1} x {nNy-1})")
        
        # Create simple grid
        nodeX = np.zeros(nNodes, dtype=np.float64)
        nodeY = np.zeros(nNodes, dtype=np.float64)
        
        for j in range(nNy):
            for i in range(nNx):
                idx = j * nNx + i
                nodeX[idx] = i * 1.0  # 1m spacing
                nodeY[idx] = j * 1.0
        
        print(f"  Node coordinates created: {len(nodeX)} nodes")
        print(f"  X range: {nodeX.min()} to {nodeX.max()}")
        print(f"  Y range: {nodeY.min()} to {nodeY.max()}")
        
        # Homogeneous conductivity
        cond = np.full(nElem, 0.01, dtype=np.float64)  # 100 ohm-m
        print(f"  Conductivity: {len(cond)} elements, value = {cond[0]} S/m")
        
        # Electrode mapping (1-based for Fortran)
        elecNodeID = np.arange(1, nElec + nInf + 1, dtype=np.int32)
        elecNodeID[-1] = 1  # Infinity electrode maps to node 1
        print(f"  Electrode mapping: {elecNodeID}")
        
        # Minimal survey: 2 measurements (Wenner-like)
        # A-B-M-N configuration
        stingCMD = np.array([
            1, 2, 3, 4,  # First measurement: A=1, B=2, M=3, N=4
            1, 3, 2, 4   # Second measurement: A=1, B=3, M=2, N=4
        ], dtype=np.int32)
        nData = 2
        
        print(f"  Survey commands: {stingCMD.reshape(-1, 4)}")
        print(f"  Measurements: {nData}")
        
        # CRITICAL: Parameter regions - must match Fortran expectations
        # Based on Fortran signature: ParamX1/X2, ParamY1/Y2 are INTEGER arrays
        # They define parameter regions, not individual elements
        
        # SINGLE PARAMETER REGION covering the entire element grid
        nParamX = 1  # One region in X direction
        nParamY = 1  # One region in Y direction
        
        # Parameter region bounds (1-based, element indices)
        p1 = np.array([1], dtype=np.int32)              # Region starts at element X=1
        p2 = np.array([nNx - 1], dtype=np.int32)        # Region ends at element X=3 (for 4 nodes)
        q1 = np.array([1], dtype=np.int32)              # Region starts at element Y=1  
        q2 = np.array([nNy - 1], dtype=np.int32)        # Region ends at element Y=2 (for 3 nodes)
        
        print(f"  Parameter setup:")
        print(f"    nParamX={nParamX}, nParamY={nParamY}")
        print(f"    ParamX1={p1}, ParamX2={p2} (elements {p1[0]} to {p2[0]})")
        print(f"    ParamY1={q1}, ParamY2={q2} (elements {q1[0]} to {q2[0]})")
        print(f"    Total parameters: {nParamX * nParamY}")
        
        # Infinity electrodes
        inf = np.array([nElec + nInf], dtype=np.int32)
        print(f"  Infinity electrodes: {inf}")
        
        return {
            'nData': nData,
            'nElec': nElec,
            'nInf': nInf,
            'nNodes': nNodes,
            'nElem': nElem,
            'nNx': nNx,
            'nNy': nNy,
            'nodeX': nodeX,
            'nodeY': nodeY,
            'cond': cond,
            'elecNodeID': elecNodeID,
            'stingCMD': stingCMD,
            'nParamX': nParamX,
            'nParamY': nParamY,
            'p1': p1, 'p2': p2,
            'q1': q1, 'q2': q2,
            'inf': inf
        }
    
    def test_getjacobian_0(self):
        """
        Test GetJacobian=0 (should be simpler case)
        This should calculate V/I without sensitivity matrix
        """
        print(f"\n" + "="*60)
        print(f"TESTING GetJacobian=0 (V/I only, no sensitivity)")
        print(f"="*60)
        
        test_data = self.create_minimal_test_case()
        
        try:
            # Step 1: Initialize Forward Globals
            print(f"\n1. Calling ei2d_InitForwGlobals...")
            self.lib.ei2d_InitForwGlobals(
                test_data['nData'],           # NumData
                test_data['nElec'] + test_data['nInf'],  # NumElectrodes  
                test_data['nInf'],            # NumInfElectrodes
                test_data['nNx'],             # NumNodeX
                test_data['nNy'],             # NumNodeY
                0,                            # ForwModMeth (0=FD)
                0,                            # ForwSolver (0=Cholesky)
                0,                            # InvMethod
                0,                            # ForwAccuracy
                50,                           # ForwCGIter
                0,                            # BCType (0=Dirichlet)
                1e-6,                         # ForwCGResid
                0.0,                          # MinTxRxSep
                0.0                           # MaxTxRxSep
            )
            print("   ✓ InitForwGlobals completed")
            
            # Step 2: Set Parameter Dimensions
            print(f"\n2. Calling ei2d_SetNumParamForward...")
            self.lib.ei2d_SetNumParamForward(test_data['nParamX'], test_data['nParamY'])
            print(f"   ✓ SetNumParamForward({test_data['nParamX']}, {test_data['nParamY']}) completed")
            
            # Step 3: Prepare output arrays
            VI = np.zeros(test_data['nData'], dtype=np.float64)
            jacobian = np.zeros(test_data['nData'] * test_data['nParamX'] * test_data['nParamY'], dtype=np.float64)
            
            print(f"\n3. Calling ei2d_ForwardFD with GetJacobian=0...")
            print(f"   Output arrays: VI={len(VI)}, Jacobian={len(jacobian)}")
            
            # Step 4: Call ForwardFD with GetJacobian=0
            self.lib.ei2d_ForwardFD(
                test_data['nodeX'].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                test_data['nodeY'].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                test_data['cond'].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                VI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                jacobian.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                test_data['elecNodeID'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['stingCMD'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['p1'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['p2'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['q1'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['q2'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['inf'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                0,                            # GetJacobian=0
                test_data['nNodes'],          # nNodes
                test_data['nElem'],           # nElem
                test_data['nData']            # nData
            )
            
            print("   ✓ ForwardFD completed without crash")
            
            # Step 5: Analyze results
            print(f"\n4. RESULTS ANALYSIS:")
            print(f"   VI array length: {len(VI)}")
            print(f"   VI values: {VI}")
            print(f"   VI sum: {np.sum(VI)}")
            print(f"   VI non-zero count: {np.count_nonzero(VI)}")
            print(f"   VI NaN count: {np.sum(np.isnan(VI))}")
            print(f"   VI min/max: {np.min(VI)} / {np.max(VI)}")
            
            # Check if we got valid results
            valid_count = np.sum(np.isfinite(VI) & (VI != 0.0))
            if valid_count > 0:
                print(f"   🎉 SUCCESS: Got {valid_count} valid V/I values!")
                return True
            else:
                print(f"   ❌ PROBLEM: No valid V/I values (all zeros or NaN)")
                return False
                
        except Exception as e:
            print(f"   ❌ EXCEPTION: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    def test_getjacobian_1(self):
        """
        Test GetJacobian=1 (with sensitivity matrix)
        This may trigger the array bounds error
        """
        print(f"\n" + "="*60)
        print(f"TESTING GetJacobian=1 (V/I + Sensitivity)")
        print(f"="*60)
        
        test_data = self.create_minimal_test_case()
        
        try:
            # Step 1: Initialize Forward Globals  
            print(f"\n1. Calling ei2d_InitForwGlobals...")
            self.lib.ei2d_InitForwGlobals(
                test_data['nData'],
                test_data['nElec'] + test_data['nInf'],
                test_data['nInf'],
                test_data['nNx'],
                test_data['nNy'],
                0,  # ForwModMeth (0=FD)
                0,  # ForwSolver (0=Cholesky)
                0,  # InvMethod
                0,  # ForwAccuracy
                50, # ForwCGIter
                0,  # BCType (0=Dirichlet)
                1e-6, 0.0, 0.0
            )
            print("   ✓ InitForwGlobals completed")
            
            # Step 2: Set Parameter Dimensions
            print(f"\n2. Calling ei2d_SetNumParamForward...")
            self.lib.ei2d_SetNumParamForward(test_data['nParamX'], test_data['nParamY'])
            print("   ✓ SetNumParamForward completed")
            
            # Step 3: Prepare output arrays
            VI = np.zeros(test_data['nData'], dtype=np.float64)
            jacobian = np.zeros(test_data['nData'] * test_data['nParamX'] * test_data['nParamY'], dtype=np.float64)
            
            print(f"\n3. Calling ei2d_ForwardFD with GetJacobian=1...")
            print(f"   This may trigger the array bounds error in Sensitivity.f90")
            
            # Step 4: Call ForwardFD with GetJacobian=1
            self.lib.ei2d_ForwardFD(
                test_data['nodeX'].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                test_data['nodeY'].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                test_data['cond'].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                VI.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                jacobian.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                test_data['elecNodeID'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['stingCMD'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['p1'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['p2'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['q1'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['q2'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                test_data['inf'].ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                1,                            # GetJacobian=1
                test_data['nNodes'],
                test_data['nElem'], 
                test_data['nData']
            )
            
            print("   🎉 ForwardFD with GetJacobian=1 completed without crash!")
            
            # Step 5: Analyze results
            print(f"\n4. RESULTS ANALYSIS:")
            print(f"   VI values: {VI}")
            print(f"   VI non-zero: {np.count_nonzero(VI)}")
            print(f"   Jacobian sum: {np.sum(jacobian)}")
            print(f"   Jacobian non-zero: {np.count_nonzero(jacobian)}")
            
            return True
            
        except Exception as e:
            print(f"   ❌ EXCEPTION: {e}")
            import traceback
            traceback.print_exc()
            return False

def main():
    """Run minimal Fortran debugging tests"""
    print("MINIMAL FORTRAN DEBUGGING TEST")
    print("Focus: Isolate and fix core Fortran engine issues")
    print()
    
    try:
        tester = MinimalEI2DTest()
        
        # Test GetJacobian=0 first (simpler case)
        success_0 = tester.test_getjacobian_0()
        
        # Test GetJacobian=1 (may have array bounds issue)
        success_1 = tester.test_getjacobian_1()
        
        print(f"\n" + "="*60)
        print(f"SUMMARY:")
        print(f"  GetJacobian=0 (V/I only): {'✓ PASS' if success_0 else '❌ FAIL'}")
        print(f"  GetJacobian=1 (V/I+Jacobian): {'✓ PASS' if success_1 else '❌ FAIL'}")
        
        if success_0 and success_1:
            print(f"\n🎉 All tests passed! Fortran engine working correctly.")
        elif success_0:
            print(f"\n⚠️  GetJacobian=0 works, but GetJacobian=1 has issues.")
            print(f"   Focus debugging on Sensitivity.f90 array bounds.")
        else:
            print(f"\n💥 Both tests failed. Need to debug basic parameter setup.")
            
    except Exception as e:
        print(f"❌ Test setup failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()