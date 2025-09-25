#include <cstdio>
#include <cstdint>
#include <vector>
#include <algorithm>
#include "ei2d_api.h"

int main() {
  // --- geometry / mesh ---
  const int32_t nElec = 8;   // number of physical electrodes along the surface
  const int32_t nInf  = 1;   // one “infinite” electrode slot
  const double  a     = 1.0; // spacing
  const int32_t nNx   = nElec; // top row matches electrodes (nodes 1..nNx)
  const int32_t nNy   = 6;     // depth layers

  const int32_t nNodes = nNx * nNy;
  const int32_t nElem  = (nNx - 1) * (nNy - 1);

  // node coordinates (row-major; top row is nodes 1..nNx in Fortran)
  std::vector<double> nodeX(nNodes), nodeY(nNodes);
  for (int32_t j = 0; j < nNy; ++j) {
    for (int32_t i = 0; i < nNx; ++i) {
      const int32_t id = j * nNx + i;
      nodeX[id] = i * a;
      nodeY[id] = j * a;
    }
  }

  // homogeneous half-space: rho=100 ohm-m => cond = 0.01 S/m
  std::vector<double> cond(nElem, 0.01);

  // --- electrodes (Fortran 1-based node indices) ---
  // size = nElec + nInf; physical electrodes map to nodes 1..nElec; the "inf" slot gets a dummy node index (1)
  std::vector<int32_t> elecNodeID(nElec + nInf);
  for (int i = 0; i < nElec; ++i) elecNodeID[i] = i + 1; // 1..8
  elecNodeID[nElec] = 1; // dummy for the ∞ slot

  // which entry (1-based) is the infinite electrode? the last one
  std::vector<int32_t> inf(nInf);
  inf[0] = nElec + nInf; // = 9

  // --- build dipole–dipole ABMN for n=1..2 (A=i, B=i+1, M=i+1+n, N=i+2+n) ---
  std::vector<int32_t> stingCMD; stingCMD.reserve(4 * 32);
  for (int32_t n = 1; n <= 2; ++n) {
    for (int32_t i = 1; i <= nElec - (2 + n); ++i) {
      const int32_t A = i;
      const int32_t B = i + 1;
      const int32_t M = i + 1 + n;
      const int32_t N = i + 2 + n;
      stingCMD.push_back(A);
      stingCMD.push_back(B);
      stingCMD.push_back(M);
      stingCMD.push_back(N);
    }
  }
  const int32_t nData = static_cast<int32_t>(stingCMD.size() / 4);

  // --- one global parameter window covering entire mesh ---
  const int32_t nParamX = 1, nParamY = 1;
  std::vector<int32_t> p1{1}, p2{nNx}, q1{1}, q2{nNy};

  // --- init engine ---
  ei2d_InitForwGlobals(
    /*NumData*/          nData,
    /*NumElectrodes*/    nElec + nInf,
    /*NumInfElectrodes*/ nInf,
    /*NumNodeX*/         nNx,
    /*NumNodeY*/         nNy,
    /*ForwModMeth*/      0,   // FD
    /*ForwSolver*/       0,   // Cholesky
    /*InvMethod*/        0,
    /*ForwAccuracy*/     0,
    /*ForwCGIter*/       50,
    /*BCType*/           0,
    /*ForwCGResid*/      1e-6,
    /*MinTxRxSep*/       0.0,
    /*MaxTxRxSep*/       0.0
  );
  ei2d_SetNumParamForward(nParamX, nParamY);

  // --- outputs ---
  std::vector<double> VI(nData, 0.0);
  std::vector<double> J(nData * (nParamX * nParamY), 0.0);

  // sanity guards
  int maxElecID = 0;
  for (int v : elecNodeID) if (v > maxElecID) maxElecID = v;
  int maxCmd = 0;
  for (size_t i = 0; i < stingCMD.size(); ++i) if (stingCMD[i] > maxCmd) maxCmd = stingCMD[i];
  std::printf("DBG: nElec=%d, nData=%d, max(ElecNodeID)=%d, max(StingCMD)=%d\n",
              nElec, nData, maxElecID, maxCmd);

  // --- run forward (no Jacobian) ---
  ei2d_ForwardFD(
    nodeX.data(), nodeY.data(), cond.data(),
    VI.data(), J.data(),
    elecNodeID.data(), stingCMD.data(),
    p1.data(), p2.data(), q1.data(), q2.data(),
    inf.data(), /*GetJacobian*/ 0,
    nNodes, nElem, nData
  );

  std::printf("SYNTH OK. nData=%d, first few V/I: ", nData);
for (int i = 0; i < std::min<int>(5, nData); ++i) std::printf("%.6g ", VI[i]);
std::printf("\n");

// Save all results to CSV
FILE* f = std::fopen("build/forward_out.csv", "w");
if (f) {
  std::fprintf(f, "i,A,B,M,N,VI\n");
  for (int i = 0; i < nData; ++i) {
    int A = stingCMD[4*i+0], B = stingCMD[4*i+1];
    int M = stingCMD[4*i+2], N = stingCMD[4*i+3];
    std::fprintf(f, "%d,%d,%d,%d,%d,%.17g\n", i+1, A,B,M,N, VI[i]);
  }
  std::fclose(f);
  std::puts("Wrote build/forward_out.csv");
}

}
