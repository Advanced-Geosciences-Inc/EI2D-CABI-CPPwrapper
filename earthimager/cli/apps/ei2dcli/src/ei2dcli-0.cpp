#include <cstdio>
#include <vector>
#include <cstdint>
#include "ei2d_api.h"

int main() {
  const int32_t nNx=2, nNy=3;                 // 2 x 3 grid
  const int32_t nNodes = nNx*nNy;             // 6
  const int32_t nElem  = (nNx-1)*(nNy-1);     // 2
  const int32_t nData  = 1;                   // single datum is simpler
  const int32_t nParamX=1, nParamY=1;

  ei2d_InitForwGlobals(nData, 6, 0, nNx, nNy,
                       0, 0, 0, 0, 50, 0,
                       1e-6, 0.0, 0.0);
  ei2d_SetNumParamForward(nParamX, nParamY);

  // nodes (row-major): (x,y): (0,0),(1,0),(0,1),(1,1),(0,2),(1,2)
  std::vector<double> nodeX{0,1, 0,1, 0,1};
  std::vector<double> nodeY{0,0, 1,1, 2,2};
  std::vector<double> cond(nElem, 0.01);

  std::vector<double> VI(nData, 0.0), Jac(nData * (nParamX*nParamY), 0.0);

  // 4 electrodes mapped to valid nodes (1-based indices)
  std::vector<int32_t> elecNodeID{1,2,3,4,4,4};

  // A=1, B=2, M=3, N=4
  std::vector<int32_t> stingCMD{1,2,3,4};

  // param windows (trivial)
  std::vector<int32_t> p1(nParamX,1), p2(nParamX,nNx);
  std::vector<int32_t> q1(nParamY,1), q2(nParamY,nNy);
  std::vector<int32_t> inf(1,0);

  ei2d_ForwardFD(nodeX.data(), nodeY.data(), cond.data(),
                 VI.data(), Jac.data(),
                 elecNodeID.data(), stingCMD.data(),
                 p1.data(), p2.data(), q1.data(), q2.data(),
                 inf.data(), /*GetJacobian*/ 0,
                 nNodes, nElem, nData);

  std::printf("SELFTEST OK. VI[0]=%.6g\n", VI[0]);
  return 0;
}
