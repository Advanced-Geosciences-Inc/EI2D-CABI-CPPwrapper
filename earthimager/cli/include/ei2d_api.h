#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void ei2d_CheckKinds(void);

void ei2d_InitForwGlobals(int32_t NumData, int32_t NumElectrodes, int32_t NumInfElectrodes,
                          int32_t NumNodeX, int32_t NumNodeY,
                          int32_t ForwModMeth, int32_t ForwSolver, int32_t InvMethod,
                          int32_t ForwAccuracy, int32_t ForwCGIter, int32_t BCType,
                          double ForwCGResid, double MinTxRxSep, double MaxTxRxSep);

void ei2d_SetNumParamForward(int32_t NumParamX, int32_t NumParamY);

void ei2d_ForwardFD(const double* NodeX, const double* NodeY, const double* Cond,
                    double* VIcalc, double* Jacobian,
                    const int32_t* ElecNodeID,
                    const int32_t* StingCMD,
                    const int32_t* ParamX1, const int32_t* ParamX2, const int32_t* ParamY1, const int32_t* ParamY2,
                    const int32_t* InfElec,
                    int32_t GetJacobian,
                    int32_t nNodes, int32_t nElem, int32_t nData);

void ei2d_ForwardFE(const double* NodeX, const double* NodeY, const double* Cond,
                    double* VIcalc, double* Jacobian,
                    const int32_t* ElecNodeID,
                    const int32_t* StingCMD,
                    const int32_t* ParamX1, const int32_t* ParamX2, const int32_t* ParamY1, const int32_t* ParamY2,
                    const int32_t* InfElec,
                    const double* CenterNodeX, const double* CenterNodeY, const double* ElemArea,
                    int32_t GetJacobian,
                    int32_t nNodes, int32_t nElem, int32_t nData);

void ei2d_FreeMemForward(void);

void ei2d_InitInvGlobals(int32_t NumData, int32_t NumElemX, int32_t NumElemY,
                         int32_t NumParamX, int32_t NumParamY,
                         int32_t InvMethod, int32_t IPInvMethod, int32_t MaxNumIterInvCG, int32_t IPPosMeth,
                         double ModResoFactor, double EpsilonD, double EpsilonM);

void ei2d_InvPCGLS(const double* ObsData, const double* CalcData, const double* DataWeight,
                   const double* ModelParam, const double* PriorModel, double* ModelUpdate,
                   const double* Jacobian, double DampingFactor,
                   int32_t ResIPFlag, int32_t IterNum,
                   int32_t nData, int32_t nParam);

void ei2d_InvPCGOC(const double* ObsData, const double* CalcData, const double* DataWeight,
                   const double* ModelParam, const double* PriorModel, double* ModelUpdate,
                   const double* Jacobian,
                   const double* RoughX0, const double* RoughY0, double* RoughX, double* RoughY,
                   double Lagrange, double DampFactor,
                   int32_t ResIPFlag, int32_t IterNum,
                   int32_t nData, int32_t nParam);

void ei2d_InvPCGRobust(const double* ObsData, const double* CalcData, const double* DataWeight,
                       const double* ModelParam, const double* PriorModel, double* ModelUpdate,
                       const double* Jacobian,
                       const double* RoughX0, const double* RoughY0, double* RoughX, double* RoughY,
                       double Lagrange,
                       int32_t ResIPFlag, int32_t IterNum,
                       int32_t nData, int32_t nParam);

void ei2d_UpdateJacobian(const double* p, const double* DataDiff, double* Jacobian,
                         int32_t nParam, int32_t nData);

#ifdef __cplusplus
}
#endif
