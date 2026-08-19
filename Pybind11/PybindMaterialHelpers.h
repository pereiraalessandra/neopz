/// @file PybindMaterialHelpers.h
/// @brief Templates shared by every concrete material's binding file
/// (BindMaterialDarcy.cpp, BindMaterialElasticity.cpp, ...): CreateBC,
/// SetForcingFunction and SetExactSol (each with callback and analytic-object
/// overloads), HasExactSol, NEvalErrors.
/// @note SetExactSol/HasExactSol/NEvalErrors sit behind a virtual base
/// (TPZMatError): these templates call them directly, not through a
/// pointer to member.

#ifndef PZ_PYBIND_MATERIAL_HELPERS_H
#define PZ_PYBIND_MATERIAL_HELPERS_H

#include "PybindCallbacks.h"

#include "TPZBndCondT.h"
#include "pzfmatrix.h"

/// TPZMaterialT::CreateBC returns TPZBndCondT<STATE>*, not the more generic
/// TPZBndCond*. Matching that return type here is what makes pybind resolve
/// the object as a TPZBndCondT — needed for SetForcingFunctionBC.
template <class TMaterial>
TPZBndCondT<STATE> *CreateBoundaryCondition(TMaterial &material, int id, int type, const TPZFMatrix<STATE> &val1, const TPZVec<STATE> &val2) {
  return material.CreateBC(&material, id, type, val1, val2);
}

template <class TMaterial>
void SetForcingFunction(TMaterial &material, py::function f, int pOrder) {
  material.SetForcingFunction(MakeForcingFunction<STATE>(f), pOrder);
}

/// analytic.ForceFunc() captures &analytic: same lifetime requirement as
/// SetExactSolAnalytic below, enforced the same way (keep_alive<1, 2>).
template <class TMaterial, class TAnalytic>
void SetForcingFunctionAnalytic(TMaterial &material, TAnalytic &analytic, int pOrder) {
  material.SetForcingFunction(analytic.ForceFunc(), pOrder);
}

template <class TMaterial>
void SetExactSolCallback(TMaterial &material, py::function f, int pOrder) {
  material.SetExactSol(MakeExactSolution<STATE>(f), pOrder);
}

/// analytic.ExactSolution() captures &analytic: the material must not
/// outlive it, enforced by keep_alive<1, 2> on the caller's .def().
/// Templated on the analytic type: TLaplaceExample1 for Darcy,
/// TElasticity2DAnalytic/3DAnalytic for elasticity.
template <class TMaterial, class TAnalytic>
void SetExactSolAnalytic(TMaterial &material, TAnalytic &analytic, int pOrder) {
  material.SetExactSol(analytic.ExactSolution(), pOrder);
}

template <class TMaterial>
bool HasExactSol(const TMaterial &material) {
  return material.HasExactSol();
}

template <class TMaterial>
int NEvalErrors(const TMaterial &material) {
  return material.NEvalErrors();
}

#endif // PZ_PYBIND_MATERIAL_HELPERS_H
