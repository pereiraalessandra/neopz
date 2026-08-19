/// @file PybindCallbacks.h
/// @brief Adapts Python callables to the std::function signatures NeoPZ expects.
/// @note NeoPZ calls these from assembly, which may be multi-threaded and hence
/// running without the GIL; every wrapper here reacquires it before calling Python.

#ifndef PZ_PYBIND_CALLBACKS_H
#define PZ_PYBIND_CALLBACKS_H

#include <pybind11/pybind11.h>

#include "PybindCasters.h"
#include "TPZMatError.h"
#include "TPZMatTypes.h"

namespace py = pybind11;

/// wraps loc -> result into a ForcingFunctionType; loc and result convert
/// through the TPZVec caster on the way in and out.
template <class TVar>
ForcingFunctionType<TVar> MakeForcingFunction(py::function func) {
  return [func](const TPZVec<REAL> &loc, TPZVec<TVar> &result) {
    py::gil_scoped_acquire gil;
    result = func(loc).template cast<TPZVec<TVar>>();
  };
}

/// wraps loc -> (result, deriv) into an ExactSolType; func must return a
/// 2-tuple, converted through the TPZVec and TPZFMatrix casters. The shape
/// matches what TLaplaceExample1.Solution(loc) already returns, so an analytic
/// solution's own Solution method can be passed here unchanged.
template <class TVar>
ExactSolType<TVar> MakeExactSolution(py::function func) {
  return [func](const TPZVec<REAL> &loc, TPZVec<TVar> &result, TPZFMatrix<TVar> &deriv) {
    py::gil_scoped_acquire gil;
    py::tuple ret = func(loc);
    result = ret[0].template cast<TPZVec<TVar>>();
    deriv = ret[1].template cast<TPZFMatrix<TVar>>();
  };
}

#endif // PZ_PYBIND_CALLBACKS_H
