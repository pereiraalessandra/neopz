/// @file Module.cpp
/// @brief Entry point of the pyneopz extension module
/// @note the only file carrying PYBIND11_MODULE; the others expose init functions

#include "Bindings.h"

PYBIND11_MODULE(pyneopz, m) {
  m.doc() = "Python bindings for NeoPZ";

  /// dependency order: bases and enums before whatever mentions them.
  /// Each AnalyticSolutionX must precede InitMaterial (whose TPZBndCondT
  /// binds SetForcingFunctionBC for every analytic type) and its own
  /// InitMaterialX (whose SetExactSol overloads take the analytic type).
  InitEnums(m);
  InitGeom(m);
  InitAnalyticSolutionDarcy(m);
  InitAnalyticSolutionElasticity(m);
  InitMaterial(m);
  InitMaterialDarcy(m);
  InitMaterialElasticity(m);
  InitMesh(m);
  InitApproxCreator(m);
  InitAnalysis(m);
  InitPost(m);
}