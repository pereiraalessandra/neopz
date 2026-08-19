/// @file BindMaterialElasticity.cpp
/// @brief Bindings for the elasticity material (Material/Elasticity/)

#include "Bindings.h"
#include "PybindMaterialHelpers.h"

#include "Elasticity/TPZElasticity2D.h"
#include "Elasticity/TPZElasticity3D.h"
#include "Elasticity/TPZMixedElasticityND.h"
#include "TPZAnalyticSolution.h"
#include "TPZMaterial.h"

void InitMaterialElasticity(py::module_ &m) {
  py::class_<TPZElasticity2D, TPZMaterial, std::unique_ptr<TPZElasticity2D, py::nodelete>>(
      m, "TPZElasticity2D", "Elasticity material for a single H1 space, 2D.")
      .def(py::init<int>(), py::arg("id"))
      .def("SetElasticity", &TPZElasticity2D::SetElasticity, py::arg("E"), py::arg("nu"))
      .def("SetPlaneStrain", &TPZElasticity2D::SetPlaneStrain)
      .def("SetPlaneStress", &TPZElasticity2D::SetPlaneStress)
      .def("SetForcingFunction", &SetForcingFunction<TPZElasticity2D>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> list[float]")
      .def("SetForcingFunction", &SetForcingFunctionAnalytic<TPZElasticity2D, TElasticity2DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("SetExactSol", &SetExactSolCallback<TPZElasticity2D>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> tuple[list[float], list[list[float]]]")
      .def("SetExactSol", &SetExactSolAnalytic<TPZElasticity2D, TElasticity2DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("HasExactSol", &HasExactSol<TPZElasticity2D>)
      .def("NEvalErrors", &NEvalErrors<TPZElasticity2D>)
      .def("CreateBC", &CreateBoundaryCondition<TPZElasticity2D>, py::arg("id"), py::arg("type"), py::arg("val1"), py::arg("val2"), py::return_value_policy::reference, py::keep_alive<0, 1>());

  // No SetPlaneStrain/SetPlaneStress: that concept is 2D only.
  py::class_<TPZElasticity3D, TPZMaterial, std::unique_ptr<TPZElasticity3D, py::nodelete>>(
      m, "TPZElasticity3D", "Elasticity material for a single H1 space, 3D.")
      .def(py::init<int>(), py::arg("id"))
      // Kept as SetMaterialDataHook, matching the C++ name exactly (unlike
      // every other elasticity material's SetElasticity(E, nu)) so code
      // ported from C++ NeoPZ translates one to one.
      .def("SetMaterialDataHook", &TPZElasticity3D::SetMaterialDataHook, py::arg("E"), py::arg("nu"))
      .def("SetForcingFunction", &SetForcingFunction<TPZElasticity3D>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> list[float]")
      .def("SetForcingFunction", &SetForcingFunctionAnalytic<TPZElasticity3D, TElasticity3DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("SetExactSol", &SetExactSolCallback<TPZElasticity3D>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> tuple[list[float], list[list[float]]]")
      .def("SetExactSol", &SetExactSolAnalytic<TPZElasticity3D, TElasticity3DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("HasExactSol", &HasExactSol<TPZElasticity3D>)
      .def("NEvalErrors", &NEvalErrors<TPZElasticity3D>)
      .def("CreateBC", &CreateBoundaryCondition<TPZElasticity3D>, py::arg("id"), py::arg("type"), py::arg("val1"), py::arg("val2"), py::return_value_policy::reference, py::keep_alive<0, 1>());

  py::class_<TPZMixedElasticityND, TPZMaterial, std::unique_ptr<TPZMixedElasticityND, py::nodelete>>(
      m, "TPZMixedElasticityND", "Elasticity material for combined spaces (stress, displacement, rotation).")
      .def(py::init<int, int>(), py::arg("id"), py::arg("dim") = 2)
      .def("SetElasticity", &TPZMixedElasticityND::SetElasticity, py::arg("E"), py::arg("nu"))
      .def("SetPlaneStrain", &TPZMixedElasticityND::SetPlaneStrain)
      .def("SetPlaneStress", &TPZMixedElasticityND::SetPlaneStress)
      .def("SetForcingFunction", &SetForcingFunction<TPZMixedElasticityND>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> list[float]")
      .def("SetForcingFunction", &SetForcingFunctionAnalytic<TPZMixedElasticityND, TElasticity2DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("SetForcingFunction", &SetForcingFunctionAnalytic<TPZMixedElasticityND, TElasticity3DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("SetExactSol", &SetExactSolCallback<TPZMixedElasticityND>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> tuple[list[float], list[list[float]]]")
      .def("SetExactSol", &SetExactSolAnalytic<TPZMixedElasticityND, TElasticity2DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("SetExactSol", &SetExactSolAnalytic<TPZMixedElasticityND, TElasticity3DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("HasExactSol", &HasExactSol<TPZMixedElasticityND>)
      .def("NEvalErrors", &NEvalErrors<TPZMixedElasticityND>)
      .def("CreateBC", &CreateBoundaryCondition<TPZMixedElasticityND>, py::arg("id"), py::arg("type"), py::arg("val1"), py::arg("val2"), py::return_value_policy::reference, py::keep_alive<0, 1>());
}
