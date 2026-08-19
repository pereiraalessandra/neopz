/// @file BindMaterialDarcy.cpp
/// @brief Bindings for the Darcy materials (Material/DarcyFlow/)

#include "Bindings.h"
#include "PybindMaterialHelpers.h"

#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZAnalyticSolution.h"
#include "TPZMaterial.h"

namespace {
/// TPZIsotropicPermeability is a virtual base: no pointer to member
template <class TMaterial>
void SetConstantPermeability(TMaterial &material, STATE permeability) {
  material.SetConstantPermeability(permeability);
}
} // namespace

void InitMaterialDarcy(py::module_ &m) {
  py::class_<TPZDarcyFlow, TPZMaterial, std::unique_ptr<TPZDarcyFlow, py::nodelete>>(
      m, "TPZDarcyFlow", "Darcy material for a single H1 space.")
      .def(py::init<int, int>(), py::arg("id"), py::arg("dim"))
      .def("SetConstantPermeability", &SetConstantPermeability<TPZDarcyFlow>, py::arg("permeability"))
      .def("SetForcingFunction", &SetForcingFunction<TPZDarcyFlow>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> list[float]")
      .def("SetForcingFunction", &SetForcingFunctionAnalytic<TPZDarcyFlow, TLaplaceExample1>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("SetExactSol", &SetExactSolCallback<TPZDarcyFlow>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> tuple[list[float], list[list[float]]]")
      .def("SetExactSol", &SetExactSolAnalytic<TPZDarcyFlow, TLaplaceExample1>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("HasExactSol", &HasExactSol<TPZDarcyFlow>)
      .def("NEvalErrors", &NEvalErrors<TPZDarcyFlow>)
      .def("CreateBC", &CreateBoundaryCondition<TPZDarcyFlow>, py::arg("id"), py::arg("type"), py::arg("val1"), py::arg("val2"), py::return_value_policy::reference, py::keep_alive<0, 1>());

  py::class_<TPZMixedDarcyFlow, TPZMaterial, std::unique_ptr<TPZMixedDarcyFlow, py::nodelete>>(
      m, "TPZMixedDarcyFlow", "Darcy material for combined spaces (HDiv flux, L2 pressure).")
      .def(py::init<int, int>(), py::arg("id"), py::arg("dim"))
      .def("SetConstantPermeability", &SetConstantPermeability<TPZMixedDarcyFlow>, py::arg("permeability"))
      .def("SetForcingFunction", &SetForcingFunction<TPZMixedDarcyFlow>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> list[float]")
      .def("SetForcingFunction", &SetForcingFunctionAnalytic<TPZMixedDarcyFlow, TLaplaceExample1>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("SetExactSol", &SetExactSolCallback<TPZMixedDarcyFlow>, py::arg("f"), py::arg("pOrder") = 1, "f(loc: list[float]) -> tuple[list[float], list[list[float]]]")
      .def("SetExactSol", &SetExactSolAnalytic<TPZMixedDarcyFlow, TLaplaceExample1>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("HasExactSol", &HasExactSol<TPZMixedDarcyFlow>)
      .def("NEvalErrors", &NEvalErrors<TPZMixedDarcyFlow>)
      .def("CreateBC", &CreateBoundaryCondition<TPZMixedDarcyFlow>, py::arg("id"), py::arg("type"), py::arg("val1"), py::arg("val2"), py::return_value_policy::reference, py::keep_alive<0, 1>());
}
