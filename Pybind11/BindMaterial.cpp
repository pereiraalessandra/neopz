/// @file BindMaterial.cpp
/// @brief Bindings for the material interfaces shared by every concrete
/// material (Material/TPZMaterial.h, TPZBndCond.h, TPZBndCondT.h). Concrete
/// materials live in one file per physics family: BindMaterialDarcy.cpp,
/// BindMaterialElasticity.cpp, ...
/// @note materials inserted into a mesh are deleted by TPZCompMesh::CleanUp,
/// hence the py::nodelete holder used throughout those files. TPZBndCond is a
/// sibling of TPZMaterial, not derived from it, which is why
/// InsertMaterialObject is overloaded on both

#include "Bindings.h"
#include "PybindCallbacks.h"

#include "TPZAnalyticSolution.h"
#include "TPZBndCond.h"
#include "TPZBndCondT.h"
#include "TPZMaterial.h"

namespace {
/// ForcingFunctionBCType has the same (loc, result, deriv) shape as
/// ExactSolType, so MakeExactSolution is reused unchanged.
void SetForcingFunctionBCCallback(TPZBndCondT<STATE> &bc, py::function f, int pOrder) {
  bc.SetForcingFunctionBC(MakeExactSolution<STATE>(f), pOrder);
}

template <class TAnalytic>
void SetForcingFunctionBCAnalytic(TPZBndCondT<STATE> &bc, TAnalytic &analytic, int pOrder) {
  bc.SetForcingFunctionBC(analytic.ExactSolution(), pOrder);
}
} // namespace

void InitMaterial(py::module_ &m) {
  py::class_<TPZMaterial>(m, "TPZMaterial", "Interface of every NeoPZ material.")
      .def("Id", &TPZMaterial::Id)
      .def("NStateVariables", &TPZMaterial::NStateVariables)
      .def("Dimension", &TPZMaterial::Dimension)
      .def("HasForcingFunction", &TPZMaterial::HasForcingFunction)
      .def("__str__", [](const TPZMaterial &material) { return PrintToString(material); });

  py::class_<TPZBndCond>(m, "TPZBndCond", "Type agnostic interface of a boundary condition.")
      .def("Id", &TPZBndCond::Id)
      .def("Type", &TPZBndCond::Type, "0 is Dirichlet, 1 is Neumann.")
      .def("__str__", [](const TPZBndCond &bndcond) { return PrintToString(bndcond); });

  // What CreateBC actually returns at runtime: pybind resolves the dynamic
  // type via RTTI, so this class need not be named explicitly at the call
  // site for its methods to become available on the returned object.
  py::class_<TPZBndCondT<STATE>, TPZBndCond>(
      m, "TPZBndCondT", "Boundary condition with a concrete state type.")
      .def("SetForcingFunctionBC", &SetForcingFunctionBCCallback, py::arg("f"), py::arg("pOrder") = 1, "Position-dependent boundary values, needed whenever the exact "
                                                                                                       "solution is not constant. "
                                                                                                       "f(loc: list[float]) -> tuple[list[float], list[list[float]]]")
      .def("SetForcingFunctionBC", &SetForcingFunctionBCAnalytic<TLaplaceExample1>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("SetForcingFunctionBC", &SetForcingFunctionBCAnalytic<TElasticity2DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("SetForcingFunctionBC", &SetForcingFunctionBCAnalytic<TElasticity3DAnalytic>, py::arg("analytic"), py::arg("pOrder") = 1, py::keep_alive<1, 2>())
      .def("Val1", &TPZBndCondT<STATE>::Val1)
      .def("Val2", &TPZBndCondT<STATE>::Val2);
}
