/// @file BindApproxCreator.cpp
/// @brief Bindings for TPZApproxCreator and its derived classes (Pre/)
/// @note keep_alive guards the pointers the creators store, and every method
/// returning a mesh keeps its creator alive. TPZApproxCreator.h only forward
/// declares the material and mesh types

#include "Bindings.h"

#include "TPZApproxCreator.h"
#include "TPZBndCond.h"
#include "TPZH1ApproxCreator.h"
#include "TPZHDivApproxCreator.h"
#include "TPZMHMApproxCreator.h"
#include "TPZMHMHDivApproxCreator.h"
#include "TPZMaterial.h"
#include "TPZMultiphysicsCompMesh.h"

void InitApproxCreator(py::module_ &m) {
  py::class_<TPZHybrid>(m, "TPZHybrid")
      .def(py::init<>())
      .def(py::init<int64_t, int64_t, int64_t>(), py::arg("lagrange"), py::arg("left"), py::arg("right"))
      .def_readwrite("fLagrange", &TPZHybrid::fLagrange)
      .def_readwrite("fLeft", &TPZHybrid::fLeft)
      .def_readwrite("fRight", &TPZHybrid::fRight);

  auto approxCreator = py::class_<TPZApproxCreator>(m, "TPZApproxCreator");

  py::class_<TPZApproxCreator::HybridizationData>(
      approxCreator, "HybridizationData", "Material ids and multipliers describing a hybridization. The list "
                                          "and dict fields convert by copy: read, change, assign back.")
      .def(py::init<>())
      .def(py::init<const TPZApproxCreator::HybridizationData &>(), py::arg("other"))
      .def("__copy__", [](const TPZApproxCreator::HybridizationData &data) { return TPZApproxCreator::HybridizationData(data); })
      .def_readwrite("fWrapMatId", &TPZApproxCreator::HybridizationData::fWrapMatId)
      .def_readwrite("fLeftInterfaceMatId", &TPZApproxCreator::HybridizationData::fLeftInterfaceMatId)
      .def_readwrite("fRightInterfaceMatId", &TPZApproxCreator::HybridizationData::fRightInterfaceMatId)
      .def_readwrite("fLagrangeMatId", &TPZApproxCreator::HybridizationData::fLagrangeMatId)
      .def_readwrite("fSecondLeftInterfaceMatId", &TPZApproxCreator::HybridizationData::fSecondLeftInterfaceMatId)
      .def_readwrite("fSecondRightInterfaceMatId", &TPZApproxCreator::HybridizationData::fSecondRightInterfaceMatId)
      .def_readwrite("fSecondLagrangeMatId", &TPZApproxCreator::HybridizationData::fSecondLagrangeMatId)
      .def_readwrite("fHybridizeBCLevel", &TPZApproxCreator::HybridizationData::fHybridizeBCLevel)
      .def_readwrite("fMultipliers", &TPZApproxCreator::HybridizationData::fMultipliers)
      .def_readwrite("fInterfaces", &TPZApproxCreator::HybridizationData::fInterfaces)
      .def("SetProblemHybridH1", &TPZApproxCreator::HybridizationData::SetProblemHybridH1, py::arg("prob"), py::arg("hybrid"))
      .def("SetProblemHybridHDiv", &TPZApproxCreator::HybridizationData::SetProblemHybridHDiv, py::arg("prob"), py::arg("hybrid"))
      .def("__str__", [](const TPZApproxCreator::HybridizationData &data) { return PrintToString(data); });

  approxCreator
      .def("SetHybridType", &TPZApproxCreator::SetHybridType, py::arg("hybrid"))
      .def("HybridType", &TPZApproxCreator::HybridType)
      .def("SetProbType", &TPZApproxCreator::SetProbType, py::arg("prob"))
      .def("ProbType", &TPZApproxCreator::ProbType)
      .def("SetNState", &TPZApproxCreator::SetNState, py::arg("nstate"))
      .def("NState", &TPZApproxCreator::NState)
      .def("NumMeshes", &TPZApproxCreator::NumMeshes)
      .def("SetDefaultOrder", &TPZApproxCreator::SetDefaultOrder, py::arg("order"))
      .def("GetDefaultOrder", &TPZApproxCreator::GetDefaultOrder)
      .def("SetExtraInternalOrder", &TPZApproxCreator::SetExtraInternalOrder, py::arg("order"), "0 is HDiv, 1 HDiv+, 2 HDiv++.")
      .def("GetExtraInternalOrder", &TPZApproxCreator::GetExtraInternalOrder)
      .def("SetShouldCondense", &TPZApproxCreator::SetShouldCondense, py::arg("isCondensed"))
      /// GetShouldCondense takes an argument it ignores; kept in the
      /// signature to match the C++ call site exactly.
      .def("GetShouldCondense", &TPZApproxCreator::GetShouldCondense, py::arg("isCondensed"))
      // IsRigidBodySpaces() is a C++ reference accessor (creator.IsRigidBodySpaces()
      // = value); Python has no equivalent assignable-call-result syntax. The
      // getter keeps the identical call-with-parens read (creator.IsRigidBodySpaces()),
      // and a second, value-taking overload serves as the setter.
      .def("IsRigidBodySpaces", [](const TPZApproxCreator &creator) { return creator.IsRigidBodySpaces(); }, "Enriches the space with constant fields.")
      .def("IsRigidBodySpaces", [](TPZApproxCreator &creator, bool value) { creator.IsRigidBodySpaces() = value; }, py::arg("value"))
      .def("HybridData", &TPZApproxCreator::HybridData, "Copy of the hybridization data, filled by SetHybridType.")
      .def("SetHybridData", &TPZApproxCreator::SetHybridData, py::arg("data"))
      .def("SetHybridizeBoundary", &TPZApproxCreator::SetHybridizeBoundary, "Requires a hybridization type to have been set.")
      .def("SetMeshElementType", &TPZApproxCreator::SetMeshElementType)
      .def("InsertMaterialObject", py::overload_cast<TPZMaterial *>(&TPZApproxCreator::InsertMaterialObject), py::arg("material"), py::keep_alive<1, 2>())
      .def("InsertMaterialObject", py::overload_cast<TPZBndCond *>(&TPZApproxCreator::InsertMaterialObject), py::arg("bndcond"), py::keep_alive<1, 2>())
      .def("GeoMesh", py::overload_cast<>(&TPZApproxCreator::GeoMesh), py::return_value_policy::reference)
      .def("CreateApproximationSpace", &TPZApproxCreator::CreateApproximationSpace, py::keep_alive<0, 1>(), "Creates the atomic meshes and the multiphysics mesh.")
      .def("__str__", [](TPZApproxCreator &creator) { return PrintToString(creator); });

  py::class_<TPZH1ApproxCreator, TPZApproxCreator>(m, "TPZH1ApproxCreator")
      .def(py::init<>())
      .def(py::init<TPZGeoMesh *>(), py::arg("gmesh"), py::keep_alive<1, 2>())
      .def("CreateClassicH1ApproximationSpace", &TPZH1ApproxCreator::CreateClassicH1ApproximationSpace, py::keep_alive<0, 1>(), "Single H1 mesh. Only valid without hybridization.")
      .def("CreateAtomicMeshes", [](TPZH1ApproxCreator &self) {
                 TPZManVector<TPZCompMesh *> meshvec;
                 self.CreateAtomicMeshes(meshvec);
                 return meshvec; }, py::return_value_policy::reference, "Atomic meshes. Requires a hybridization type.")
      .def("GroupElements", &TPZH1ApproxCreator::GroupElements, py::arg("cmesh"))
      .def("CondenseElements", &TPZH1ApproxCreator::CondenseElements, py::arg("cmesh"))
      .def("GroupAndCondenseElements", &TPZH1ApproxCreator::GroupAndCondenseElements, py::arg("cmesh"));

  py::class_<TPZHDivApproxCreator, TPZApproxCreator>(m, "TPZHDivApproxCreator")
      .def(py::init<>())
      .def(py::init<TPZGeoMesh *>(), py::arg("gmesh"), py::keep_alive<1, 2>())
      // Same reference-accessor situation as IsRigidBodySpaces above: two
      // overloads instead of a property, keeping the getter's call syntax
      // (with parens) identical to the C++ read.
      .def("HdivFamily", [](const TPZHDivApproxCreator &self) { return self.HdivFamily(); })
      .def("HdivFamily", [](TPZHDivApproxCreator &self, HDivFamily family) { self.HdivFamily() = family; }, py::arg("family"))
      .def_static("PrintAllMeshes", &TPZHDivApproxCreator::PrintAllMeshes, py::arg("cmesh"))
      .def("SetShouldCondensePressure", &TPZHDivApproxCreator::SetShouldCondensePressure, py::arg("isCondensed"), "Darcy only.")
      .def("ShouldCondensePressure", &TPZHDivApproxCreator::ShouldCondensePressure)
      .def("CheckSetupConsistency", &TPZHDivApproxCreator::CheckSetupConsistency, "Raises RuntimeError if the configuration is invalid.")
      .def("AddInterfaceComputationalElements", &TPZHDivApproxCreator::AddInterfaceComputationalElements, py::arg("cmesh"))
      .def("AddInterfaceComputationalElementsBackup", &TPZHDivApproxCreator::AddInterfaceComputationalElementsBackup, py::arg("cmesh"), "Earlier implementation, kept for debugging.")
      .def("PrintMeshElementsConnectInfo", &TPZHDivApproxCreator::PrintMeshElementsConnectInfo, py::arg("cmesh"))
      .def("CreateAtomicMeshes", [](TPZHDivApproxCreator &self, int lagLevelCounter) {
                 TPZManVector<TPZCompMesh *, 7> meshvec;
                 self.CreateAtomicMeshes(meshvec, lagLevelCounter);
                 return std::make_pair(meshvec, lagLevelCounter); }, py::arg("lagLevelCounter") = 1, py::return_value_policy::reference, "Returns (meshes, lagLevelCounter).")
      .def("CreateMultiPhysicsMesh", [](TPZHDivApproxCreator &self, TPZManVector<TPZCompMesh *, 7> meshvec, int lagLevelCounter) {
                 TPZMultiphysicsCompMesh *cmeshmulti = nullptr;
                 self.CreateMultiPhysicsMesh(meshvec, lagLevelCounter, cmeshmulti);
                 return cmeshmulti; }, py::arg("meshvec"), py::arg("lagLevelCounter") = 1, py::keep_alive<0, 1>(), "Builds the multiphysics mesh from CreateAtomicMeshes.");

  /// taken by value: a caster cannot propagate through a non const ref
  py::class_<TPZMHMApproxCreator>(m, "TPZMHMApproxCreator", "Element partition used by the MHM creators.")
      .def(py::init<>())
      .def(py::init([](TPZVec<int64_t> elementPartition) { return new TPZMHMApproxCreator(elementPartition); }), py::arg("elementPartition"));

  py::class_<TPZMHMHDivApproxCreator, TPZHDivApproxCreator, TPZMHMApproxCreator>(
      m, "TPZMHMHDivApproxCreator")
      .def(py::init<TPZGeoMesh *>(), py::arg("gmesh"), py::keep_alive<1, 2>())
      .def(py::init([](TPZGeoMesh *gmesh, TPZVec<int64_t> elementPartition) { return new TPZMHMHDivApproxCreator(gmesh, elementPartition); }), py::arg("gmesh"), py::arg("elementPartition"), py::keep_alive<1, 2>())
      .def("BuildMultiphysicsCMesh", &TPZMHMHDivApproxCreator::BuildMultiphysicsCMesh, py::keep_alive<0, 1>(), "MHM driver. Call SetPOrderSkeleton first: the C++ constructors "
                                                                                                               "leave that member uninitialized.")
      .def("PutinSubstructures", &TPZMHMHDivApproxCreator::PutinSubstructures, py::arg("cmesh"))
      .def("CondenseElements", &TPZMHMHDivApproxCreator::CondenseElements, py::arg("cmesh"))
      .def("SetPOrderSkeleton", &TPZMHMHDivApproxCreator::SetPOrderSkeleton, py::arg("order"))
      .def("GetPOrderSkeleton", &TPZMHMHDivApproxCreator::GetPOrderSkeleton);
}
