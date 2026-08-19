/// @file BindAnalyticSolutionElasticity.cpp
/// @brief Bindings for the built-in 2D/3D elasticity exact solutions
/// (Pre/TPZAnalyticSolution.h)

#include "Bindings.h"

#include "TPZAnalyticSolution.h"

void InitAnalyticSolutionElasticity(py::module_ &m) {
  auto elasticity2D = py::class_<TElasticity2DAnalytic>(
      m, "TElasticity2DAnalytic", "Built-in exact solutions for 2D elasticity problems.");

  py::enum_<TElasticity2DAnalytic::EDefState>(elasticity2D, "EDefState")
      .value("ENone", TElasticity2DAnalytic::ENone)
      .value("EDispx", TElasticity2DAnalytic::EDispx)
      .value("EDispy", TElasticity2DAnalytic::EDispy)
      .value("EDispxy", TElasticity2DAnalytic::EDispxy)
      .value("ERot", TElasticity2DAnalytic::ERot)
      .value("EStretchx", TElasticity2DAnalytic::EStretchx)
      .value("EUniAxialx", TElasticity2DAnalytic::EUniAxialx)
      .value("EStretchy", TElasticity2DAnalytic::EStretchy)
      .value("EShear", TElasticity2DAnalytic::EShear)
      .value("EHomogeneous", TElasticity2DAnalytic::EHomogeneous)
      .value("EBend", TElasticity2DAnalytic::EBend)
      .value("ELoadedBeam", TElasticity2DAnalytic::ELoadedBeam)
      .value("Etest1", TElasticity2DAnalytic::Etest1)
      .value("Etest2", TElasticity2DAnalytic::Etest2)
      .value("EThiago", TElasticity2DAnalytic::EThiago)
      .value("EPoly", TElasticity2DAnalytic::EPoly)
      .value("ESquareRootUpper", TElasticity2DAnalytic::ESquareRootUpper)
      .value("ESquareRootLower", TElasticity2DAnalytic::ESquareRootLower)
      .value("ESquareRoot", TElasticity2DAnalytic::ESquareRoot);

  elasticity2D.def(py::init<>())
      .def_readwrite("fProblemType", &TElasticity2DAnalytic::fProblemType)
      .def_readwrite("fPlaneStress", &TElasticity2DAnalytic::fPlaneStress, "1 for plane stress, 0 for plane strain.")
      .def_readwrite("fSignConvention", &TElasticity2DAnalytic::fSignConvention)
      // gE and gPoisson are static: shared by every TElasticity2DAnalytic instance
      .def_readwrite_static("gE", &TElasticity2DAnalytic::gE)
      .def_readwrite_static("gPoisson", &TElasticity2DAnalytic::gPoisson)
      .def("Name", &TElasticity2DAnalytic::Name)
      // Force/Solution write into force[0..1] / u[0..1] without resizing.
      .def("Force", [](const TElasticity2DAnalytic &self, const TPZVec<REAL> &x) {
                 TPZManVector<STATE> force(2);
                 self.Force(x, force);
                 return force; }, py::arg("x"), "Source term consistent with the exact solution.")
      .def("Solution", [](const TElasticity2DAnalytic &self, const TPZVec<REAL> &x) {
                 TPZManVector<STATE> u(2);
                 TPZFMatrix<STATE> gradu;
                 self.Solution(x, u, gradu);
                 return std::make_pair(u, gradu); }, py::arg("x"), "Returns (u, gradu) at the given point.");

  auto elasticity3D = py::class_<TElasticity3DAnalytic>(
      m, "TElasticity3DAnalytic", "Built-in exact solutions for 3D elasticity problems.");

  py::enum_<TElasticity3DAnalytic::EDefState>(elasticity3D, "EDefState")
      .value("ENone", TElasticity3DAnalytic::ENone)
      .value("EDispx", TElasticity3DAnalytic::EDispx)
      .value("EDispy", TElasticity3DAnalytic::EDispy)
      .value("EDispxyz", TElasticity3DAnalytic::EDispxyz)
      .value("ERot", TElasticity3DAnalytic::ERot)
      .value("ERotXYZ", TElasticity3DAnalytic::ERotXYZ)
      .value("EStretchx", TElasticity3DAnalytic::EStretchx)
      .value("EUniAxialx", TElasticity3DAnalytic::EUniAxialx)
      .value("EStretchy", TElasticity3DAnalytic::EStretchy)
      .value("EStretchz", TElasticity3DAnalytic::EStretchz)
      .value("EHomogeneous", TElasticity3DAnalytic::EHomogeneous)
      .value("EShear", TElasticity3DAnalytic::EShear)
      .value("EBend", TElasticity3DAnalytic::EBend)
      .value("ELoadedBeam", TElasticity3DAnalytic::ELoadedBeam)
      .value("Etest1", TElasticity3DAnalytic::Etest1)
      .value("Etest2", TElasticity3DAnalytic::Etest2)
      .value("ETestShearMoment", TElasticity3DAnalytic::ETestShearMoment)
      .value("ESphere", TElasticity3DAnalytic::ESphere)
      .value("EYotov", TElasticity3DAnalytic::EYotov);

  // fE/fPoisson are instance members here, unlike TElasticity2DAnalytic's
  // static gE/gPoisson.
  elasticity3D.def(py::init<>())
      .def_readwrite("fProblemType", &TElasticity3DAnalytic::fProblemType)
      .def_readwrite("fE", &TElasticity3DAnalytic::fE)
      .def_readwrite("fPoisson", &TElasticity3DAnalytic::fPoisson)
      .def_readwrite("fSignConvention", &TElasticity3DAnalytic::fSignConvention)
      .def("Force", [](const TElasticity3DAnalytic &self, const TPZVec<REAL> &x) {
                 TPZManVector<STATE> force(3);
                 self.Force(x, force);
                 return force; }, py::arg("x"), "Source term consistent with the exact solution.")
      .def("Solution", [](const TElasticity3DAnalytic &self, const TPZVec<REAL> &x) {
                 TPZManVector<STATE> u(3);
                 TPZFMatrix<STATE> gradu;
                 self.Solution(x, u, gradu);
                 return std::make_pair(u, gradu); }, py::arg("x"), "Returns (u, gradu) at the given point.");
}
