/// @file BindAnalyticSolutionDarcy.cpp
/// @brief Bindings for the built-in Darcy/Laplace exact solutions
/// (Pre/TPZAnalyticSolution.h)

#include "Bindings.h"

#include "TPZAnalyticSolution.h"

void InitAnalyticSolutionDarcy(py::module_ &m) {
  auto laplaceExample = py::class_<TLaplaceExample1>(
      m, "TLaplaceExample1", "Built-in exact solutions for Darcy/Laplace problems.");

  py::enum_<TLaplaceExample1::EExactSol>(laplaceExample, "EExactSol")
      .value("ENone", TLaplaceExample1::ENone)
      .value("EConst", TLaplaceExample1::EConst)
      .value("EX", TLaplaceExample1::EX)
      .value("EY", TLaplaceExample1::EY)
      .value("EZ", TLaplaceExample1::EZ)
      .value("EXpY", TLaplaceExample1::EXpY)
      .value("EX2", TLaplaceExample1::EX2)
      .value("ESinSin", TLaplaceExample1::ESinSin)
      .value("ECosCos", TLaplaceExample1::ECosCos)
      .value("EArcTan", TLaplaceExample1::EArcTan)
      .value("EArcTanSingular", TLaplaceExample1::EArcTanSingular)
      .value("ESteepWave", TLaplaceExample1::ESteepWave)
      .value("ESteepWave2", TLaplaceExample1::ESteepWave2)
      .value("ESinDist", TLaplaceExample1::ESinDist)
      .value("E10SinSin", TLaplaceExample1::E10SinSin)
      .value("E2SinSin", TLaplaceExample1::E2SinSin)
      .value("ESinSinDirNonHom", TLaplaceExample1::ESinSinDirNonHom)
      .value("ESinMark", TLaplaceExample1::ESinMark)
      .value("ESinMark2", TLaplaceExample1::ESinMark2)
      .value("ECosMark", TLaplaceExample1::ECosMark)
      .value("ESteklovNonConst", TLaplaceExample1::ESteklovNonConst)
      .value("ESteklovNonConst2", TLaplaceExample1::ESteklovNonConst2)
      .value("EPerpendicularSteklovNonConst", TLaplaceExample1::EPerpendicularSteklovNonConst)
      .value("EGalvisNonConst", TLaplaceExample1::EGalvisNonConst)
      .value("EBoundaryLayer", TLaplaceExample1::EBoundaryLayer)
      .value("EBubble", TLaplaceExample1::EBubble)
      .value("EBubble2D", TLaplaceExample1::EBubble2D)
      .value("ESinCosCircle", TLaplaceExample1::ESinCosCircle)
      .value("EHarmonic", TLaplaceExample1::EHarmonic)
      .value("EHarmonic2", TLaplaceExample1::EHarmonic2)
      .value("ESquareRootUpper", TLaplaceExample1::ESquareRootUpper)
      .value("ESquareRootLower", TLaplaceExample1::ESquareRootLower)
      .value("ESquareRoot", TLaplaceExample1::ESquareRoot)
      .value("ELaplace2D", TLaplaceExample1::ELaplace2D)
      .value("EHarmonic3", TLaplaceExample1::EHarmonic3)
      .value("EHarmonicPoly", TLaplaceExample1::EHarmonicPoly)
      .value("ESharpGaussian2D", TLaplaceExample1::ESharpGaussian2D);

  laplaceExample.def(py::init<>())
      .def_readwrite("fExact", &TLaplaceExample1::fExact)
      .def_readwrite("fDimension", &TLaplaceExample1::fDimension)
      .def_readwrite("fSignConvention", &TLaplaceExample1::fSignConvention)
      .def_readwrite("fCenter", &TLaplaceExample1::fCenter, "Center point used by solutions defined around one, "
                                                            "such as EArcTan.")
      .def("Name", &TLaplaceExample1::Name)
      // Force and Solution write into force[0] / u[0] without resizing: the
      // caller is expected to pre-size the output, one state variable here.
      .def("Force", [](const TLaplaceExample1 &self, const TPZVec<REAL> &x) {
                 TPZManVector<STATE> force(1);
                 self.Force(x, force);
                 return force; }, py::arg("x"), "Source term consistent with the exact solution.")
      .def("Solution", [](const TLaplaceExample1 &self, const TPZVec<REAL> &x) {
                 TPZManVector<STATE> u(1);
                 TPZFMatrix<STATE> gradu;
                 self.Solution(x, u, gradu);
                 return std::make_pair(u, gradu); }, py::arg("x"), "Returns (u, gradu) at the given point.");
}
