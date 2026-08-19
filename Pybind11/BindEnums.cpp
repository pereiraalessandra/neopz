/// @file BindEnums.cpp
/// @brief Bindings for the enumerations used by the other binding files
/// @note export_values() is deliberately not used: ENone appears in three of
/// these enums and exporting would make one shadow the others

#include "Bindings.h"

#include "MMeshType.h"
#include "TPZApproxCreator.h"
#include "TPZEnumApproxFamily.h"
#include "pzbasematrix.h"

void InitEnums(py::module_ &m) {
  py::enum_<HybridizationType>(m, "HybridizationType")
      .value("ENone", HybridizationType::ENone)
      .value("EStandard", HybridizationType::EStandard)
      .value("EStandardSquared", HybridizationType::EStandardSquared)
      .value("ESemi", HybridizationType::ESemi);

  py::enum_<ProblemType>(m, "ProblemType")
      .value("ENone", ProblemType::ENone)
      .value("EElastic", ProblemType::EElastic)
      .value("EDarcy", ProblemType::EDarcy)
      .value("EStokes", ProblemType::EStokes);

  py::enum_<HDivFamily>(m, "HDivFamily")
      .value("EHDivStandard", HDivFamily::EHDivStandard)
      .value("EHDivConstant", HDivFamily::EHDivConstant)
      .value("EHDivKernel", HDivFamily::EHDivKernel)
      .value("EHDivOptimized", HDivFamily::EHDivOptimized);

  py::enum_<H1Family>(m, "H1Family")
      .value("EH1Standard", H1Family::EH1Standard)
      .value("EH1WidePrism", H1Family::EH1WidePrism);

  py::enum_<MMeshType>(m, "MMeshType")
      .value("EQuadrilateral", MMeshType::EQuadrilateral)
      .value("ETriangular", MMeshType::ETriangular)
      .value("EHexahedral", MMeshType::EHexahedral)
      .value("ETetrahedral", MMeshType::ETetrahedral)
      .value("EPyramidal", MMeshType::EPyramidal)
      .value("EPrismatic", MMeshType::EPrismatic)
      .value("EHexaPyrMixed", MMeshType::EHexaPyrMixed)
      .value("ENoType", MMeshType::ENoType);

  py::enum_<DecomposeType>(m, "DecomposeType")
      .value("ENoDecompose", DecomposeType::ENoDecompose)
      .value("ELU", DecomposeType::ELU)
      .value("ELUPivot", DecomposeType::ELUPivot)
      .value("ECholesky", DecomposeType::ECholesky)
      .value("ELDLt", DecomposeType::ELDLt);
}
