/// @file Bindings.h
/// @brief Shared includes and the declaration of every init function

#ifndef PZ_PYBIND_BINDINGS_H
#define PZ_PYBIND_BINDINGS_H

#include "PybindCasters.h"
#include "PybindUtils.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

/// one init function per binding file, called from Module.cpp in dependency order
void InitEnums(py::module_ &m);
void InitGeom(py::module_ &m);
void InitAnalyticSolutionDarcy(py::module_ &m);
void InitAnalyticSolutionElasticity(py::module_ &m);
void InitMaterial(py::module_ &m);
void InitMaterialDarcy(py::module_ &m);
void InitMaterialElasticity(py::module_ &m);
void InitMesh(py::module_ &m);
void InitApproxCreator(py::module_ &m);
void InitAnalysis(py::module_ &m);
void InitPost(py::module_ &m);

#endif // PZ_PYBIND_BINDINGS_H