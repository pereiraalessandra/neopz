#
# Created by Alessandra on 10/08/2026.
#
"""Regression tests for the NeoPZ Python bindings.

Uses unittest from the standard library: no external dependency is needed.
To run it:

    export PYTHONPATH=path/to/build/python:$PYTHONPATH
    python3 -m unittest TestPyNeoPZ -v
"""

# ----- Unit test includes -----
import gc
import math
import os
import tempfile
import unittest

# ----- 3rd party includes -----
import numpy as np

# ----- PZ includes -----
import pyneopz


# ==========> Helpers <==========
# ===============================

# Minimal Gmsh 4.1 ASCII mesh: unit square, one triangulated domain and
# four boundary edges.
SQUARE_MSH = """$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
5
1 1 "bottom"
1 2 "right"
1 3 "top"
1 4 "left"
2 5 "domain"
$EndPhysicalNames
$Entities
4 4 1 0
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
1 0 0 0 1 0 0 1 1 2 1 -2
2 1 0 0 1 1 0 1 2 2 2 -3
3 0 1 0 1 1 0 1 3 2 3 -4
4 0 0 0 0 1 0 1 4 2 4 -1
1 0 0 0 1 1 0 1 5 4 1 2 3 4
$EndEntities
$Nodes
9 5 1 5
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
1 1 0
0 4 0 1
4
0 1 0
1 1 0 0
1 2 0 0
1 3 0 0
1 4 0 0
2 1 0 1
5
0.5 0.5 0
$EndNodes
$Elements
5 8 1 8
1 1 1 1
1 1 2
1 2 1 1
2 2 3
1 3 1 1
3 3 4
1 4 1 1
4 4 1
2 1 2 4
5 1 2 5
6 4 1 5
7 2 3 5
8 3 4 5
$EndElements
"""


def WriteSquareMsh():
    """Writes SQUARE_MSH to a temp file and returns its path."""
    fd, path = tempfile.mkstemp(suffix=".msh")
    with os.fdopen(fd, "w") as f:
        f.write(SQUARE_MSH)
    return path


def CreateGeoMesh3D():
    """Unit cube, one tetrahedron subdivision, six boundary faces (-1..-6)."""
    gen = pyneopz.TPZGenGrid3D(
        [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1, 1, 1], pyneopz.MMeshType.ETetrahedral
    )
    gen.BuildVolumetricElements(1)
    return gen.BuildBoundaryElements(-1, -2, -3, -4, -5, -6)


def EvaluateConstantSolutionErrorOnMesh(gmesh, darcy, boundaryMatIds, pOrder=2):
    """Solves a Darcy problem with constant exact solution 1 (ExactSol
    assumed already set) and Dirichlet BCs on boundaryMatIds; returns the
    error vector. For meshes not from TPZGenGrid2D (Gmsh, TPZGenGrid3D)."""
    creator = pyneopz.TPZHDivApproxCreator(gmesh)
    creator.SetProbType(pyneopz.ProblemType.EDarcy)
    creator.SetDefaultOrder(pOrder)
    creator.InsertMaterialObject(darcy)
    for matid in boundaryMatIds:
        creator.InsertMaterialObject(darcy.CreateBC(matid, 0, [[0.0]], [1.0]))
    cmesh = creator.CreateApproximationSpace()

    an = pyneopz.TPZLinearAnalysis(cmesh)
    an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
    solver = pyneopz.TPZStepSolver()
    solver.SetDirect(pyneopz.DecomposeType.ELDLt)
    an.SetSolver(solver)
    an.Assemble()
    an.Solve()

    pyneopz.TransferFromMultiPhysics(cmesh.MeshVector(), cmesh)
    cmesh.LoadReferences()
    return cmesh.EvaluateError(False, darcy.NEvalErrors())


def Create2DGeoMesh(nDiv=2):
    """Unit square with nDiv divisions per direction and Dirichlet boundary."""
    gmesh = pyneopz.TPZGeoMesh()
    gen = pyneopz.TPZGenGrid2D([nDiv, nDiv], [0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
    gen.SetElementType(pyneopz.MMeshType.EQuadrilateral)
    gen.Read(gmesh, 1)
    for side in [4, 5, 6, 7]:
        gen.SetBC(gmesh, side, -1)
    gmesh.BuildConnectivity()
    return gmesh


def CreateH1ApproximationSpace(gmesh, pOrder=1):
    """Classic H1 space built on the given geometric mesh."""
    darcy = pyneopz.TPZDarcyFlow(1, 2)
    darcy.SetConstantPermeability(1.0)
    bc = darcy.CreateBC(-1, 0, [[0.0]], [1.0])

    creator = pyneopz.TPZH1ApproxCreator(gmesh)
    creator.SetProbType(pyneopz.ProblemType.EDarcy)
    creator.SetDefaultOrder(pOrder)
    creator.InsertMaterialObject(darcy)
    creator.InsertMaterialObject(bc)
    return creator.CreateClassicH1ApproximationSpace()


def CreateHDivApproximationSpace(gmesh, pOrder=1):
    """Mixed Darcy approximation space built on the given geometric mesh."""
    darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
    darcy.SetConstantPermeability(1.0)
    bc = darcy.CreateBC(-1, 0, [[0.0]], [1.0])

    creator = pyneopz.TPZHDivApproxCreator(gmesh)
    creator.SetProbType(pyneopz.ProblemType.EDarcy)
    creator.SetDefaultOrder(pOrder)
    creator.InsertMaterialObject(darcy)
    creator.InsertMaterialObject(bc)
    return creator.CreateApproximationSpace()


def SolveDirect(cmesh, decomposeType=pyneopz.DecomposeType.ELDLt):
    an = pyneopz.TPZLinearAnalysis(cmesh)
    an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
    solver = pyneopz.TPZStepSolver()
    solver.SetDirect(decomposeType)
    an.SetSolver(solver)
    an.Assemble()
    an.Solve()
    return an


def SolveIteratively(cmesh, useCG, decomposeType=pyneopz.DecomposeType.ELDLt):
    """Warms up with a direct solve, builds a preconditioner from a snapshot
    of that matrix (see MakePreconditionerFrom), then re-solves with CG or
    GMRES."""
    an = pyneopz.TPZLinearAnalysis(cmesh)
    an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
    warmup = pyneopz.TPZStepSolver()
    warmup.SetDirect(decomposeType)
    an.SetSolver(warmup)
    an.Assemble()

    precond = pyneopz.MakePreconditionerFrom(an, decomposeType)
    solver = pyneopz.TPZStepSolver()
    if useCG:
        solver.SetCG(500, precond, 1e-12, 0)
    else:
        solver.SetGMRES(200, 30, precond, 1e-12, 0)
    an.SetSolver(solver)
    an.Assemble()
    an.Solve()
    return an


def EvaluateConstantSolutionError(darcy):
    """Solves a Darcy problem with constant exact solution 1 (ExactSol
    assumed already set); returns the error vector from darcy.NEvalErrors()."""
    bc = darcy.CreateBC(-1, 0, [[0.0]], [1.0])
    creator = pyneopz.TPZHDivApproxCreator(Create2DGeoMesh())
    creator.SetProbType(pyneopz.ProblemType.EDarcy)
    creator.SetDefaultOrder(2)
    creator.InsertMaterialObject(darcy)
    creator.InsertMaterialObject(bc)
    cmesh = creator.CreateApproximationSpace()

    an = pyneopz.TPZLinearAnalysis(cmesh)
    an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
    solver = pyneopz.TPZStepSolver()
    solver.SetDirect(pyneopz.DecomposeType.ELDLt)
    an.SetSolver(solver)
    an.Assemble()
    an.Solve()

    pyneopz.TransferFromMultiPhysics(cmesh.MeshVector(), cmesh)
    cmesh.LoadReferences()
    return cmesh.EvaluateError(False, darcy.NEvalErrors())


def EvaluateStretchXError(setExact):
    """Solves a 2D plane-strain elasticity problem with exact solution
    EStretchx (u = (x, 0), constant strain); returns the error vector from
    mat.NEvalErrors(). setExact(mat, bc) calls SetExactSol/SetForcingFunctionBC,
    letting the same problem run through either the TElasticity2DAnalytic
    overload or a Python callback."""
    mat = pyneopz.TPZMixedElasticityND(1, 2)
    mat.SetElasticity(1.0, 0.0)
    mat.SetPlaneStrain()
    bc = mat.CreateBC(-1, 0, [[0.0, 0.0], [0.0, 0.0]], [0.0, 0.0])
    setExact(mat, bc)

    creator = pyneopz.TPZHDivApproxCreator(Create2DGeoMesh())
    creator.SetProbType(pyneopz.ProblemType.EElastic)
    creator.SetDefaultOrder(2)
    creator.InsertMaterialObject(mat)
    creator.InsertMaterialObject(bc)
    cmesh = creator.CreateApproximationSpace()

    an = pyneopz.TPZLinearAnalysis(cmesh)
    an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
    solver = pyneopz.TPZStepSolver()
    solver.SetDirect(pyneopz.DecomposeType.ELDLt)
    an.SetSolver(solver)
    an.Assemble()
    an.Solve()

    pyneopz.TransferFromMultiPhysics(cmesh.MeshVector(), cmesh)
    cmesh.LoadReferences()
    return mat.NEvalErrors(), cmesh.EvaluateError(False, mat.NEvalErrors())


# ==========> Test cases <==========
# ==================================

class TestEnums(unittest.TestCase):
    def test_do_not_leak_into_the_module(self):
        # without export_values(), ENone only exists qualified by its enum
        self.assertFalse(hasattr(pyneopz, "ENone"))
        self.assertTrue(hasattr(pyneopz.ProblemType, "ENone"))
        self.assertTrue(hasattr(pyneopz.HybridizationType, "ENone"))

    def test_values(self):
        self.assertEqual(int(pyneopz.HDivFamily.EHDivConstant), 1)
        self.assertEqual(
            list(pyneopz.H1Family.__members__), ["EH1Standard", "EH1WidePrism"]
        )


class TestGeometry(unittest.TestCase):
    def test_generated_mesh(self):
        gmesh = Create2DGeoMesh()
        # 4 domain quadrilaterals + 8 boundary edges (2 per side)
        self.assertEqual(gmesh.NElements(), 12)
        self.assertEqual(gmesh.NNodes(), 9)
        self.assertEqual(gmesh.Dimension(), 2)

    def test_tpzvec_caster(self):
        gen = pyneopz.TPZGenGrid2D([3, 3], [0.0, 0.0, 0.0], [2.0, 2.0, 0.0])
        gmesh = pyneopz.TPZGeoMesh()
        gen.Read(gmesh, 1)
        self.assertEqual(gmesh.NElements(), 9)


class TestGmshReader(unittest.TestCase):
    def test_loads_physical_names(self):
        reader = pyneopz.TPZGmshReader()
        gmesh = reader.GeometricGmshMesh(WriteSquareMsh())
        # 4 boundary edges + 4 triangles (Gmsh's own triangulation of the square)
        self.assertEqual(gmesh.NElements(), 8)
        self.assertEqual(gmesh.NNodes(), 5)
        self.assertEqual(gmesh.Dimension(), 2)

        names = reader.GetDimPhysicalTagName()
        self.assertEqual(names[1], {1: "bottom", 2: "right", 3: "top", 4: "left"})
        self.assertEqual(names[2], {5: "domain"})

    def test_error_is_near_zero_on_a_mesh_read_from_file(self):
        # same verification style as TestErrorEvaluation, applied to a mesh
        # that did not come from TPZGenGrid2D
        reader = pyneopz.TPZGmshReader()
        gmesh = reader.GeometricGmshMesh(WriteSquareMsh())

        example = pyneopz.TLaplaceExample1()
        example.fExact = pyneopz.TLaplaceExample1.EExactSol.EConst
        darcy = pyneopz.TPZMixedDarcyFlow(5, 2)  # 5 = "domain"
        darcy.SetConstantPermeability(1.0)
        darcy.SetExactSol(example, 2)

        error = EvaluateConstantSolutionErrorOnMesh(gmesh, darcy, [1, 2, 3, 4])
        for component in error:
            self.assertLess(component, 1e-9)


class TestGenGrid3D(unittest.TestCase):
    def test_generated_mesh(self):
        gmesh = CreateGeoMesh3D()
        self.assertEqual(gmesh.Dimension(), 3)
        # one cube split into tets (volume) plus its six boundary faces
        self.assertEqual(gmesh.NElements(), 17)
        self.assertEqual(gmesh.NNodes(), 8)

    def test_error_is_near_zero_in_3d(self):
        gmesh = CreateGeoMesh3D()
        example = pyneopz.TLaplaceExample1()
        example.fExact = pyneopz.TLaplaceExample1.EExactSol.EConst
        darcy = pyneopz.TPZMixedDarcyFlow(1, 3)
        darcy.SetConstantPermeability(1.0)
        darcy.SetExactSol(example, 2)

        error = EvaluateConstantSolutionErrorOnMesh(
            gmesh, darcy, [-1, -2, -3, -4, -5, -6], pOrder=1
        )
        for component in error:
            self.assertLess(component, 1e-9)


class TestCheckGeom(unittest.TestCase):
    def test_uniform_refine_element_counts(self):
        # 4 domain quads + 8 boundary edges = 12; each pass subdivides every
        # leaf element only (quads into 4, edges into 2); counts include
        # the now-divided parents that stay in the element vector
        gmesh = Create2DGeoMesh()
        self.assertEqual(gmesh.NElements(), 12)

        pyneopz.TPZCheckGeom(gmesh).UniformRefine(1)
        self.assertEqual(gmesh.NElements(), 44)

        pyneopz.TPZCheckGeom(gmesh).UniformRefine(1)
        self.assertEqual(gmesh.NElements(), 140)

    def test_refined_mesh_still_solves(self):
        gmesh = Create2DGeoMesh()
        pyneopz.TPZCheckGeom(gmesh).UniformRefine(1)

        example = pyneopz.TLaplaceExample1()
        example.fExact = pyneopz.TLaplaceExample1.EExactSol.EConst
        darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
        darcy.SetConstantPermeability(1.0)
        darcy.SetExactSol(example, 2)

        error = EvaluateConstantSolutionErrorOnMesh(gmesh, darcy, [-1], pOrder=1)
        for component in error:
            self.assertLess(component, 1e-9)


class TestMaterial(unittest.TestCase):
    def test_abstract_classes_have_no_constructor(self):
        with self.assertRaises(TypeError):
            pyneopz.TPZMaterial()
        with self.assertRaises(TypeError):
            pyneopz.TPZApproxCreator()

    def test_create_bc_returns_tpzbndcond(self):
        darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
        bc = darcy.CreateBC(-1, 0, [[0.0]], [1.0])
        self.assertIsInstance(bc, pyneopz.TPZBndCond)
        self.assertEqual(bc.Id(), -1)

    def test_create_bc_returns_the_concrete_type(self):
        # CreateBoundaryCondition must declare TPZBndCondT<STATE>* and not the
        # more generic TPZBndCond*, or pybind loses the concrete type on
        # return and methods like SetForcingFunctionBC become unreachable
        darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
        bc = darcy.CreateBC(-1, 0, [[0.0]], [1.0])
        self.assertIsInstance(bc, pyneopz.TPZBndCondT)

    def test_has_forcing_function_reflects_state(self):
        darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
        self.assertFalse(darcy.HasForcingFunction())
        darcy.SetForcingFunction(lambda loc: [1.0], 1)
        self.assertTrue(darcy.HasForcingFunction())

    def test_forcing_function_affects_the_solution(self):
        # homogeneous Dirichlet with no source: the solution is identically
        # zero, and any change in norm traces back to the Python callback alone
        def SolutionNorm(darcy):
            bc = darcy.CreateBC(-1, 0, [[0.0]], [0.0])
            creator = pyneopz.TPZHDivApproxCreator(Create2DGeoMesh())
            creator.SetProbType(pyneopz.ProblemType.EDarcy)
            creator.SetDefaultOrder(2)
            creator.InsertMaterialObject(darcy)
            creator.InsertMaterialObject(bc)
            cmesh = creator.CreateApproximationSpace()
            an = pyneopz.TPZLinearAnalysis(cmesh)
            an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
            solver = pyneopz.TPZStepSolver()
            solver.SetDirect(pyneopz.DecomposeType.ELDLt)
            an.SetSolver(solver)
            an.Assemble()
            an.Solve()
            return an.SolutionNorm()

        withoutSource = pyneopz.TPZMixedDarcyFlow(1, 2)
        withoutSource.SetConstantPermeability(1.0)

        withSource = pyneopz.TPZMixedDarcyFlow(1, 2)
        withSource.SetConstantPermeability(1.0)
        withSource.SetForcingFunction(lambda loc: [loc[0] * loc[1]], 2)

        self.assertEqual(SolutionNorm(withoutSource), 0.0)
        self.assertGreater(SolutionNorm(withSource), 0.0)


class TestMixedElasticityND(unittest.TestCase):
    def test_has_exact_sol_reflects_state(self):
        mat = pyneopz.TPZMixedElasticityND(1, 2)
        mat.SetElasticity(1.0, 0.0)
        self.assertFalse(mat.HasExactSol())
        example = pyneopz.TElasticity2DAnalytic()
        mat.SetExactSol(example, 2)
        self.assertTrue(mat.HasExactSol())

    def test_n_eval_errors(self):
        mat = pyneopz.TPZMixedElasticityND(1, 2)
        # pressure/stress, sigma-eps residual, div(sigma), displacement,
        # rotation, symmetry, and the exact solution's own energy norm
        self.assertEqual(mat.NEvalErrors(), 7)

    def test_create_bc_returns_tpzbndcondt(self):
        mat = pyneopz.TPZMixedElasticityND(1, 2)
        bc = mat.CreateBC(-1, 0, [[0.0, 0.0], [0.0, 0.0]], [0.0, 0.0])
        self.assertIsInstance(bc, pyneopz.TPZBndCondT)


class TestElasticity2D(unittest.TestCase):
    # single-element mesh: sidesteps a NeoPZ core limitation in classic H1
    # elasticity on multi-element meshes.
    @staticmethod
    def SolveSingleElementRigidTranslation(ux):
        gmesh = pyneopz.TPZGeoMesh()
        gen = pyneopz.TPZGenGrid2D([1, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        gen.SetElementType(pyneopz.MMeshType.EQuadrilateral)
        gen.Read(gmesh, 1)
        for side in [4, 5, 6, 7]:
            gen.SetBC(gmesh, side, -1)
        gmesh.BuildConnectivity()

        mat = pyneopz.TPZElasticity2D(1)
        mat.SetElasticity(1.0, 0.0)
        mat.SetPlaneStrain()
        bc = mat.CreateBC(-1, 0, [[0.0, 0.0], [0.0, 0.0]], [ux, 0.0])

        creator = pyneopz.TPZH1ApproxCreator(gmesh)
        creator.SetProbType(pyneopz.ProblemType.EElastic)
        creator.SetNState(2)
        creator.SetDefaultOrder(1)
        creator.InsertMaterialObject(mat)
        creator.InsertMaterialObject(bc)
        cmesh = creator.CreateClassicH1ApproximationSpace()

        an = pyneopz.TPZLinearAnalysis(cmesh)
        an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
        solver = pyneopz.TPZStepSolver()
        solver.SetDirect(pyneopz.DecomposeType.ELDLt)
        an.SetSolver(solver)
        an.Assemble()
        an.Solve()
        return an

    def test_rigid_translation_solution_is_exact(self):
        an = self.SolveSingleElementRigidTranslation(1.0)
        # four nodes, u=(1,0) at every one: norm = sqrt(4*1^2) = 2
        self.assertAlmostEqual(an.SolutionNorm(), 2.0, places=10)

    def test_has_exact_sol_reflects_state(self):
        mat = pyneopz.TPZElasticity2D(1)
        mat.SetElasticity(1.0, 0.0)
        self.assertFalse(mat.HasExactSol())
        example = pyneopz.TElasticity2DAnalytic()
        mat.SetExactSol(example, 2)
        self.assertTrue(mat.HasExactSol())

    def test_n_eval_errors(self):
        mat = pyneopz.TPZElasticity2D(1)
        self.assertEqual(mat.NEvalErrors(), 6)

    def test_create_bc_returns_tpzbndcondt(self):
        mat = pyneopz.TPZElasticity2D(1)
        bc = mat.CreateBC(-1, 0, [[0.0, 0.0], [0.0, 0.0]], [0.0, 0.0])
        self.assertIsInstance(bc, pyneopz.TPZBndCondT)


class TestElasticity3D(unittest.TestCase):
    def test_set_material_data_hook(self):
        # Kept as SetMaterialDataHook (not renamed to SetElasticity like the
        # rest of the family) to match the C++ name exactly.
        mat = pyneopz.TPZElasticity3D(1)
        mat.SetMaterialDataHook(1.0, 0.3)
        self.assertFalse(mat.HasExactSol())
        example = pyneopz.TElasticity3DAnalytic()
        mat.SetExactSol(example, 2)
        self.assertTrue(mat.HasExactSol())

    def test_n_eval_errors(self):
        mat = pyneopz.TPZElasticity3D(1)
        self.assertEqual(mat.NEvalErrors(), 3)

    def test_create_bc_returns_tpzbndcondt(self):
        mat = pyneopz.TPZElasticity3D(1)
        bc = mat.CreateBC(-1, 0, [[0.0, 0.0], [0.0, 0.0]], [0.0, 0.0])
        self.assertIsInstance(bc, pyneopz.TPZBndCondT)


class TestAnalyticSolution(unittest.TestCase):
    def test_name_reflects_the_selected_case(self):
        example = pyneopz.TLaplaceExample1()
        example.fExact = pyneopz.TLaplaceExample1.EExactSol.ESinSin
        self.assertEqual(example.Name(), "SinSin")

    def test_linear_solution_has_no_source(self):
        # u = x has zero Laplacian: the consistent source term is zero
        example = pyneopz.TLaplaceExample1()
        example.fExact = pyneopz.TLaplaceExample1.EExactSol.EX
        u, _ = example.Solution([0.5, 0.3, 0.0])
        force = example.Force([0.5, 0.3, 0.0])
        self.assertAlmostEqual(u[0], 0.5, places=10)
        self.assertAlmostEqual(force[0], 0.0, places=10)

    def test_sinsin_solution_matches_its_laplacian(self):
        # u = sin(pi x) sin(pi y); force = -Laplacian(u) = 2 pi^2 sin(pi x) sin(pi y)
        example = pyneopz.TLaplaceExample1()
        example.fExact = pyneopz.TLaplaceExample1.EExactSol.ESinSin
        u, _ = example.Solution([0.5, 0.5, 0.0])
        force = example.Force([0.5, 0.5, 0.0])
        self.assertAlmostEqual(u[0], 1.0, places=10)
        self.assertAlmostEqual(force[0], 2.0 * math.pi**2, places=6)

    def test_elasticity_2d_stretch_x(self):
        # u = (x, 0): constant strain, consistent source term is zero
        pyneopz.TElasticity2DAnalytic.gE = 1.0
        pyneopz.TElasticity2DAnalytic.gPoisson = 0.0
        example = pyneopz.TElasticity2DAnalytic()
        example.fProblemType = pyneopz.TElasticity2DAnalytic.EDefState.EStretchx
        self.assertEqual(example.Name(), "StretchingX")
        u, _ = example.Solution([0.5, 0.3, 0.0])
        force = example.Force([0.5, 0.3, 0.0])
        self.assertAlmostEqual(u[0], 0.5, places=10)
        self.assertAlmostEqual(u[1], 0.0, places=10)
        self.assertAlmostEqual(force[0], 0.0, places=9)
        self.assertAlmostEqual(force[1], 0.0, places=9)

    def test_elasticity_3d_construction(self):
        example = pyneopz.TElasticity3DAnalytic()
        example.fProblemType = pyneopz.TElasticity3DAnalytic.EDefState.EDispx
        example.fE = 1.0
        example.fPoisson = 0.0
        u, _ = example.Solution([0.5, 0.3, 0.2])
        self.assertEqual(len(u), 3)


class TestHDivApproxSpaceCreator(unittest.TestCase):
    @staticmethod
    def MakeCreator(gmesh, pOrder=1):
        darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
        darcy.SetConstantPermeability(1.0)
        bc = darcy.CreateBC(-1, 0, [[0.0]], [1.0])
        creator = pyneopz.TPZHDivApproxCreator(gmesh)
        creator.SetProbType(pyneopz.ProblemType.EDarcy)
        creator.SetDefaultOrder(pOrder)
        creator.InsertMaterialObject(darcy)
        creator.InsertMaterialObject(bc)
        return creator

    def test_default_constructor(self):
        creator = pyneopz.TPZHDivApproxCreator()
        self.assertIsInstance(creator, pyneopz.TPZApproxCreator)

    def test_atomic_meshes_returned_with_the_lagrange_counter(self):
        creator = self.MakeCreator(Create2DGeoMesh())
        meshes, lagLevelCounter = creator.CreateAtomicMeshes()
        self.assertEqual(len(meshes), 2)
        self.assertEqual(lagLevelCounter, 2)
        for mesh in meshes:
            self.assertIsInstance(mesh, pyneopz.TPZCompMesh)

    def test_step_by_step_matches_the_driver(self):
        # CreateAtomicMeshes plus CreateMultiPhysicsMesh must reproduce what
        # CreateApproximationSpace does in a single call
        creator = self.MakeCreator(Create2DGeoMesh())
        meshes, lagLevelCounter = creator.CreateAtomicMeshes()
        cmesh = creator.CreateMultiPhysicsMesh(meshes, lagLevelCounter)
        self.assertIsInstance(cmesh, pyneopz.TPZMultiphysicsCompMesh)
        self.assertEqual(cmesh.NEquations(), 40)

    def test_inherits_from_abstract_base(self):
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        self.assertIsInstance(creator, pyneopz.TPZApproxCreator)

    def test_virtual_dispatch_through_base(self):
        # SetProbType is pure virtual: bound on the base, it runs the override
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        creator.SetProbType(pyneopz.ProblemType.EDarcy)
        self.assertEqual(creator.ProbType(), pyneopz.ProblemType.EDarcy)

    def test_hdiv_family_can_be_read_and_set(self):
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        creator.HdivFamily(pyneopz.HDivFamily.EHDivConstant)
        self.assertEqual(creator.HdivFamily(), pyneopz.HDivFamily.EHDivConstant)

    def test_approximation_space(self):
        cmesh = CreateHDivApproximationSpace(Create2DGeoMesh())
        self.assertEqual(cmesh.NEquations(), 40)
        self.assertIsInstance(cmesh, pyneopz.TPZCompMesh)

    def test_invalid_setup_raises(self):
        # with no material inserted, the NeoPZ DebugStop arrives as RuntimeError
        creator = pyneopz.TPZHDivApproxCreator(Create2DGeoMesh())
        creator.SetProbType(pyneopz.ProblemType.EDarcy)
        with self.assertRaises(RuntimeError):
            creator.CreateApproximationSpace()


class TestApproxCreator(unittest.TestCase):
    """Members of the abstract base, exercised through a concrete creator."""

    def test_rigid_body_spaces_can_be_read_and_set(self):
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        self.assertFalse(creator.IsRigidBodySpaces())
        creator.IsRigidBodySpaces(True)
        self.assertTrue(creator.IsRigidBodySpaces())

    def test_should_condense_can_be_read(self):
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        # GetShouldCondense ignores its argument in NeoPZ itself; any bool works.
        self.assertTrue(creator.GetShouldCondense(False))
        creator.SetShouldCondense(False)
        self.assertFalse(creator.GetShouldCondense(False))

    def test_hybridization_data_is_nested(self):
        self.assertTrue(hasattr(pyneopz.TPZApproxCreator, "HybridizationData"))

    def test_hybridization_data_is_filled_by_set_hybrid_type(self):
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        creator.SetProbType(pyneopz.ProblemType.EDarcy)
        creator.SetHybridType(pyneopz.HybridizationType.EStandard)
        data = creator.HybridData()
        self.assertIsInstance(data, pyneopz.TPZApproxCreator.HybridizationData)
        self.assertEqual(data.fMultipliers, [1.0, 1.0, 1.0, -1.0])
        self.assertIsInstance(data.fInterfaces, dict)

    def test_hybridization_data_round_trip(self):
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        data = creator.HybridData()
        data.fWrapMatId = 999
        creator.SetHybridData(data)
        self.assertEqual(creator.HybridData().fWrapMatId, 999)

    def test_tpzhybrid_fields(self):
        hybrid = pyneopz.TPZHybrid(1, 2, 3)
        self.assertEqual((hybrid.fLagrange, hybrid.fLeft, hybrid.fRight), (1, 2, 3))


class TestH1ApproxSpaceCreator(unittest.TestCase):
    def test_default_constructor(self):
        # unlike TPZCompMesh, the default constructor here is public
        creator = pyneopz.TPZH1ApproxCreator()
        self.assertIsInstance(creator, pyneopz.TPZApproxCreator)

    def test_classic_h1_space(self):
        # one degree of freedom per node for order 1 on a 3x3 node grid
        cmesh = CreateH1ApproximationSpace(Create2DGeoMesh())
        self.assertEqual(cmesh.NEquations(), 9)

    def test_atomic_meshes_are_returned_as_a_list(self):
        # hybridization is required, and it demands an extra internal order
        darcy = pyneopz.TPZDarcyFlow(1, 2)
        darcy.SetConstantPermeability(1.0)
        bc = darcy.CreateBC(-1, 0, [[0.0]], [1.0])

        creator = pyneopz.TPZH1ApproxCreator(Create2DGeoMesh())
        creator.SetProbType(pyneopz.ProblemType.EDarcy)
        creator.SetHybridType(pyneopz.HybridizationType.EStandard)
        creator.SetDefaultOrder(1)
        creator.SetExtraInternalOrder(2)
        creator.InsertMaterialObject(darcy)
        creator.InsertMaterialObject(bc)

        meshes = creator.CreateAtomicMeshes()
        self.assertIsInstance(meshes, list)
        self.assertEqual(len(meshes), 2)
        for mesh in meshes:
            self.assertIsInstance(mesh, pyneopz.TPZCompMesh)

    def test_constant_pressure_solution(self):
        cmesh = CreateH1ApproximationSpace(Create2DGeoMesh())
        an = pyneopz.TPZLinearAnalysis(cmesh)
        an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
        solver = pyneopz.TPZStepSolver()
        # a classic H1 Darcy system is symmetric positive definite
        solver.SetDirect(pyneopz.DecomposeType.ECholesky)
        an.SetSolver(solver)
        an.Assemble()
        an.Solve()
        # sqrt(9) for the nine nodal pressures equal to one
        self.assertAlmostEqual(an.SolutionNorm(), 3.0, places=9)


class TestMHMHDivApproxSpaceCreator(unittest.TestCase):
    @staticmethod
    def MakeCreator(gmesh, elementPartition=None):
        darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
        darcy.SetConstantPermeability(1.0)
        bc = darcy.CreateBC(-1, 0, [[0.0]], [1.0])
        if elementPartition is None:
            creator = pyneopz.TPZMHMHDivApproxCreator(gmesh)
        else:
            creator = pyneopz.TPZMHMHDivApproxCreator(gmesh, elementPartition)
        creator.SetProbType(pyneopz.ProblemType.EDarcy)
        creator.SetDefaultOrder(1)
        creator.SetPOrderSkeleton(1)
        creator.InsertMaterialObject(darcy)
        creator.InsertMaterialObject(bc)
        return creator

    def test_multiple_inheritance(self):
        creator = pyneopz.TPZMHMHDivApproxCreator(pyneopz.TPZGeoMesh())
        self.assertIsInstance(creator, pyneopz.TPZHDivApproxCreator)
        self.assertIsInstance(creator, pyneopz.TPZMHMApproxCreator)
        self.assertIsInstance(creator, pyneopz.TPZApproxCreator)

    def test_inherited_members(self):
        creator = pyneopz.TPZMHMHDivApproxCreator(pyneopz.TPZGeoMesh())
        # from TPZHDivApproxCreator
        self.assertEqual(creator.HdivFamily(), pyneopz.HDivFamily.EHDivStandard)
        # from TPZApproxCreator
        creator.SetNState(1)
        self.assertEqual(creator.NState(), 1)

    def test_skeleton_order(self):
        creator = pyneopz.TPZMHMHDivApproxCreator(pyneopz.TPZGeoMesh())
        creator.SetPOrderSkeleton(2)
        self.assertEqual(creator.GetPOrderSkeleton(), 2)

    def test_inherited_driver(self):
        creator = self.MakeCreator(Create2DGeoMesh())
        cmesh = creator.CreateApproximationSpace()
        self.assertEqual(cmesh.NEquations(), 40)

    def test_mhm_driver_with_element_partition(self):
        creator = self.MakeCreator(Create2DGeoMesh(), [0, 0, 1, 1])
        cmesh = creator.BuildMultiphysicsCMesh()
        self.assertIsInstance(cmesh, pyneopz.TPZMultiphysicsCompMesh)
        self.assertEqual(cmesh.NEquations(), 64)


class TestLinearAnalysis(unittest.TestCase):
    def test_constant_pressure_solution(self):
        cmesh = CreateHDivApproximationSpace(Create2DGeoMesh())
        an = pyneopz.TPZLinearAnalysis(cmesh)
        an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
        solver = pyneopz.TPZStepSolver()
        # mixed formulation: symmetric indefinite system, LDLt and not Cholesky
        solver.SetDirect(pyneopz.DecomposeType.ELDLt)
        an.SetSolver(solver)
        an.Assemble()
        an.Solve()
        self.assertAlmostEqual(an.SolutionNorm(), 4.0, places=9)


class TestIterativeSolvers(unittest.TestCase):
    def test_gmres_matches_direct_on_the_mixed_formulation(self):
        # symmetric indefinite saddle point: only GMRES is guaranteed to work
        gmesh = Create2DGeoMesh()
        pyneopz.TPZCheckGeom(gmesh).UniformRefine(1)
        direct = SolveDirect(CreateHDivApproximationSpace(gmesh))

        gmesh2 = Create2DGeoMesh()
        pyneopz.TPZCheckGeom(gmesh2).UniformRefine(1)
        iterative = SolveIteratively(CreateHDivApproximationSpace(gmesh2), useCG=False)

        self.assertAlmostEqual(direct.SolutionNorm(), iterative.SolutionNorm(), places=6)

    def test_cg_matches_direct_on_an_spd_system(self):
        # classic (non-mixed) H1 Darcy: genuinely SPD, unlike the formulation above
        cholesky = pyneopz.DecomposeType.ECholesky
        direct = SolveDirect(CreateH1ApproximationSpace(Create2DGeoMesh(3), 2), cholesky)
        iterative = SolveIteratively(
            CreateH1ApproximationSpace(Create2DGeoMesh(3), 2), useCG=True, decomposeType=cholesky
        )
        self.assertAlmostEqual(direct.SolutionNorm(), iterative.SolutionNorm(), places=6)

    def test_set_num_threads_does_not_crash(self):
        cmesh = CreateHDivApproximationSpace(Create2DGeoMesh())
        strmat = pyneopz.TPZSkylineStructMatrix(cmesh)
        strmat.SetNumThreads(2)
        an = pyneopz.TPZLinearAnalysis(cmesh)
        an.SetStructuralMatrix(strmat)
        solver = pyneopz.TPZStepSolver()
        solver.SetDirect(pyneopz.DecomposeType.ELDLt)
        an.SetSolver(solver)
        an.Assemble()
        an.Solve()
        self.assertAlmostEqual(an.SolutionNorm(), 4.0, places=9)


class TestNumpyInterop(unittest.TestCase):
    @staticmethod
    def SolveConstantPressureProblem():
        cmesh = CreateHDivApproximationSpace(Create2DGeoMesh())
        an = pyneopz.TPZLinearAnalysis(cmesh)
        an.SetStructuralMatrix(pyneopz.TPZSkylineStructMatrix(cmesh))
        solver = pyneopz.TPZStepSolver()
        solver.SetDirect(pyneopz.DecomposeType.ELDLt)
        an.SetSolver(solver)
        an.Assemble()
        an.Solve()
        return an

    def test_solution_is_zero_copy(self):
        an = self.SolveConstantPressureProblem()
        sol = an.Solution()
        self.assertEqual(sol.shape, (40, 1))
        # same buffer read through two independent paths: C++ Norm() inside
        # SolutionNorm, and numpy reading the array it was handed
        self.assertAlmostEqual(float(np.linalg.norm(sol)), an.SolutionNorm(), places=9)
        self.assertTrue(np.shares_memory(sol, np.asarray(sol)))

    def test_solution_mutation_is_visible_through_the_same_buffer(self):
        an = self.SolveConstantPressureProblem()
        sol = an.Solution()
        sol[0, 0] = 999.0
        self.assertEqual(an.Solution()[0, 0], 999.0)

    def test_solution_survives_the_analysis_going_out_of_scope(self):
        # base=self on the returned array keeps the TPZLinearAnalysis (and the
        # buffer it owns) alive even after the last Python reference to it,
        # other than the array itself, is gone
        sol = self.SolveConstantPressureProblem().Solution()
        gc.collect()
        self.assertAlmostEqual(float(np.linalg.norm(sol)), 4.0, places=9)


class TestErrorEvaluation(unittest.TestCase):
    def test_has_exact_sol_reflects_state(self):
        darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
        self.assertFalse(darcy.HasExactSol())
        example = pyneopz.TLaplaceExample1()
        example.fExact = pyneopz.TLaplaceExample1.EExactSol.EConst
        darcy.SetExactSol(example, 2)
        self.assertTrue(darcy.HasExactSol())

    def test_error_is_near_zero_for_an_analytic_solution(self):
        # u = 1 is exactly representable: every error component should be
        # at machine precision, not just "small"
        darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
        darcy.SetConstantPermeability(1.0)
        example = pyneopz.TLaplaceExample1()
        example.fExact = pyneopz.TLaplaceExample1.EExactSol.EConst
        darcy.SetExactSol(example, 2)

        self.assertEqual(darcy.NEvalErrors(), 5)
        for component in EvaluateConstantSolutionError(darcy):
            self.assertLess(component, 1e-9)

    def test_error_is_near_zero_for_a_python_callback(self):
        # same problem as above, but the exact solution is a Python callable
        # instead of a TLaplaceExample1; both paths should agree
        darcy = pyneopz.TPZMixedDarcyFlow(1, 2)
        darcy.SetConstantPermeability(1.0)
        darcy.SetExactSol(lambda loc: ([1.0], [[0.0], [0.0], [0.0]]), 2)

        for component in EvaluateConstantSolutionError(darcy):
            self.assertLess(component, 1e-9)

    def test_elasticity_error_is_near_zero_for_an_analytic_solution(self):
        # errors[6] is the exact solution's own energy norm (strain*stress),
        # not a discrepancy: for EStretchx with E=1, nu=0, strain=stress=
        # [[1,0],[0,0]], giving exactly 1. Only errors[0:6] must vanish.
        pyneopz.TElasticity2DAnalytic.gE = 1.0
        pyneopz.TElasticity2DAnalytic.gPoisson = 0.0
        example = pyneopz.TElasticity2DAnalytic()
        example.fProblemType = pyneopz.TElasticity2DAnalytic.EDefState.EStretchx

        def setExact(mat, bc):
            mat.SetExactSol(example, 2)
            bc.SetForcingFunctionBC(example, 2)

        nErrors, error = EvaluateStretchXError(setExact)
        self.assertEqual(nErrors, 7)
        for component in error[:6]:
            self.assertLess(component, 1e-9)
        self.assertAlmostEqual(error[6], 1.0, places=9)

    def test_elasticity_error_is_near_zero_for_a_python_callback(self):
        # same problem via a Python callable forwarding to
        # TElasticity2DAnalytic.Solution; both paths must agree, digit for digit
        pyneopz.TElasticity2DAnalytic.gE = 1.0
        pyneopz.TElasticity2DAnalytic.gPoisson = 0.0
        example = pyneopz.TElasticity2DAnalytic()
        example.fProblemType = pyneopz.TElasticity2DAnalytic.EDefState.EStretchx

        def setExact(mat, bc):
            mat.SetExactSol(lambda loc: example.Solution(loc), 2)
            bc.SetForcingFunctionBC(lambda loc: example.Solution(loc), 2)

        _, error = EvaluateStretchXError(setExact)
        for component in error[:6]:
            self.assertLess(component, 1e-9)
        self.assertAlmostEqual(error[6], 1.0, places=9)


if __name__ == "__main__":
    unittest.main()
