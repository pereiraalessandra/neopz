#
# Created by Alessandra on 10/08/2026.
#
"""Regression tests for the NeoPZ Python bindings.

Uses unittest from the standard library, so no external dependency is needed.
The expected values come from manually verified runs.

To run it:

    export PYTHONPATH=<build>/python:$PYTHONPATH
    python3 -m unittest TestPyNeoPZ -v
"""

# ----- Unit test includes -----
import unittest

# ----- PZ includes -----
import pyneopz


# ==========> Helpers <==========
# ===============================

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

    def test_hdiv_family_property(self):
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        creator.HdivFamily = pyneopz.HDivFamily.EHDivConstant
        self.assertEqual(creator.HdivFamily, pyneopz.HDivFamily.EHDivConstant)

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

    def test_rigid_body_spaces_is_a_property(self):
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        self.assertFalse(creator.IsRigidBodySpaces)
        creator.IsRigidBodySpaces = True
        self.assertTrue(creator.IsRigidBodySpaces)

    def test_should_condense_can_be_read(self):
        creator = pyneopz.TPZHDivApproxCreator(pyneopz.TPZGeoMesh())
        self.assertTrue(creator.ShouldCondense())
        creator.SetShouldCondense(False)
        self.assertFalse(creator.ShouldCondense())

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
        self.assertEqual(creator.HdivFamily, pyneopz.HDivFamily.EHDivStandard)
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


if __name__ == "__main__":
    unittest.main()
