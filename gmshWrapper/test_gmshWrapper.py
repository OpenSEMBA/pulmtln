from src.mesher import *

import os
import gmsh

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
testdata_path = dir_path + '/../testData/'


def test_partially_filled_coax():
    meshFromStep(testdata_path, 'partially_filled_coax')


def test_empty_coax():
    meshFromStep(testdata_path, 'empty_coax')


def test_two_wires_coax():
    meshFromStep(testdata_path, 'two_wires_coax')


def test_stepShapes_for_partially_filled_coax():
    case_name = 'partially_filled_coax'

    gmsh.initialize()
    gmsh.model.add(case_name)

    stepShapes = StepShapes(
        gmsh.model.occ.importShapes(
            dir_path + case_name + '/' + case_name + '.step', highestDimOnly=False
        )
    )

    gmsh.finalize()

    assert (len(stepShapes.conductors_bdr) == 2)
    assert (len(stepShapes.conductors_surface) == 2)
    assert (len(stepShapes.dielectric_surface) == 1)
