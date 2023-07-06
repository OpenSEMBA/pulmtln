from src.mesher import meshFromStep

import os

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
testdata_path = dir_path + '/../testData/'

def test_partially_filled_coax():
    meshFromStep(testdata_path, 'partially_filled_coax')

def test_empty_coax():
    meshFromStep(testdata_path, 'empty_coax')

def test_two_wires_coax():
    meshFromStep(testdata_path, 'two_wires_coax')