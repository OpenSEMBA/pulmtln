import gmsh
import os
import sys
from collections import defaultdict

CASE_NAME = 'empty_coax'

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'

gmsh.initialize()
gmsh.model.add(CASE_NAME)

# Importing from FreeCAD generated steps. STEP default units are mm.
allShapes = gmsh.model.occ.importShapes(
    dir_path + 'empty_coax.step', 
    highestDimOnly=False)
gmsh.model.occ.synchronize()

# Build map of names to entities.
conductors_surface = dict()
conductors_bdr = dict()
for s in allShapes:
    entity_name = gmsh.model.get_entity_name(*s)
    label = "Conductor_"
    if label in entity_name:
        ini = entity_name.index(label) + len(label)
        end = ini+1
        num = int(entity_name[ini:end])
        if s[0] == 2:
            conductors_surface[num] = s
        elif s[0] == 1:
            conductors_bdr[num] = s
        else:
            raise ValueError("Invalid shape dimension for conductor.")


# --- Geometry manipulation ---
# Creates domain.
region = conductors_surface[0]
for i in range(1,len(conductors_surface)):
    gmsh.model.occ.cut(
        [region], 
        [conductors_surface[i]], removeTool=True)
gmsh.model.occ.synchronize()

# --- Physical groups ---
# Boundaries.
for conductor_num, bdr in conductors_bdr.items():
    name = "Conductor_" + str(conductor_num)
    tag_id = gmsh.model.addPhysicalGroup(1, [bdr[1]], name=name)

# Domains.
domain_tag = { 'Vacuum': 1 }
gmsh.model.addPhysicalGroup(2, [region[1]], name='Vacuum')

# Meshing.
gmsh.option.setNumber("Mesh.MeshSizeMin", 10)
gmsh.option.setNumber("Mesh.MeshSizeMax", 10)
gmsh.option.setNumber("Mesh.ElementOrder", 3)
gmsh.option.setNumber("Mesh.ScalingFactor", 1e-3)
gmsh.model.mesh.setSize([(0,5)], 0.01)
gmsh.model.mesh.generate(2)

# Exporting
# gmsh.write(dir_path + CASE_NAME + '.vtk')

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(dir_path + CASE_NAME + '.msh')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
