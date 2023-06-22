import gmsh
import os
import sys
from collections import defaultdict

CASE_NAME = 'empty_coax'

dir_path = os.path.dirname(os.path.realpath(__file__))

gmsh.initialize()
gmsh.model.add(CASE_NAME)

# Importing from FreeCAD generated steps.
allShapes = gmsh.model.occ.importShapes(
    dir_path + '\\' + 'empty_coax.step', 
    highestDimOnly=False)
gmsh.model.occ.synchronize()

# Build map of names to entities.
conductors_surface = dict()
conductors_boundary = dict()
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
            conductors_boundary[num] = s
        else:
            raise ValueError("Invalid shape dimension for conductor.")

# Geometry manipulation.
# Creates domain.
ents = gmsh.model.occ.get_entities()
region = conductors_surface[0]
for i in range(1,len(conductors_surface)):
    gmsh.model.occ.cut(
        [region], 
        [conductors_surface[i]], removeTool=True)

gmsh.model.occ.synchronize()
ents = gmsh.model.occ.get_entities()

# Ensures boundaries are embedded.
# embeddedFace, embFaceMap = gmsh.model.occ.fragment(region, boundary)
# gmsh.model.occ.remove(gmsh.model.occ.getEntities(1), True)
# gmsh.model.occ.synchronize()

# Physical groups.
conductor_tag = {
    0: 1,
    1: 2
}


# Creates physical groups
gmsh.model.addPhysicalGroup(2, [region[0][1]], material_tag['Vacuum'])

pec_bdr = []
for line in boundary:
    pec_bdr.append(line[1])

gmsh.model.addPhysicalGroup(1, pec_bdr, material_tag['PEC'])
gmsh.model.setPhysicalName(1, material_tag['PEC'], 'PEC')

# Meshing.
gmsh.option.setNumber("Mesh.MeshSizeMin", 5)
gmsh.option.setNumber("Mesh.MeshSizeMax", 10)
gmsh.option.setNumber("Mesh.ElementOrder", 2)
gmsh.model.mesh.setSize([(0,5)], 0.01)
gmsh.model.mesh.generate(2)



# Exporting
gmsh.write(CASE_NAME + '.vtk')

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(CASE_NAME + '.msh')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
