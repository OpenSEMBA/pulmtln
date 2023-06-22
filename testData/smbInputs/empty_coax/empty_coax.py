import gmsh
import os
import sys

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
conductors = dict()
for s in allShapes:
    entity_name = gmsh.model.get_entity_name(*s)
    label = "Conductor_"
    if label in entity_name:
        ini = entity_name.index(label) + len(label)
        end = ini+1
        num = int(entity_name[ini:end])
        conductors[num] = s
print(conductors)

# Geometry manipulation.
ents = gmsh.model.occ.get_entities()
res = []
for i in range(1,len(conductors)):
    cut_res = gmsh.model.occ.cut([conductors[0]], [conductors[i]], removeTool=True)
    res.append(cut_res)

gmsh.model.occ.synchronize()
ents = gmsh.model.occ.get_entities()

# Meshing.
gmsh.option.setNumber("Mesh.MeshSizeMin", 1)
gmsh.option.setNumber("Mesh.MeshSizeMax", 5)
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
