import gmsh
import os
import sys
from collections import defaultdict

def meshFromStep(dir_path: str, case_name: str):
    gmsh.initialize()
    gmsh.model.add(case_name)

    # Importing from FreeCAD generated steps. STEP default units are mm.
    allShapes = gmsh.model.occ.importShapes(
        dir_path + case_name + '/' + case_name + '.step', 
        highestDimOnly=False)
    gmsh.model.occ.synchronize()

    # Build map of names to entities.
    conductors_surface = dict()
    conductors_bdr = dict()
    for s in allShapes:
        entity_name = gmsh.model.get_entity_name(*s)
        label = "Conductor_"
        if label in entity_name and s[0] == 2:
            ini = entity_name.index(label) + len(label)
            end = ini+1
            num = int(entity_name[ini:end])
            conductors_surface[num] = s
            conductors_bdr[num] = gmsh.model.get_boundary([s])


    # --- Geometry manipulation ---
    # Creates domain.
    region = conductors_surface[0]
    for k, v in conductors_surface.items():
        if k == 0:
            continue
        gmsh.model.occ.cut([region], [v], removeTool=True)
    gmsh.model.occ.synchronize()

    # --- Physical groups ---
    # Boundaries.
    for conductor_num, bdrs in conductors_bdr.items():
        name = "Conductor_" + str(conductor_num)
        for bdr in bdrs:
            gmsh.model.addPhysicalGroup(1, [bdr[1]], name=name)

    # Domains.
    domain_tag = { 'Vacuum': 1 }
    gmsh.model.addPhysicalGroup(2, [region[1]], name='Vacuum')

    # Meshing.
    gmsh.option.setNumber("Mesh.MeshSizeMin",   0.1)
    gmsh.option.setNumber("Mesh.MeshSizeMax",   10)
    gmsh.option.setNumber("Mesh.ElementOrder",  3)
    gmsh.option.setNumber("Mesh.ScalingFactor", 1e-3)
    gmsh.model.mesh.setSize([(0,5)], 0.01)
    gmsh.model.mesh.generate(2)

    # Exporting
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(case_name + '.msh')

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()
