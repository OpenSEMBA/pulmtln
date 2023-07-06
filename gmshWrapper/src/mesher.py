import gmsh
import os
import sys
from collections import defaultdict

DEFAULT_MESHING_OPTIONS = {
    "Mesh.MeshSizeFromCurvature": 10,
    "Mesh.ElementOrder": 3,
    "Mesh.ScalingFactor": 1e-3,
    "General.DrawBoundingBoxes": 1,
    "Mesh.SurfaceFaces": 1
}

class StepShapes:
    def __init__(self, shapes):
        gmsh.model.occ.synchronize()

        self.allShapes = shapes

        self.pecs = self. get_surfaces(shapes, "Conductor_")
        self.dielectrics = self.get_surfaces(shapes, "Dielectric_")

    @staticmethod
    def get_surfaces(shapes, label: str):
        surfaces = dict()
        for s in shapes:
            entity_name = gmsh.model.get_entity_name(*s)
            if s[0] != 2 or label not in entity_name:
                continue
            ini = entity_name.index(label) + len(label)
            end = ini+1 # BUGGY FOR TWO DIGIT NUMBERS !!!!
            num = int(entity_name[ini:end])
            surfaces[num] = s

        return surfaces

def meshFromStep(
        folder: str,
        case_name: str,
        meshing_options=DEFAULT_MESHING_OPTIONS):

    gmsh.initialize()
    gmsh.model.add(case_name)

    # Importing from FreeCAD generated steps. STEP default units are mm.
    stepShapes = StepShapes(
        gmsh.model.occ.importShapes(
            folder + case_name + '/' + case_name + '.step',
            highestDimOnly=False
        )
    )
    
    # --- Geometry manipulation ---
    # Creates global domain.
    region = stepShapes.pecs[0] 

    for num, surf in stepShapes.pecs.items():
        if num == 0:
            continue
        for _, dielectric_surf in stepShapes.dielectrics.items():
            gmsh.model.occ.cut([dielectric_surf], [surf], removeTool=False)    
        gmsh.model.occ.cut([region], [surf], removeTool=False)

    for _, surf in stepShapes.dielectrics.items():
        gmsh.model.occ.cut([region], [surf], removeTool=False)
    
    gmsh.model.occ.synchronize()

    # Prepares PEC boundaries and removes surfaces.
    pec_bdrs = dict()
    for num, surf in stepShapes.pecs.items():
        pec_bdrs[num] = gmsh.model.getBoundary([surf])
        if num != 0:
            gmsh.model.occ.remove([surf])
    
    
            
    gmsh.model.occ.synchronize()

    # --- Physical groups ---
    # Boundaries.
    for num, bdrs in pec_bdrs.items():
        name = "Conductor_" + str(num)
        for bdr in bdrs:
            if bdr[1] > 0:
                gmsh.model.addPhysicalGroup(1, [bdr[1]], name=name)

    # Domains.
    gmsh.model.addPhysicalGroup(2, [region[1]], name='Vacuum')

    for num, surf in stepShapes.dielectrics.items():
        name = "Dielectric_" + str(num)
        gmsh.model.addPhysicalGroup(2, [surf[1]], name=name)

    # Meshing.
    for [opt, val] in meshing_options.items():
        gmsh.option.setNumber(opt, val)

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    # Exporting
    gmsh.write(case_name + '.msh')

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()
