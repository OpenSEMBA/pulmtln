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
        self.allShapes = shapes

        self.pec_surfaces = self. get_surfaces(shapes, "Conductor_")
        self.dielectric_surface = self.get_surfaces(shapes, "Dielectric_")
        
        self.pec_bdrs = self.get_boundaries(self.pec_surfaces)
        self.dielectric_bdr = self.get_boundaries(self.dielectric_surface)


    def get_surfaces(shapes, label: str):
        surfaces = dict()
        for s in shapes:
            entity_name = gmsh.model.get_entity_name(*s)
            if s[0] != 2:
                continue
            ini = entity_name.index(label) + len(label)
            end = ini+1
            num = int(entity_name[ini:end])
            surfaces[num] = s
        
        return surfaces
    
    
    def get_boundaries(surfaces):
        boundaries = dict()
        for [num, s] in surfaces.items():
            boundaries[num] = gmsh.model.get_boundary([s])

        return boundaries


def meshFromStep(
        folder: str, 
        case_name: str, 
        meshing_options = DEFAULT_MESHING_OPTIONS):
    
    gmsh.initialize()
    gmsh.model.add(case_name)

    # Importing from FreeCAD generated steps. STEP default units are mm.
    stepShapes = StepShapes(
        gmsh.model.occ.importShapes(
        folder + case_name + '/' + case_name + '.step', highestDimOnly=False)
    )
    gmsh.model.occ.synchronize()

    # --- Geometry manipulation ---
    # Creates global domain.
    region = stepShapes.conductors_surface[0]
    for k, v in stepShapes.conductors_surface.items():
        if k == 0:
            continue
        gmsh.model.occ.cut([region], [v], removeTool=True)
    gmsh.model.occ.synchronize()

    # --- Physical groups ---
    # Boundaries.
    for conductor_num, bdrs in stepShapes.conductors_bdr.items():
        name = "Conductor_" + str(conductor_num)
        for bdr in bdrs:
            gmsh.model.addPhysicalGroup(1, [bdr[1]], name=name)

    # Domains.
    domain_tag = { 'Vacuum': 1 }
    gmsh.model.addPhysicalGroup(2, [region[1]], name='Vacuum')

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
