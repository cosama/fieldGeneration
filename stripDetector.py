from fenics import *
from dolfin import *
from mshr import *
from math import *
import numpy as np

def BuildPlanarMesh(x, y, thickness):
    mesh = Mesh()
    # Define 3D geometry - for now, just a box
    box = Box(Point(0, 0, 0), Point(x, y, thickness))

    print("\nBuilding mesh for planar geometry ({} mm thick).".format(thickness))
    generator = CSGCGALMeshGenerator3D()
    mesh = generator.generate(CSGCGALDomain3D(box))
    
    # Refine mesh at the surface to get the electrodes better...
    cell_markers = MeshFunction("bool", mesh, 0)
    cell_markers.set_all(False)
    surface = CompiledSubDomain("on_boundary")
    surface.mark(cell_markers, True)
    mesh = refine(mesh, cell_markers)
    mesh = refine(mesh, cell_markers)
    return mesh
