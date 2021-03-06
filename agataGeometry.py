from fenics import *
from dolfin import *
from mshr import *
from math import *
import numpy as np

class OuterContact(SubDomain):
    def inside(self, x, on_boundary):
        r = sqrt(x[0]*x[0] + x[1]*x[1])
        return on_boundary and not near(x[2],90.) and (near(x[2], 0.) or r > 20.)

class InnerContact(SubDomain):
    def inside(self, x, on_boundary):
        r = sqrt(x[0]*x[0] + x[1]*x[1])
        return on_boundary and r < 7. and x[2] > 0.1

def BuildAGATAMesh(type):
    mesh = Mesh()
    # Define 3D geometry
    cylinder = Cylinder(dolfin.Point(0,0,0), dolfin.Point(0,0,90), 40., 40., 50);
    cc = Cylinder(Point(0,0,12), Point(0,0,90), 6., 6., 50)
    if (type == 'A'):
        quad = Surface3D("agataA.off")
        print("\nBuilding mesh for AGATA geometry (type {} crystal).".format(type))
        
    # Still need to make .off files for AGATA B and C type...
    if (type == 'B'):
        quad = Surface3D("agataA.off")
        print("\nBuilding mesh for AGATA geometry (type {} crystal).".format(type))
    if (type == 'C'):
        quad = Surface3D("agataA.off")
        print("\nBuilding mesh for AGATA geometry (type {} crystal).".format(type))
    g3d = (cylinder * quad)-cc
    
    generator = CSGCGALMeshGenerator3D()
    generator.parameters["edge_size"] = 0.01
    generator.parameters["facet_angle"] = 25.0
    generator.parameters["facet_size"] = 0.01
    
    mesh = generator.generate(CSGCGALDomain3D(g3d))
    return mesh

