from fenics import *
from dolfin import *
from mshr import *
from math import *
import numpy as np

def BuildGRETINAMesh(type):
    mesh = Mesh()
    # Define 3D geometry
    cylinder = Cylinder(dolfin.Point(0,0,0), dolfin.Point(0,0,90), 40., 40., 50);
    cc = Cylinder(Point(0,0,12), Point(0,0,90), 6., 6., 50)
    if (type == 'A'):
        quad = Surface3D("gretinaA.off")
        print("\nBuilding mesh for GRETINA geometry (type {} crystal).".format(type))
    if (type == 'B'):
        quad = Surface3D("gretinaB.off")
        print("\nBuilding mesh for GRETINA geometry (type {} crystal).".format(type))
    g3d = (cylinder * quad)-cc
    
    generator = CSGCGALMeshGenerator3D()
    generator.parameters["edge_size"] = 0.01
    generator.parameters["facet_angle"] = 25.0
    generator.parameters["facet_size"] = 0.01
    
    mesh = generator.generate(CSGCGALDomain3D(g3d))
    return mesh

def GetGRETINAOuterWP(type, seg):
    if (type == 'A'):
        theta01 = 1.50953286
        theta12 = 0.49530276
        theta23 = -0.5196095
        theta34 = -1.6157653
        theta45 = -2.6875067
        theta50 = 2.58606836
        minZ = 0.
        maxZ = 90.

        if (seg//6 == 0):
            minZ = 0.
            maxZ = 8.
        if (seg//6 == 1):
            minZ = 8.
            maxZ = 22.
        if (seg//6 == 2):
            minZ = 22.
            maxZ = 38.
        if (seg//6 == 3):
            minZ = 38.
            maxZ = 55.
        if (seg//6 == 4):
            minZ = 55.
            maxZ = 76.
        if (seg//6 == 5):
            minZ = 76.
            maxZ = 90.

            if (seg%6 == 1):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta12, ang2=theta01)
            if (seg%6 == 2):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta23, ang2=theta12)
            if (seg%6 == 3):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta34, ang2=theta23)
            if (seg%6 == 4):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta45, ang2=theta34)
            if (seg%6 == 5):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (((atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<pi)) or ((atan2(x[1],x[0])<ang2) and (atan2(x[1],x[0])>-pi)))? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta50, ang2=theta45)
            if (seg%6 == 0):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta01, ang2=theta50)

    if (type=='B'):
        theta01 = -0.0007962
        theta12 = -0.9791537
        theta23 = -2.0230593
        theta34 = -3.1120207
        theta45 = 2.13766026
        theta50 = 1.03990010
        minZ = 0.
        maxZ = 90.

        if (seg//6 == 0):
            minZ = 0.
            maxZ = 8.
        if (seg//6 == 1):
            minZ = 8.
            maxZ = 22.
        if (seg//6 == 2):
            minZ = 22.
            maxZ = 38.
        if (seg//6 == 3):
            minZ = 38.
            maxZ = 55.
        if (seg//6 == 4):
            minZ = 55.
            maxZ = 76.
        if (seg//6 == 5):
            minZ = 76.
            maxZ = 90.

            if (seg%6 == 1):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta12, ang2=theta01)
            if (seg%6 == 2):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta23, ang2=theta12)
            if (seg%6 == 3):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta34, ang2=theta23)
            if (seg%6 == 4):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (((atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<pi)) or ((atan2(x[1],x[0])<ang2) and (atan2(x[1],x[0])>-pi)))? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta45, ang2=theta34)
            if (seg%6 == 5):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta50, ang2=theta45)               
            if (seg%6 == 0):
                u_WP = Expression('(x[2]<maxZ) and (x[2]>=minZ) and (atan2(x[1],x[0])>ang1) and (atan2(x[1],x[0])<ang2) ? 1. : 0.', degree=2, maxZ=maxZ, minZ=minZ, ang1=theta01, ang2=theta50)

        
    return u_WP
