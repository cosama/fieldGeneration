from fenics import *
from dolfin import *
from mshr import *
from math import *
import numpy as np

from gretinaGeometry import *
from agataGeometry import *

class OuterContact(SubDomain):
    def inside(self, x, on_boundary):
        r = sqrt(x[0]*x[0] + x[1]*x[1])
        return on_boundary and x[2] < 90.

class InnerContact(SubDomain):
    def inside(self, x, on_boundary):
        r = sqrt(x[0]*x[0] + x[1]*x[1])
        return on_boundary and r < 7. and x[2] > 0.1

class TopContact(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2],0)

class BottomContact(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], 5.)
    
def Solver(mesh, V, f, uOut, uCore, cOut, cCore):
    epsilon0 = Constant(8.8541878E-15)
    epsilonGe = 16*epsilon0

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(epsilonGe*grad(u), grad(v))*dx
    L = f*v*dx

    bc_Out = DirichletBC(V, uOut, cOut)
    bc_Core = DirichletBC(V, uCore, cCore)
    bcs = [bc_Out, bc_Core]
    
    # Compute solution
    u = Function(V)
    problem = LinearVariationalProblem(a, L, u, bcs)
    solver = LinearVariationalSolver(problem)
    solver.parameters["linear_solver"] = "mumps"
    solver.solve()
    return u

def BuildPlanar(thickness):
    mesh = Mesh()
    # Define 3D geometry
    box = Box(Point(0, 0, 0), Point(80., 80., thickness))

    print("\nBuilding mesh for simple planar ({} mm thick).".format(thickness))
    generator = CSGCGALMeshGenerator3D()
    generator.parameters["edge_size"] = 0.05
    generator.parameters["facet_angle"] = 25.0
    generator.parameters["facet_size"] = 0.05
    
    mesh = generator.generate(CSGCGALDomain3D(box))
    return mesh


def SolveField(fieldType, detType, xtalType=None, xtalHV=None, rho0=None, drho=None, xtalDepletion=None):

    # Define geometry
    if (detType == "planar"):
        mesh = BuildPlanar(5.)

    if (detType == "coaxG"):
        xtaltype = 'A'
        mesh = BuildGRETINAMesh(xtaltype)

    if (detType == "coaxA"):
        xtaltype = 'A'
        mesh = BuildAGATAMesh(xtaltype)


    sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    sub_domains.set_all(0)

    if (detType == "coaxG" or detType == "coaxA"):
        outerContact = OuterContact()
        outerContact.mark(sub_domains,1)
        innerContact = InnerContact()
        innerContact.mark(sub_domains, 7)

    if (detType == "planar"):
        outerContact = BottomContact()
        outerContact.mark(sub_domains, 1)
        innerContact = TopContact()
        innerContact.mark(sub_domains, 7)

    meshfile = File('mesh.pvd')
    meshfile << sub_domains
        
    # Defining function space and constants
    V = FunctionSpace(mesh, 'CG', 1)

    e = 1.602176E-19
    ePerCubic_mm = e*1.0E10/1000. # Assume impurity in units of 10^10 cm^3

    print("Mesh initialized.  Beginning calculations...")

    # Scale impurities to reproduce depletion voltage
    if (fieldType == "impurity" or fieldType == "all"):
        fScale = 1.0
        iScale = 0.5

        if xtalDepletion is None:
            xtalDepletion = 2000.
        if rho0 is None:
            rho0 = 0.0
        if drho is None:
            drho = 0.0
    
        while (iScale > 0.0001):
            irho0 = rho0*fScale
            idrho = drho*fScale
    
            u_Out = Constant(0.0)
            u_Core = Constant(xtalDepletion)
    
            f = Expression("{}*({} + {}*x[2])".format(ePerCubic_mm, irho0, idrho*0.1), degree=1)
            u = Solver(mesh, V, f, u_Out, u_Core, outerContact, innerContact)
            vertex_values_u = u.compute_vertex_values()
            maxV = np.amax(vertex_values_u)
            print("f={} -- max = {}, depleted = {}".format(fScale, maxV, near(maxV, xtalDepletion, 1E-8)))
            if (near(maxV, xtalDepletion, 1E-8)):
                fScale = fScale + iScale
            if not near(maxV, xtalDepletion, 1E-8):
                fScale = fScale - iScale
            iScale = iScale/2.
            
        print("\nTo match depletion voltage ({} V), impurities scaled by {} -- rho0 = {}, drho = {}.".format(xtalDepletion, fScale, irho0, idrho))
        rho0 = irho0
        drho = idrho
        
    if (fieldType == "eField" or fieldType == "all"):
        
        # Calculate electric potential first
        if xtalHV is None:
            xtalHV = 4000.
        if rho0 is None:
            rho0 = 0.0
        if drho is None:
            drho = 0.0

        u_Out = Constant(0.0)
        u_Core = Constant(xtalHV)

        f = Expression("{}*({} + {}*x[2])".format(ePerCubic_mm, irho0, idrho*0.1), degree=1)
        u = Solver(mesh, V, f, u_Out, u_Core, outerContact, innerContact)
                
        potentialfile = File('potential.pvd')
        potentialfile << u

        fieldVec = VectorFunctionSpace(mesh, "CG", 1)
        gradu = project(grad(u), fieldVec)
    
        fieldFile = File('field.pvd')
        fieldFile << gradu  

    if (fieldType == "wPot" or fieldType == "all"):

        # Now solve for weighting potentials (all segments + CC)
        u_Ground = Constant(0.0)
        u_Unity = Constant(1.0)

        f = Constant(0.0)

        # Segments only in properly for GRETINA right now, so only Gretina for WP for the moment
        if (detType == "coaxG"):
            for segment in range(36):
                print("Solving for weighting potential (segment {}).".format(segment))
                if (segment<36):
                    u = Solver(mesh, V, f, GetGRETINAOuterWP(xtaltype, segment), u_Ground, outerContact, innerContact)
                if (segment==36):
                    u = Solver(mesh, V, f, u_Ground, u_Unity, outerContact, innerContact)
                name = 'WP' + str(segment) + '.pvd'
                wpfile = File(name)
                wpfile << u


