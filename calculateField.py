from fenics import *
from dolfin import *
from mshr import *
from math import *
import numpy as np

from gretinaGeometry import *
from agataGeometry import *
from stripDetector import *

parameters['allow_extrapolation'] = True
    
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

def SolveField(fieldType, detType, xtalType=None, xtalHV=None, rho0=None, drho=None, xtalDepletion=None):

    # Define geometry
    if (detType == "planar"):
        mesh = BuildPlanarMesh(20., 20., 10.)

    if (detType == "coaxG"):
        if xtalType is None:
            xtalType = 'A'
        mesh = BuildGRETINAMesh(xtalType)

    if (detType == "coaxA"):
        if xtalType is None:
            xtalType = 'A'
        mesh = BuildAGATAMesh(xtalType)

    sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    sub_domains.set_all(0)

    if (detType == "coaxG" or detType == "coaxA"):
        outerContact = CompiledSubDomain("on_boundary && !near(x[2], depth) && (near(x[2], 0.) || sqrt(x[0]*x[0]+x[1]*x[1]) > r)", depth = 90., r = 20.)
        outerContact.mark(sub_domains,2)
        innerContact = CompiledSubDomain("on_boundary && !near(x[2], 0.) && sqrt(x[0]*x[0]+x[1]*x[1]) < r", r = 7.)
        innerContact.mark(sub_domains, 7)

    if (detType == "planar"):
        outerContact = CompiledSubDomain("on_boundary && near(x[2],thickness) && fmod(x[0],pitch) < gap", thickness=5., pitch = 1., gap = 0.5)
        outerContact.mark(sub_domains, 2)
        innerContact = CompiledSubDomain("on_boundary && near(x[2], 0.) && fmod(x[1],pitch) < gap", pitch = 1., gap = 0.5)
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
        fScaleLastDepleted = 1.0
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
                fScaleLastDepleted = fScale
                fScale = fScale + iScale
            if not near(maxV, xtalDepletion, 1E-8):
                fScale = fScale - iScale
            iScale = iScale/2.
            
        rho0 = rho0*fScaleLastDepleted
        drho = drho*fScaleLastDepleted
        print("\nTo match depletion voltage ({} V), impurities scaled by {} -- rho0 = {}, drho = {}.".format(xtalDepletion, fScaleLastDepleted, irho0, idrho))
    
    if (fieldType == "eField" or fieldType == "all"):
        
        # Calculate electric potential first
        if xtalHV is None:
            xtalHV = 4000.
        if rho0 is None:
            rho0 = 0.0
        if drho is None:
            drho = 0.0

        if fieldType == "eField":
            irho0 = rho0
            idrho = drho
            
        u_Out = Constant(0.0)
        u_Core = Constant(xtalHV)

        f = Expression("{}*({} + {}*x[2])".format(ePerCubic_mm, irho0, idrho*0.1), degree=1)
        u = Solver(mesh, V, f, u_Out, u_Core, outerContact, innerContact)
                
        potentialfile = File('potential.pvd')
        potentialfile << u

        print("\nSolved potential, now getting field vector...")
        
        fieldVec = VectorFunctionSpace(mesh, "CG", 1)
        gradu = project(grad(u), fieldVec, solver_type = "mumps")

        fieldFile = File('field.pvd')
        fieldFile << gradu  

        tree = mesh.bounding_box_tree()

        fileOut = open('./fieldOutput.dat', 'w+')
        print("# x y z Ex Ey Ez |E| V", file=fileOut)
        
        for i,x in enumerate(np.linspace(-43., 43., 87.)):
            for j,y in enumerate(np.linspace(-43., 43., 87.)):
                for k,z in enumerate(np.linspace(-1., 93., 95.)):
                    if tree.compute_first_collision(Point(x,y,z)) < (2 ** 32 -1):
                        voltage = u(Point(x,y,z))
                        gradV = gradu(Point(x,y,z))
                        Emag = np.sqrt(gradV[0]*gradV[0] + gradV[1]*gradV[1] + gradV[2]*gradV[2])
                        print("{:.5f}\t{:.5f}\t{:.5f}\t{:.6E}\t{:.6E}\t{:.6E}\t{:.6E}\t{:.6E}".format(x/1000.,y/1000.,z/1000.,-gradV[0]*1000.,-gradV[1]*1000.,-gradV[2]*1000.,Emag*1000.,voltage*1000.), file=fileOut)
                    else:
                        print("{:.5f}\t{:.5f}\t{:.5f}\t{:.6E}\t{:.6E}\t{:.6E}\t{:.6E}\t{:.6E}".format(x/1000.,y/1000.,z/1000.,0.0,0.0,0.0,0.0,0.0), file=fileOut)
        
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


