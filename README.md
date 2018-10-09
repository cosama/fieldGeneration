# fieldGeneration
Calculation of electric potentials, electric fields and weighting potentials for GRETINA and other HPGe detector geometries.
This code uses the FEniCS package in Python, which I've been running in a Docker container.  
To run, start the Docker container, then start Python.  I've been testing by running inside Python after loading everything from "calculateField".

from calculateField import *

SolveField(fieldType, detType, xtalType, xtalHV, rho0, drho, xtalDepletion)
>  fieldType: <"impurity"|"eField"|"wPot"|"all">  
   detType: <"coaxG"|"coaxA"|"planar"> for GRETINA, AGATA and planar respectively  
   xtalType: <"A"|"B"> for A-type or B-type  
   xtalHV:  bias voltage in V
   rho0: impurity concentration at front (higher value) in units of 10<sup>10</sup>/cm<sup>3</sup>  
   drho: impurity concentration slope, in units of 10<sup>10</sup>/cm  
   xtalDepletion: depletion voltage in V  
