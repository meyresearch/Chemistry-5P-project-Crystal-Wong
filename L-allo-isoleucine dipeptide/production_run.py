import numpy as np
from openmm.app import *
from openmm.unit import *
from openmm import *
from openmm.app.modeller import Modeller
from pdbfixer import PDBFixer
from mdtraj.reporters import DCDReporter

solvated_tmp = PDBFile('L_allo_solvated_tmp.pdb')


forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
system = forcefield.createSystem(solvated_tmp.topology,
                                 nonbondedMethod=PME,
                                 nonbondedCutoff=1*nanometer)

integrator = LangevinMiddleIntegrator(293*kelvin,
                                      1.0/picoseconds,
                                      2.0*femtoseconds)

platform = openmm.Platform.getPlatformByName('CUDA')
platform_properties = {'DeviceIndex':'0'}
simulation = Simulation(solvated_tmp.topology, system, integrator, platform=platform,
                        platformProperties=platform_properties)
simulation.loadCheckpoint('L_allo_NVT_equilibrated.chk')

# Trajectory saved every 1 ps
simulation.reporters.append(DCDReporter('L_allo_prd_equilibration.dcd', 500, np.arange(0,41,1)))
simulation.reporters.append(CheckpointReporter('L_allo_prd_equilibrated.chk', reportInterval=1000))
simulation.reporters.append(StateDataReporter('L_allo_prd_file.txt', 500, step=True, potentialEnergy=True, temperature=True))

# NVT equilibration for 2 us 
simulation.step(1000000000)