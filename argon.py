"""argon.py"""

# Constant-temperature constant-pressure MD simulation of Argon

from MMTK import *
from MMTK.ForceFields import LennardJonesForceField
from MMTK.Environment import NoseThermostat, AndersenBarostat
from MMTK.Trajectory import Trajectory, TrajectoryOutput, LogOutput
from MMTK.Dynamics import VelocityVerletIntegrator, Velocity Scaler, \
                            TranslationRemover, BarostatReset
import string
from Scientific.IO.TextFile import TextFile

# Open the config file and read box size
conf_file = TextFile('argon.conf.gz')
lx, ly, lz = map(string.atof, string.split(conf_file.readline()))

# Construct periodic universe using Lennard-Jones (noble gas) force field
# with cutoff of 15 Angstroms
universe = OrthorhombicPeriodicUniverse((lx*Units.Ang, ly*Units.Ang, lz*Units.Ang),
                                        LennardJonesForceField(15.*Units.Ang))

# Read the atom positions and construct the atoms
while 1:
    line = conf_file.readline()
    if not line: 
        break
    x, y, z, = map(string.atof, string.split(line))
    universe.addObject(Atom('Ar', position=Vector(x*Units.Ang, y*Units.Ang, z*Units.Ang)))

# Define thermodynamic parameters
temperature = 94.4*Units.K
pressure = 1.*Units.atm

# Add thermostat and barostat
universe.thermostat = NoseThermostat(temperature)
universe.barostat = AndersenBarostat(pressure)

# Initialise velocities
universe.initializeVelocitiesToTemperature(temperature)

# Create trajectory and integrator
trajectory = Trajectory(universe, "argon_npt.nc", "w", "Argon NPT test")
integrator = VelocityVerletIntegrator(universe, delta_t=10*Units.fs)

# Periodical actions for trajectory output and text log output
output_actions = [TrajectoryOutput(trajectory, ('configuration', 'energy', 'thermodynamic',
                                                'time', 'auxiliary'), 0, None, 20),
                    LogOutput("argon.log", ('time', 'energy'), 0, None, 100)]

# Do some equilibration steps, rescaling velocities and resetting barostat in regular intervals
integrator(steps = 2000, actions = [TranslationRemover(0, None, 100), 
                                    VelocityScaler(temperature, 0.1*temperature, 0, None, 100),
                                    BarostatReset(100)] + output_actions)

# Do some "production" steps
integrator(steps = 2000, actions = [TranslationRemover(0, None, 100)] + output_actions)

# Close trajectory
trajectory.close()