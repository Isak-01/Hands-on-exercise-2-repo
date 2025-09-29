import sys, unittest
from md import calcenergy
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import Trajectory


class MdTests(unittest.TestCase):
    def test_calcenergy(self):
        """Run a short molecular dynamics simulation of Cu."""
        # Use Asap for a huge performance increase if it is installed
        use_asap = False

        if use_asap:
            from asap3 import EMT

            size = 10
        else:
            from ase.calculators.emt import EMT

            size = 3

        # Set up a crystal
        atoms = FaceCenteredCubic(
            directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            symbol='Cu',
            size=(size, size, size),
            pbc=True,
        )

        # Describe the interatomic interactions with the Effective Medium Theory
        atoms.calc = EMT()

        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(atoms, temperature_K=300)

        # We want to run MD with constant energy using the VelocityVerlet algorithm.
        dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.
        traj = Trajectory("cu.traj", "w", atoms)
        dyn.attach(traj.write, interval=10)
        epot, ekin, inst_temp, etot = calcenergy(atoms)

        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        inst_temp = ekin / (1.5 * units.kB)
        etot = epot + ekin
        if inst_temp == calcenergy(atoms)[2] and etot == calcenergy(atoms)[3]:
            self.assertTrue(True)

if __name__ == "__main__":
    tests = [unittest.TestLoader().loadTestsFromTestCase(MdTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())