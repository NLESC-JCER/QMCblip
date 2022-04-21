from ase import Atoms
from qmcblip.gamess.utils import WavefunctionCreator
from qmcblip.champ import CHAMP
from qmcblip.champio import Settings

atoms = Atoms('C2', [(0,0,-0.7), (0,0, 0.7)])
wf = WavefunctionCreator(atoms, '../../champ')

wf.setup_rhf()
wf.setup_cas(system=dict(mwords=500), drt=dict(nmcc=2, ndoc=2, nval=2))
wf.setup_ci(system=dict(mwords=500), cidrt=dict(nfzc=2, ndoc=2, nval=2))
wf.convert_to_champ()
input = wf.create_champ_input()

input.optwf.nopt_iter = 100

atoms.calc = CHAMP(champ_loc='../../champ/bin/vmc.mov1', settings = input)

print(atoms.get_total_energy())
