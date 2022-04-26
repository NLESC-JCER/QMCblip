from qmcblip.gamess.utils import Presets
from qmcblip.champ import CHAMP

atoms = Presets.C2.atoms

input = Presets.C2().build('../../champ')

input.optwf.nopt_iter = 100

atoms.calc = CHAMP(champ_loc='../../champ/bin/vmc.mov1', settings = input)

print(atoms.get_total_energy())
