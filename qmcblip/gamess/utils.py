"""Tools for simulations using GAMESS."""

import glob
import math
import os
import subprocess
from pathlib import Path
from shutil import copyfile

import ase.io.gamess_us
from ase import Atoms
from pydantic import BaseModel
from qmcblip.champio import Settings


def _write_block(name, args):
    out = [f' ${name.upper()}']
    if isinstance(args, dict):
        for key, val in args.items():
            out.append(f'  {key.upper()}={ase.io.gamess_us._format_value(val)}')
    else:
        out.append(args)
    out.append(' $END')
    return '\n'.join(out)

class WavefunctionCreator:
    """Tool for creating CHAMP-compatible wavefunctions using GAMESS.

    Args:
        atoms (:obj:`ase.Atoms`): atoms object.
        champ_loc (:obj:`Path`): CHAMP basedir location.
    """

    def __init__(self, atoms, champ_loc):
        self.atoms = atoms
        self.champ_path = Path(champ_loc).resolve()
        if not self.champ_path.is_dir():
            raise ValueError("Not a valid CHAMP path!")

        ase.io.gamess_us._write_block = _write_block

        from ase.calculators.gamess_us import GAMESSUS
        self.gamess = GAMESSUS

        self.basename = atoms.symbols.get_chemical_formula()
        self.userscr = None

        if not Path('pool').is_dir():
            os.mkdir('pool')
        if not Path('gamess').is_dir():
            os.mkdir('gamess')

        self.symbols = set(atoms.get_chemical_symbols())
        for atom in self.symbols:
            copyfile(self.champ_path.joinpath('pool/BFD/ECP_gamess/' + atom), 'gamess/' + atom)
            copyfile(self.champ_path.joinpath('pool/BFD/ECP_champ/BFD.gauss_ecp.dat.' + atom),
                    'pool/BFD.gauss_ecp.dat.' + atom)

        self._get_ecp()

        self._create_jastrow()

        self._create_jastrow_der()

        self.vec = None

    def _create_jastrow(self):
        """Create a jastrow file.
        """
        with open('pool/jastrow', 'w', encoding='utf-8') as file:
            file.write("jastrow_parameter   1\n")
            file.write("  5  5  0           norda,nordb,nordc\n")
            file.write("   0.60000000   0.00000000     scalek,a21\n")

            for _, _ in enumerate(self.symbols):
                file.write("   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000" +\
                "   0.00000000 (a(iparmj),iparmj=1,nparma)\n")

            file.write("   0.50000000   0.50000000   0.00000000   0.00000000" +\
                "   0.00000000   0.00000000 (b(iparmj),iparmj=1,nparmb)\n")

            for _, _ in enumerate(self.symbols):
                file.write("   0.00000000   0.00000000   0.00000000   0.00000000" +\
                    "   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000" +\
                    "   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000" +\
                    "   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000" +\
                    "   0.00000000   0.00000000   0.00000000   0.00000000" +\
                    "   (c(iparmj),iparmj=1,nparmc)\n")

            file.write("end\n")

    def _create_jastrow_der(self):
        """Create a jastrow_der file.
        """

        with open('pool/jastrow_der', 'w', encoding='utf-8') as file:
            file.write("jasderiv\n")

            for _, _ in enumerate(self.symbols):
                file.write("4 ")

            file.write("5 0 0 0 0 0 0  nparma,nparmb,nparmc,nparmf\n")

            for _, _ in enumerate(self.symbols):
                file.write("  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)\n")

            file.write("2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)\n")

            for _, _ in enumerate(self.symbols):
                file.write("3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 " +\
                    "(c(iparmj),iparmj=1,nparmc)\n")

            file.write("end\n")

    def _get_ecp(self):
        """Get the ECP potential in GAMESS format.
        """
        self.ecp = {}
        self.ecp_p = {}

        for atom in self.symbols:
            with open('gamess/' + atom, 'r', encoding="utf-8") as file:
                self.ecp_p[atom] = '\n'.join(file.read().split('\n')[:-1])

        has_occured = {}
        for atom in self.symbols:
            has_occured[atom] = False
        for ind, atom in enumerate(self.atoms.get_chemical_symbols()):
            if has_occured[atom]:
                self.ecp[ind] = self.ecp_p[atom].split(" ")[0]
            else:
                self.ecp[ind] = self.ecp_p[atom]
                has_occured[atom] = True

    def _get_vec(self, simname):
        os.chdir('gamess')
        name = self.basename + "_" + simname + ".dat"

        if self.userscr is None:
            raise RuntimeError("USERSCR is not defined")
        if not Path(self.userscr).is_dir():
            raise RuntimeError("USERSCR is not found at: " + self.userscr)

        copyfile(Path(self.userscr).joinpath(name), name)

        self.vec = '\n'.join(str(subprocess.check_output(
            str(self.champ_path.joinpath('tools/interface/getvec.pl')) +\
                 ' -t VEC ' + name,
                 shell=True)).split('\\n')[1:-2]).split(" $END   \n $VEC   \n")[-1]

        os.chdir('..')

    def run_gamess(self, label, **kwargs):
        """Run a GAMESS simulation.

        Args:
            label (:obj:`str`): simulation name.
            **kwargs: other input for GAMESS.
        """
        if 'userscr' in kwargs and kwargs['userscr'] is not None:
            self.userscr = kwargs['userscr']
        if 'vec' in kwargs:
            self._get_vec(kwargs['vec'])
            kwargs['vec'] = self.vec
        if 'ecp' in kwargs and kwargs['ecp']:
            kwargs['ecp'] = self.ecp

        kwargs['label'] = self.basename + '_' + label

        os.chdir('gamess')

        calc = self.gamess(**kwargs)
        calc.calculate(self.atoms)

        if self.userscr is None:
            self.userscr = calc.userscr

        os.chdir('..')

    def convert_to_champ(self, simname):
        """Convert GAMESS output to CHAMP input.
        """
        os.chdir('gamess')

        name = self.basename + '_' + simname

        if not Path(name + ".log").is_file():
            raise ValueError("Log does not exist. Did you already run GAMESS?")

        subprocess.run([str(Path(self.champ_path).joinpath("tools/interface/gamess2qmc")),
                       '-d', '0.0', '-r', '-t', 'initial', name + '.log'], check=True)

        os.chdir('..')

        linenumber = 0
        with open('gamess/' + name + '.lcao', 'r', encoding='utf-8') as file:
            lines = file.readlines()
        with open('pool/orbitals', 'w', encoding='utf-8') as file:
            for line in lines:
                if line[0] not in ('&', '#'):
                    if linenumber == 0:
                        file.write(line.strip('\n') + " 1\n")
                    else:
                        file.write(line)
                    linenumber += 1

        linenumber = 0
        with open('gamess/' + name + '.det.1', 'r', encoding='utf-8') as file:
            lines = file.readlines()
        with open('pool/determinant', 'w', encoding='utf-8') as file:
            for line in lines:
                if line[0] not in ('&', '#'):
                    if linenumber == 0:
                        file.write(line.strip('\n') + " 1\n")
                    else:
                        file.write(line)
                    linenumber += 1

        copyfile('gamess/' + name + '.bfinfo', 'pool/BFD-D.bfinfo')

        for atom in self.symbols:
            for file in glob.glob('gamess/' + name + '.basis.' + atom + ".*"):
                copyfile(file, "pool/BFD-D.basis." + atom)

    def create_champ_input(self):
        """Create CHAMP input.

        Returns:
            :obj:`Settings <qmcblip.champio.Settings>`: CHAMP settings object.
        """
        nelec = 0
        for ind, atom in enumerate(self.atoms.get_chemical_symbols()):
            nelec += self.atoms.symbols.numbers[ind] -\
                int(self.ecp_p[atom].split('\n')[0].split(" ")[2])

        nup = math.ceil(nelec/2)

        settings = Settings(general=Settings.General(title=self.basename, basis="BFD-D",
        pseudopot="BFD"), molecule="molecule.xyz", basis_num_info="pool/BFD-D.bfinfo",
        determinants="pool/determinant", orbitals="pool/orbitals", jastrow="pool/jastrow",
        jastrow_der="pool/jastrow_der", electrons=Settings.Electrons(nup=nup, nelec=nelec))

        return settings

class Presets:
    """Presets for different GAMESS simulations.
    """
    class Template(BaseModel, arbitrary_types_allowed=True):
        """Template for preset.
        """
        userscr: Path = None

        def build(self, champ_loc):
            """Build the wavefunction.

            Args:
                champ_loc (:obj:`Path`): location of CHAMP base directory.

            Returns:
                :obj:`qmcblip.champio.Settings`: CHAMP settings object.
            """
            wf = WavefunctionCreator(self.atoms, champ_loc)
            label = None
            for key, items in self.dict().items():
                if isinstance(items, dict):
                    wf.run_gamess(label=key, userscr=self.userscr, **items)
                    label = key
            wf.convert_to_champ(label)
            return wf.create_champ_input()
    class C2(Template):
        """Preset for C2 wavefunction.
        """
        atoms = Atoms('C2', [(0, 0, -0.7), (0, 0, 0.7)])
        class RHF(BaseModel):
            """RHF simulation for C2.
            """
            contrl = dict(scftyp='RHF', runtyp="ENERGY", mult=1, normf=1, coord="UNIQUE",
                units="BOHR", ecp="READ")
            system = dict(timlim=60, memory=1000000)
            basis = dict(extfil=".TRUE.", gbasis='BFD-D')
            guess = dict(guess="HUCKEL")
            ecp = True

        class CAS(BaseModel):
            """CAS simulation for C2.
            """
            contrl = dict(scftyp='MCSCF', runtyp="ENERGY", mult=1, normf=0, coord="UNIQUE",
                units="BOHR", ecp="READ", ICHARG=0)
            system = dict(mwords=500)
            basis = dict(extfil=".TRUE.", gbasis='BFD-D')
            guess = dict(guess="MOREAD", norb=28, norder=1)
            drt = dict(mxnint=50000, nmcc=2, ndoc=2, nval=2, next=-1, group='NONE', istsym=1,
                fors='.TRUE.', nprt=0)
            mcscf=dict(cistep="GUGA", fullnr=".TRUE.", maxit=200)
            gugdia=dict(prttol=0.0, kprint=2, nstate=1, itermx=200)
            ecp = True
            vec = 'rhf'

        class CI(BaseModel):
            """CI simulation for C2.
            """
            contrl = dict(scftyp='NONE', runtyp="ENERGY", mult=1, normf=0, coord="UNIQUE",
                units="BOHR", ecp="READ", cityp='GUGA', ICHARG=0)
            system = dict(mwords=500)
            basis = dict(extfil=".TRUE.", gbasis='BFD-D')
            guess = dict(guess="MOREAD", norb=28, norder=0, prtmo='.TRUE.')
            cidrt = dict(mxnint=50000, nfzc=2, ndoc=2, nval=2, next=-1, group='NONE', istsym=1,
                fors='.TRUE.', nprt=2)
            mcscf = dict(cistep="GUGA", fullnr=".TRUE.", maxit=200)
            gugdia = dict(prttol=0.0, kprint=2, nstate=1, itermx=200)
            ecp = True
            vec = 'cas'


        rhf: RHF = RHF()
        cas: CAS = CAS()
        ci: CI = CI()
