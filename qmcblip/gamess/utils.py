"""Tools for simulations using GAMESS."""

import glob
import math
import os
import subprocess
from pathlib import Path
from shutil import copyfile

import ase.io.gamess_us
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

        if not Path('pool/jastrow').is_file():
            with open(self.champ_path.joinpath('pool/jastrow.0'), 'r', encoding='utf-8') as file:
                lines = file.readlines()
            with open('pool/jastrow', 'w', encoding='utf-8') as file:
                for line in lines:
                    if line[0] not in ('&', '#'):
                        file.write(line)

        if not Path('pool/jastrow_der').is_file():
            with open(self.champ_path.joinpath('pool/opt.inp'), 'r', encoding='utf-8') as file:
                lines = file.readlines()
            with open('pool/jastrow_der', 'w', encoding='utf-8') as file:
                for line in lines:
                    if line[0] not in ('&', '#'):
                        file.write(line)
        self.vec = None

    def _get_ecp(self):
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



    def setup_rhf(self, **kwargs):
        """Perform RHF setup.

        Args:
            kwargs: arguments for GAMESS.
        """
        name = self.basename + "_rhf"
        os.chdir('gamess')
        # Figure out units
        calc = self.gamess(label=name,
                contrl=dict(scftyp='RHF', runtyp="ENERGY", mult=1, normf=1, coord="UNIQUE",
                    units="BOHR", ecp="READ"),
                system=dict(timlim=60, memory=1000000),
                basis=dict(extfil=".TRUE.", gbasis='BFD-D'),
                guess=dict(guess="HUCKEL"),
                ecp=self.ecp)

        for key, item in kwargs.items():
            if isinstance(item, dict):
                for key2, item2 in item.items():
                    calc.parameters[key][key2] = item2
            elif key == "userscr":
                calc.userscr = item
            else:
                calc.parameters[key] = item

        calc.calculate(self.atoms)

        if calc.userscr is None:
            raise RuntimeError("USERSCR is not defined")
        if not Path(calc.userscr).is_dir():
            raise RuntimeError("USERSCR is not found at: " + calc.userscr)

        copyfile(Path(calc.userscr).joinpath(name + '.dat'), name + '.dat')

        self.vec = '\n'.join(str(subprocess.check_output(
            str(self.champ_path.joinpath('tools/interface/getvec.pl')) +\
                 ' -t VEC ' + name+ '.dat', shell=True)).split('\\n')[1:-2])

        os.chdir('..')

    def setup_cas(self, **kwargs):
        """Perform CAS setup.

        Args:
            kwargs: arguments for GAMESS.
        """
        if self.vec is None:
            raise ValueError("No vectors.")

        name = self.basename + "_cas"
        os.chdir('gamess')
        # Figure out units
        calc = self.gamess(label=name,
                contrl=dict(scftyp='MCSCF', runtyp="ENERGY", mult=1, normf=0, coord="UNIQUE",
                    units="BOHR", ecp="READ", ICHARG=0),
                system=dict(mwords=1500),
                basis=dict(extfil=".TRUE.", gbasis='BFD-D'),
                guess=dict(guess="MOREAD", norb=28, norder=1),
                drt=dict(mxnint=50000, nmcc=2, ndoc=2, nval=0, next=-1, group='NONE',
                    istsym=1, fors='.TRUE.', nprt=0),
                mcscf=dict(cistep="GUGA", fullnr=".TRUE.", maxit=200),
                gugdia=dict(prttol=0.0, kprint=2, nstate=1, itermx=200),
                ecp=self.ecp,
                vec=self.vec)

        for key, item in kwargs.items():
            if isinstance(item, dict):
                for key2, item2 in item.items():
                    calc.parameters[key][key2] = item2
            elif key == "userscr":
                calc.userscr = item
            else:
                calc.parameters[key] = item

        calc.calculate(self.atoms)

        if calc.userscr is None:
            raise RuntimeError("USERSCR is not defined")
        if not Path(calc.userscr).is_dir():
            raise RuntimeError("USERSCR is not found at: " + calc.userscr)

        copyfile(Path(calc.userscr).joinpath(name + '.dat'), name + '.dat')

        self.vec = '\n'.join(str(subprocess.check_output(
            str(self.champ_path.joinpath('tools/interface/getvec.pl')) +\
                 ' -t VEC ' + name+ '.dat',
                 shell=True)).split('\\n')[1:-2]).split(" $END   \n $VEC   \n")[1]

        os.chdir('..')

    def setup_ci(self, **kwargs):
        """Perform CI setup.

        Args:
            kwargs: arguments for GAMESS.
        """
        if self.vec is None:
            raise ValueError("No vectors.")

        name = self.basename + "_ci"
        os.chdir('gamess')
        # Figure out units
        calc = self.gamess(label=name,
                contrl=dict(scftyp='NONE', runtyp="ENERGY", mult=1, normf=0, coord="UNIQUE",
                    units="BOHR", ecp="READ", cityp='GUGA', ICHARG=0),
                system=dict(mwords=1500),
                basis=dict(extfil=".TRUE.", gbasis='BFD-D'),
                guess=dict(guess="MOREAD", norb=28, norder=0, prtmo='.TRUE.'),
                cidrt=dict(mxnint=50000, nfzc=2, ndoc=2, nval=0, next=-1, group='NONE', istsym=1,
                    fors='.TRUE.', nprt=2),
                mcscf=dict(cistep="GUGA", fullnr=".TRUE.", maxit=200),
                gugdia=dict(prttol=0.0, kprint=2, nstate=1, itermx=200),
                ecp=self.ecp,
                vec=self.vec)

        for key, item in kwargs.items():
            if isinstance(item, dict):
                for key2, item2 in item.items():
                    calc.parameters[key][key2] = item2
            elif key == "userscr":
                calc.userscr = item
            else:
                calc.parameters[key] = item

        calc.calculate(self.atoms)

        if calc.userscr is None:
            raise RuntimeError("USERSCR is not defined")
        if not Path(calc.userscr).is_dir():
            raise RuntimeError("USERSCR is not found at: " + calc.userscr)

        copyfile(Path(calc.userscr).joinpath(name + '.dat'), name + '.dat')

        os.chdir('..')

    def convert_to_champ(self):
        """Convert GAMESS output to CHAMP input.
        """
        os.chdir('gamess')

        name = self.basename + '_ci'

        if not Path(name + ".log").is_file():
            raise ValueError("CI log does not exist. Did you already run GAMESS?")

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
