"""Functions and classes for CHAMP I/O."""
from glob import glob
from os import remove
from os.path import exists
from pathlib import Path, PosixPath
from typing import Optional, Type, Union

from pydantic import BaseModel, DirectoryPath, Field, FilePath, Extra


class Settings(BaseModel, extra=Extra.allow):
    """Data class containing CHAMP configuration.

    This class can hold the neccessery configuration to run CHAMP.
    """
    class General(BaseModel, extra=Extra.allow):
        """General module class.
        """

        title: str
        """:obj:`str`: title of the simulation."""

        pool: DirectoryPath = Field(PosixPath('./pool/'), postfix="/")
        """:obj:`path <pathlib.Path>`, optional: location of the pool directory."""

        basis: str
        """:obj:`str`: basis to use (located in pool)."""

        pseudopot: Optional[str]
        """:obj:`str`, optional: pseudopotential to use."""

        mode: str = 'vmc_one_mpi1'
        """:obj:`str`, optional: QMC mode."""

        seed: int = 1837465927472523
        """int, optional: seed for CHAMP."""

        eunit: str = "Hartrees"
        """:obj:`str`, optional: energy units."""

    class Ase(BaseModel, extra=Extra.allow):
        """ASE module class.
        """

        iase: int = 1
        iforce_analy: int = 0
        node_cutoff: int = 1
        enode_cutoff: float = 0.05

    class Electrons(BaseModel, extra=Extra.allow):
        """Electrons module class.
        """

        nup: int
        """:obj:`int`: amount of electrons with upspin."""

        nelec: int
        """:obj:`int`: total amount of electrons."""

    class Optwf(BaseModel, extra=Extra.allow):
        """Wavefunction optimization module class.
        """

        ioptwf: int = 1
        ioptci: int = 1
        ioptjas: int = 1
        ioptorb: int = 1
        nextorb: int = 100
        no_active: int = 1
        nopt_iter: int = 1
        isample_cmat: int = 0
        energy_tol: float = 0.0

    class BlockingVmc(BaseModel, extra=Extra.allow):
        """VMC module class.
        """

        vmc_nstep: int = 20
        vmc_nblk: int = 400
        vmc_nblkeq: int = 1
        vmc_nconf_new: int = 0

    class Pseudo(BaseModel, extra=Extra.allow):
        """Pseudopotential module class.
        """

        nloc: int = 4
        nquad: int = 6

    general: General
    molecule: Path = Field(prefix="load ")
    """:obj:`path <pathlib.Path>`: path to molecule geometry."""

    basis_num_info: FilePath = Field(prefix="load ")
    """:obj:`path <pathlib.Path>`: path to basis info."""

    determinants: FilePath = Field(prefix="load ")
    """:obj:`path <pathlib.Path>`: path to the determinants file."""

    orbitals: FilePath = Field(prefix="load ")
    """:obj:`path <pathlib.Path>`: path to the orbitals file."""

    jastrow: FilePath = Field(prefix="load ")
    """:obj:`path <pathlib.Path>`: path to the jastrow file."""

    jastrow_der: FilePath = Field(prefix="load ")
    """:obj:`path <pathlib.Path>`: path to the jastrow derivatives file."""

    symmetry: Optional[FilePath] = Field(prefix="load ")
    """optional: path to the symmetry file."""

    ase: Ase = Ase()
    electrons: Electrons
    optwf: Optwf = Optwf()
    pseudo: Pseudo = Pseudo()
    blocking_vmc: BlockingVmc = BlockingVmc()

    def write(self, filename='vmc.inp'):
        """Write a this dataclass containing the CHAMP configuration to an input file.

        Args:
            filename (:obj:`str`, optional): input file to write to.
        """

        input_file = open(filename, 'w', encoding="utf-8")

        schema = self.schema()['properties']
        for _, item in enumerate(self):
            if isinstance(item[1], BaseModel):
                if callable(getattr(item[1], 'schema')):
                    schema2 = item[1].schema()['properties']
                else:
                    schema2 = []
                input_file.write("\n%module " + item[0] + "\n")
                for _, item2 in enumerate(item[1]):
                    if item2[0] in schema2 and 'postfix' in schema2[item2[0]]:
                        if item2[1] is not None:
                            input_file.write("\t" + item2[0] + " " + str(item2[1]) +\
                                 schema2[item2[0]]['postfix'] + "\n")
                    else:
                        if item2[1] is not None:
                            input_file.write("\t" + item2[0] + " " + str(item2[1]) + "\n")
                input_file.write("%endmodule\n\n")
            else:
                if item[0] in schema and 'prefix' in schema[item[0]]:
                    if item[1] is not None:
                        input_file.write(schema[item[0]]['prefix'] + item[0] + " " +\
                             str(item[1]) + '\n')
                else:
                    if item[1] is not None:
                        input_file.write(item[0] + " " + str(item[1]) + '\n')

        input_file.close()

    @classmethod
    def read(cls: Type['BaseModel'], filename: Union[str, Path]) -> 'BaseModel':
        """Read the CHAMP input file and convert it to a dictionary format

        Args:
            filename (:obj:`str`, optional): file of the input file to read and write to.

        Returns:
            :obj:`Settings`: Settings object containing CHAMP configuration.
        """
        # Check if the file exists
        if not exists(filename):
            raise FileNotFoundError(filename + " was not found!")

        path = PosixPath(filename).resolve().parent

        try:
            path = path.relative_to(Path.cwd()).as_posix() + "/"
        except:
            path = path.as_posix() + "/"

        output = {}

        input_file = open(filename, 'r', encoding='utf-8')

        # Set some variables to keep track of modules
        curmod = ""
        inmod = False

        for line in input_file:
            # We skip empty lines
            if not line.strip():
                continue
            # Beginning of a module
            if line.startswith("%module"):
                curmod = line[8:-1]
                inmod = True
                output[curmod] = {}
            # End of a module
            elif line.startswith("%endmodule"):
                inmod = False
            # Loading in a file
            elif line.startswith("load"):
                temp = line[5:-1]
                temp = temp.split()
                key = temp[0]
                temp.pop(0)
                value = ' '.join(temp)
                output[key] = path + value
            # All other tags
            else:
                temp = line[:-1].split()
                key = temp[0]
                temp.pop(0)
                value = ' '.join(temp)
                if inmod:
                    # If we are in a module, add the tag to that subdictionary.
                    if key != "pool":
                        output[curmod][key] = value
                    else:
                        output[curmod][key] = path + value
                else:
                    output[key] = value

        input_file.close()

        return cls(**output)


    def use_opt_wf(self):
        """Function to replace the orbitals, determinants and jastrow with the optimized files.
        """
        opt = ["det_optimal.1.iter*", "orbitals_optimal.1.iter*", "jastrow_optimal.1.iter*"]
        keys = ["determinants", "orbitals", "jastrow"]

        for ind, name in enumerate(opt):
            num = len(glob(name))
            # If there are no optimized WF files, we do not use them
            index = min(num, self.optwf.nopt_iter)
            if num >= index and index > 0:
                opt[ind] = opt[ind].strip('*') + str(index)
                setattr(self, keys[ind], opt[ind])

    def todict(self):
        """Return a dictionary of the class.

        Returns:
            :obj:`str`: json dictionary.
        """
        return self.json(exclude_none=True)

def cleanup(*args):
    """Remove files created by CHAMP.

    This function can cleanup the directory of the simulation.
    Add files as arguments to include (if not in default list)
    or exclude (if in default list) files.

    Args:
        *args: filesnames to include (if not in default list) or exclude (if in default list).
    """
    # Default list of the files that will be removed
    standard = ["force_analytic", "restart_vmc", "output.log", "parser.log", "orbitals_optimal.1.*",
                "jastrow_optimal.1.*", "det_optimal.1.*", "geo_optimal.*", "geo_optimal_final",
                "mc_configs_start", "nohup.out", "vmc_temp.inp"]

    # Add or remove file depending on the input arguments
    for item in args:
        if item in standard:
            standard.remove(item)
        else:
            standard.append(item)

    for file in standard:
        # Take care of the wildcard files
        if file.find('*') != -1:
            for file2 in glob(file):
                remove(file2)
        else:
            if exists(file):
                remove(file)
