from os.path import exists
from os import remove, rename
from glob import glob
from pydantic import BaseModel, Field, DirectoryPath, FilePath
from typing import Optional, Type, Union
from pathlib import PosixPath, Path


class Settings(BaseModel):
    class General(BaseModel):
        title: str
        pool: DirectoryPath = PosixPath('./pool/')
        basis: str
        pseudopot: Optional[str]
        mode: str = 'vmc_one_mpi1'
        nloc: Optional[int]
        nquad: Optional[int]

    class Ase(BaseModel):
        iase: int = 1
        iforce_analy: int = 0
        node_cutoff: float = 1
        enode_cutoff: float = 0.05

    class Electrons(BaseModel):
        nup: int
        nelec: int

    class Optwf(BaseModel):
        ioptwf: int = 1
        ioptci: int = 1
        ioptjas: int = 1
        ioptorb: int = 1
        nextorb: int = 100
        no_active: int = 1
        nopt_iter: int = 1
        isample_cmat: int = 0
        energy_tol: float = 0.0

    class BlockingVmc(BaseModel):
        vmc_nstep: int = 20
        vmc_nblk: int = 400
        vmc_nblkeq: int = 1
        vmc_nconf_new: int = 0

    general: General
    molecule: Path = Field(prefix="load ")
    basis_num_info: FilePath = Field(prefix="load ")
    determinants: FilePath = Field(prefix="load ")
    orbitals: FilePath = Field(prefix="load ")
    jastrow: FilePath = Field(prefix="load ")
    jastrow_der: FilePath = Field(prefix="load ")
    symmetry: Optional[FilePath] = Field(prefix="load ")
    ase: Optional[Ase]
    electrons: Electrons
    optwf: Optional[Optwf]
    blocking_vmc: Optional[BlockingVmc]

    def write(self, filename='vmc.inp'):
        """
        Write a this dataclass containing the CHAMP configuration to an input file.

        Arguments:
        filename -- input file to write to
        """
        f = open(filename, 'w')
        
        schema = self.schema()['properties']
        for ind, item in enumerate(self):
            if isinstance(item[1], BaseModel):
                f.write("\n%module " + item[0] + "\n")
                for ind2, item2 in enumerate(item[1]):
                    f.write("\t" + item2[0] + " " + str(item2[1]) + "\n")
                f.write("%endmodule\n\n")
            else:
                if 'prefix' in schema[item[0]]:
                    f.write(schema[item[0]]['prefix'] + item[0] + " " + str(item[1]) + '\n')
                else:
                    f.write(item[0] + " " + str(item[1]) + '\n')

        f.close()

    @classmethod
    def read(cls: Type['Model'], filename: Union[str, Path]) -> 'Model':
        """
        Read the CHAMP input file and convert it to a dictionary format

        Arguments:
        filename -- file of the input file to read and write to

        Output:
        Dictionary containing the CHAMP configuration
        """
        # Check if the file exists
        if not exists(filename):
            raise FileNotFoundError(filename + " was not found!")

        output = dict()

        f = open(filename, 'r')

        # Set some variables to keep track of modules
        curmod = ""
        inmod = False

        for line in f:
            # We skip empty lines
            if not line.strip():
                continue
            # Beginning of a module
            elif line.startswith("%module"):
                curmod = line.removeprefix("%module ").removesuffix('\n')
                inmod = True
                output[curmod] = dict()
            # End of a module
            elif line.startswith("%endmodule"):
                inmod = False
            # Loading in a file
            elif line.startswith("load"):
                temp = line.removeprefix("load ").removesuffix('\n')
                temp = temp.split()
                key = temp[0]
                temp.pop(0)
                value = ' '.join(temp)
                output[key] = value
            # All other tags
            else:
                temp = line.removesuffix('\n').split()
                key = temp[0]
                temp.pop(0)
                value = ' '.join(temp)
                if inmod:
                    # If we are in a module, add the tag to that subdictionary.
                    output[curmod][key] = value
                else:
                    output[key] = value

        f.close()

        return cls(**output)

   
    def use_opt_wf(self, filename="vmc.inp"):
        """
        Function to replace the orbitals, determinants and jastrow with the optimized files.

        Arguments:
        filename -- file of the input file to  write to
        """
        opt = ["det_optimal.1.iter*", "orbitals_optimal.1.iter*", "jastrow_optimal.1.iter*"]
        keys = ["determinants", "orbitals", "jastrow"]

        for ind, name in enumerate(opt):
            num = len(glob(name))
            # If there are no optimized WF files, we do not use them
            if num > 0:
                opt[ind] = opt[ind].strip('*') + str(num)
                setattr(self, keys[ind], opt[ind])

        self.write(filename)

def cleanup(*args):
    """
    Remove files created by CHAMP. Add files as arguments to include (if not in default list)
    or exclude (if in default list) files.
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
        if (file.find('*') != -1):
            for f in glob(file):
                remove(f)
        else:
            if exists(file):
                remove(file)

