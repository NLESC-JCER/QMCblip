# Import numpy and matplotlib
import numpy as np

# flare imports
from flare import otf_parser
from flare.gp import GaussianProcess
from flare.ase.calculator import FLARE_Calculator
from flare.utils.parameter_helper import ParameterHelper
import flare_pp._C_flare as flare_pp
from flare_pp.sparse_gp import SGP_Wrapper
from flare_pp.sparse_gp_calculator import SGP_Calculator

# Dataclass imports
from pydantic import BaseModel
from typing import Union, List

# CAF imports
from .otf import C_ASE_OTF as ASE_OTF


class OTFSettings(BaseModel):

    class FLARE(BaseModel):
        kernels: List[str] = ['twobody', 'threebody']
        cutoffs: List[float] = [5.0, 3.5]
        random: bool = True
        hyp_labels: List[str] = ['sig2','ls2','sig3','ls3','noise']
        opt_algorithm: str = 'L-BFGS-B'
        n_cpu: int = 1
        update_style: str = "add_n"
        update_threshold: float = None

        def get_model(self, atoms):
            """Make a FLARE calculator object.

            Args:
                atoms (:obj:`Atoms <ase.Atoms>`): atoms object.
            
            Returns:
                :obj:`FLARE_Calculator <flare.ase.calculator.FLARE_Calculator>`: 
                    FLARE calculator object.
            """
            parameters = dict()
            for ind, nbody in enumerate(self.kernels):
                parameters['cutoff_'+nbody] = self.cutoffs[ind]
            pm = ParameterHelper(
                kernels = self.kernels,   
                random = self.random,
                parameters=parameters
            )

            hm = pm.as_dict()
            hyps = hm['hyps']
            cut = hm['cutoffs']

            if self.n_cpu > 1:
                parallel = True
            else:
                parallel = False

            gp_model = GaussianProcess(
                kernels = self.kernels,
                component = 'mc',
                hyps = hyps,
                cutoffs = cut,
                hyp_labels = self.hyp_labels,
                opt_algorithm = self.opt_algorithm,
                parallel = parallel,
                n_cpus = self.n_cpu
            )

            calc = FLARE_Calculator(gp_model, 
                                    par = True, 
                                    mgp_model = None,
                                    use_mapping = False)
            return calc


    class FLAREPP(BaseModel):
        update_style: str = "threshold" 
        update_threshold: float = 0.005
        opt_algorithm: str = 'L-BFGS-B'
        max_iterations: int = 10
        variance_type: str = 'local'
        sigma_e: float = 0.12
        sigma_f: float = 0.115
        simga_s: float = 0.014
        cutoff: float = 5
        sigma: float = 2
        power: int = 2
        cutoff_function: str = "quadratic"
        radial_basis: str = "chebyshev"
        N: int = 12
        lmax: int = 3

        def get_calc(self, atoms):
            """Make a FLARE++ calculator object.

            Args:
                atoms (:obj:`Atoms <ase.Atoms>`): atoms object.
            
            Returns:
                :obj:`SGP_Calculator <flare_pp.sparse_gp_calculator.SGP_Calculator>`: 
                    FLARE++ calculator object.
            """
            species_map = dict()
            for ind, num in enumerate(set(atoms.symbols.numbers)):
                species_map[num] = ind
            cutoff = self.cutoff
            sigma = self.sigma
            power = self.power
            kernel = flare_pp.NormalizedDotProduct(sigma, power)
            many_body_cutoffs = [cutoff]
            radial_hyps = [0., cutoff]
            cutoff_hyps = []
            n_species = len(species_map)
            descriptor_settings = [n_species, self.N, self.lmax]
            descriptor_calculator = flare_pp.B2(
                                                self.radial_basis,
                                                self.cutoff_function,
                                                radial_hyps,
                                                cutoff_hyps,
                                                descriptor_settings
                                                )

            sigma_e = self.sigma_e * len(atoms)
            sigma_f = self.sigma_f 
            sigma_s = self.sigma_s

            bounds = [(None, None), (sigma_e, None), (None, None), (None, None)]

            gp_model = SGP_Wrapper([kernel], [descriptor_calculator], cutoff,
                                sigma_e, sigma_f, sigma_s, species_map,
                                variance_type=self.variance_type,
                                stress_training=False,
                                opt_method=self.opt_method,
                                bounds=bounds,
                                max_iterations=self.max_iterations)

            calc = SGP_Calculator(gp_model)

            return calc



    theory: Union[FLARE, FLAREPP] = FLARE
    output_name: str = 'OTF'
    std_tolerance_factor: float = -0.01
    min_steps_with_model: int = 0
    freeze_hyps: int = 10
    write_model: int = 0




def quicksim(atoms, timestep, steps, otfsettings = OTFSettings(), changes = []):
    """A quick-to-setup simulation using FLARE or FLARE++.

    This function allows the user to do a quick simulation using
    ab-initio calculations and FLARE(++), but at the cost of less flexibility.

    Args:
        atoms (:obj:`Atoms <ase.Atoms>`): atoms object.
        timestep (:obj:`float`): timestep of MD in ps.
        steps (:obj:`int`): amount of MD steps.
        otfsettings (:obj:`OTFSettings`, optional): OTF settings object, 
            using defaults if not given.
        changes (optional): array containing settings to update for CHAMP (advanced).
    """
    n_atoms = len(atoms)

    flare_calculator = otfsettings.theory.get_calc(atoms)

    atoms.set_calculator(flare_calculator)

    md_engine = 'CustomVerlet'
    md_kwargs = {}

    otf_params = {'init_atoms': None,
                'output_name': otfsettings.output_name,
                'std_tolerance_factor': otfsettings.std_tolerance_factor,
                'max_atoms_added' : n_atoms,
                'force_only' : False,
                'freeze_hyps': otfsettings.freeze_hyps,
                'write_model': otfsettings.write_model,
                'min_steps_with_model': otfsettings.min_steps_with_model,
                'update_style': otfsettings.theory.update_style,
                'update_threshold': otfsettings.theory.update_threshold
                }


    otf = ASE_OTF(atoms, 
                    timestep = timestep * units.fs,
                    number_of_steps = steps,
                    dft_calc = calc,
                    md_engine = md_engine,
                    md_kwargs = md_kwargs,
                    update_settings = changes,
                    **otf_params)

    otf.run()
