"""Quickly run FLARE simulations."""
from typing import Any, List

import flare_pp._C_flare as flare_pp
from ase import units
from flare_pp.sparse_gp import SGP_Wrapper
from flare_pp.sparse_gp_calculator import SGP_Calculator
from pydantic import BaseModel, Extra

from flare.ase.calculator import FLARE_Calculator
from flare.gp import GaussianProcess
from flare.utils.parameter_helper import ParameterHelper

from .otf import C_ASE_OTF as ASE_OTF


class OTFSettings(BaseModel, extra=Extra.allow):
    """Dataclass containing FLARE(++) configuration.
    """

    class Theory(BaseModel, extra=Extra.allow):
        """Abstract baseclass.
        """

    class FLARE(Theory):
        """FLARE dataclass.
        """

        kernels: List[str] = ['twobody', 'threebody']
        """:obj:`List[str]`, optional: n-body functions to use."""

        cutoffs: List[float] = [5.0, 3.5]
        """:obj:`List[float]`, optional: cutoff for the n-body functions."""

        random: bool = True
        """:obj:`bool`, optional: randomize hyperparameters."""

        hyp_labels: List[str] = ['sig2','ls2','sig3','ls3','noise']
        """:obj:`List[str]`, optional: hyperparameter labels."""

        opt_algorithm: str = 'L-BFGS-B'
        """:obj:`str`, optional: hyperparamter optimization algorithm."""

        n_cpu: int = 1
        """:obj:`int`, optional: amount of CPU cores to run FLARE on."""

        update_style: str = "add_n"
        """:obj:`str`, optional: update algorithm for GP."""

        update_threshold: float = None
        """:obj:`str`, optional: update threshold for GP, in eV/Ang."""

        gp_model: Any = None
        flare_calc: Any = None

        def get_calc(self, atoms):
            """Make a FLARE calculator object.

            Args:
                atoms (:obj:`Atoms <ase.Atoms>`): atoms object.

            Returns:
                :obj:`FLARE_Calculator <flare.ase.calculator.FLARE_Calculator>`:
                    FLARE calculator object.
            """
            parameters = {}
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

            self.gp_model = GaussianProcess(
                kernels = self.kernels,
                component = 'mc',
                hyps = hyps,
                cutoffs = cut,
                hyp_labels = self.hyp_labels,
                opt_algorithm = self.opt_algorithm,
                parallel = parallel,
                n_cpus = self.n_cpu
            )

            self.flare_calc = FLARE_Calculator(self.gp_model,
                                    par = True,
                                    mgp_model = None,
                                    use_mapping = False)
            return self.flare_calc


    class FLAREPP(Theory):
        """FLARE++ dataclass.
        """

        update_style: str = "threshold"
        """:obj:`str`, optional: update algorithm for GP."""

        update_threshold: float = 0.005
        """:obj:`str`, optional: update threshold for GP, in eV/Ang."""

        opt_algorithm: str = 'L-BFGS-B'
        """:obj:`str`, optional: hyperparamter optimization algorithm."""

        max_iterations: int = 10
        """:obj:`int`, optional: maximum number of hyperparameter optimization
            steps per MD step."""

        variance_type: str = 'local'
        """:obj:`str`, optional: uncertainty type on energy."""

        sigma_e: float = 0.12
        """:obj:`float`, optional: energy noise per atom, in kcal/mol."""

        sigma_f: float = 0.115
        """:obj:`float`, optional: force noise, in kcal/mol/Ang."""

        sigma_s: float = 0.014
        """:obj:`float`, optional: stress noise, in kcal/Ang^3."""

        cutoff: float = 5.0
        """:obj:`float`, optional: cutoff of kernel in Ang."""

        sigma: float = 2.0
        """:obj:`float`, optional: kernel thing in eV."""

        power: int = 2
        """:obj:`int`, optional: power of the kernel."""

        cutoff_function: str = "quadratic"
        """:obj:`str`, optional: cutoff function."""

        radial_basis: str = "chebyshev"
        """:obj:`str`, optional: radial basis."""

        N: int = 12
        """:obj:`int`, optional: number of radial basis functions."""

        lmax: int = 3
        """:obj:`int`, optional: largest L included in spherical harmonics."""

        kernel: Any = None
        descriptor_calculator: Any = None
        gp_model: Any = None
        flare_calc: Any = None

        def get_calc(self, atoms):
            """Make a FLARE++ calculator object.

            Args:
                atoms (:obj:`Atoms <ase.Atoms>`): atoms object.

            Returns:
                :obj:`SGP_Calculator <flare_pp.sparse_gp_calculator.SGP_Calculator>`:
                    FLARE++ calculator object.
            """
            species_map = {}
            for ind, num in enumerate(set(atoms.symbols.numbers)):
                species_map[num] = ind
            self.kernel = flare_pp.NormalizedDotProduct(self.sigma, self.power)
            radial_hyps = [0., self.cutoff]
            cutoff_hyps = []
            n_species = len(species_map)
            descriptor_settings = [n_species, self.N, self.lmax]
            self.descriptor_calculator = flare_pp.B2(
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

            self.gp_model = SGP_Wrapper([self.kernel], [self.descriptor_calculator], self.cutoff,
                                sigma_e, sigma_f, sigma_s, species_map,
                                variance_type=self.variance_type,
                                stress_training=False,
                                opt_method=self.opt_algorithm,
                                bounds=bounds,
                                max_iterations=self.max_iterations)

            self.flare_calc = SGP_Calculator(self.gp_model)

            return self.flare_calc



    theory: Theory = FLARE()
    """:obj:`Theory`: FLARE or FLAREPP."""

    output_name: str = 'OTF'
    """:obj:`str`, optional: name of the files to write to."""

    std_tolerance_factor: float = -0.01
    """:obj:`float`, optional: standard tolerance with respect to noise.
        Negative for absolute (in eV/Ang).
    """

    min_steps_with_model: int = 0
    """:obj:`int`, optional: minimum steps with model in between ab-initio calls."""

    freeze_hyps: int = 10
    """:obj:`int`, optional: freeze hyperparameters after this many ab-initio calls."""

    write_model: int = 0
    """:obj:`int`, optional: keep at 0 for FLARE++."""




def quicksim(atoms, timestep, steps, calc, otfsettings = OTFSettings(), changes = []):
    """A quick-to-setup simulation using FLARE or FLARE++.

    This function allows the user to do a quick simulation using
    ab-initio calculations and FLARE(++), but at the cost of less flexibility.

    Args:
        atoms (:obj:`Atoms <ase.Atoms>`): atoms object.
        timestep (:obj:`float`): timestep of MD in fs.
        steps (:obj:`int`): amount of MD steps.
        otfsettings (:obj:`OTFSettings`, optional): OTF settings object,
            using defaults if not given.
        changes (optional): array containing settings to update for CHAMP (advanced).
    """
    n_atoms = len(atoms)

    flare_calculator = otfsettings.theory.get_calc(atoms)

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


    sim = ASE_OTF(atoms,
                timestep = timestep * units.fs,
                number_of_steps = steps,
                dft_calc = calc,
                md_engine = md_engine,
                md_kwargs = md_kwargs,
                update_settings = changes,
                calculator = flare_calculator,
                **otf_params)

    sim.run()
