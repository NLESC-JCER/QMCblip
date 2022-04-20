"""On-the-fly training classes adapted for CHAMP."""
# MIT License

# Copyright (c) 2019 Harvard University

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

from time import time
from copy import deepcopy
import numpy as np

from flare.otf import OTF
from flare.ase.otf import ASE_OTF
from flare.ase.atoms import FLARE_Atoms
import flare.ase.dft as dft_source
from flare.utils.learner import is_std_in_bound

from ase import units

from .verlet import CustomVerlet

class C_OTF(OTF):
    """Trains a Gaussian process force field on the fly during molecular dynamics.

    Args:
        dt (float): MD timestep.
        number_of_steps (:obj:`int`): Number of timesteps in the training
            simulation.
        prev_pos_init ([type], optional): Previous positions. Defaults
            to None.
        rescale_steps (List[:obj:`int`], optional): List of frames for which the
            velocities of the atoms are rescaled. Defaults to [].
        rescale_temps (List[:obj:`int`], optional): List of rescaled temperatures.
            Defaults to [].
        gp (gp.GaussianProcess): Initial GP model.
        calculate_energy (:obj:`bool`, optional): If True, the energy of each
            frame is calculated with the GP. Defaults to False.
        calculate_efs (:obj:`bool`, optional): If True, the energy and stress of each
            frame is calculated with the GP. Defaults to False.
        write_model (:obj:`int`, optional): If 0, write never. If 1, write at
            end of run. If 2, write after each training and end of run.
            If 3, write after each time atoms are added and end of run.
            If 4, write after each training and end of run, and back up
            after each write.
        force_only (:obj:`bool`, optional): If True, only use forces for training.
            Default to False, use forces, energy and stress for training.
        std_tolerance_factor (float, optional): Threshold that determines
            when DFT is called. Specifies a multiple of the current noise
            hyperparameter. If the epistemic uncertainty on a force
            component exceeds this value, DFT is called. Defaults to 1.
        skip (:obj:`int`, optional): Number of frames that are skipped when
            dumping to the output file. Defaults to 0.
        init_atoms (List[:obj:`int`], optional): List of atoms from the input
            structure whose local environments and force components are
            used to train the initial GP model. If None is specified, all
            atoms are used to train the initial GP. Defaults to None.
        output_name (:obj:`str`, optional): Name of the output file. Defaults to
            'otf_run'.
        max_atoms_added (:obj:`int`, optional): Number of atoms added each time
            DFT is called. Defaults to 1.
        freeze_hyps (:obj:`int`, optional): Specifies the number of times the
            hyperparameters of the GP are optimized. After this many
            updates to the GP, the hyperparameters are frozen.
            Defaults to 10.
        min_steps_with_model (:obj:`int`, optional): Minimum number of steps the
            model takes in between calls to DFT. Defaults to 0.
        force_source (Union[:obj:`str`, object], optional): DFT code used to calculate
            ab initio forces during training. A custom module can be used here
            in place of the DFT modules available in the FLARE package. The
            module must contain two functions: parse_dft_input, which takes a
            file name (in string format) as input and returns the positions,
            species, cell, and masses of a structure of atoms; and run_dft_par,
            which takes a number of DFT related inputs and returns the forces
            on all atoms.  Defaults to "qe".
        npool (:obj:`int`, optional): Number of k-point pools for DFT
            calculations. Defaults to None.
        mpi (:obj:`str`, optional): Determines how mpi is called. Defaults to
            "srun".
        dft_loc (:obj:`str`): Location of DFT executable.
        dft_input (:obj:`str`): Input file.
        dft_output (:obj:`str`): Output file.
        dft_kwargs ([type], optional): Additional arguments which are
            passed when DFT is called; keyword arguments vary based on the
            program (e.g. ESPRESSO vs. VASP). Defaults to None.
        store_dft_output (Tuple[Union[:obj:`str`,List[:obj:`str`]],:obj:`str`], optional):
            After DFT calculations are called, copy the file or files
            specified in the first element of the tuple to a directory
            specified as the second element of the tuple.
            Useful when DFT calculations are expensive and want to be kept
            for later use. The first element of the tuple can either be a
            single file name, or a list of several. Copied files will be
            prepended with the date and time with the format
            'Year.Month.Day:Hour:Minute:Second:'.
        n_cpus (:obj:`int`, optional): Number of cpus used during training.
            Defaults to 1.
    """

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def run(self):
        """Performs an on-the-fly training run.

        If OTF has store_dft_output set, then the specified DFT files will
        be copied with the current date and time prepended in the format
        'Year.Month.Day:Hour:Minute:Second:'.
        """

        optional_dict = {"Restart": self.curr_step}
        self.output.write_header(
            str(self.gp),
            self.dt,
            self.number_of_steps,
            self.structure,
            self.std_tolerance,
            optional_dict,
        )

        counter = 0
        self.start_time = time()

        while self.curr_step < self.number_of_steps:
            # run DFT and train initial model if first step and DFT is on
            if (
                (self.curr_step == 0)
                and (self.std_tolerance != 0)
                and (len(self.gp.training_data) == 0)
            ):

                # Are the recorded forces from the GP or DFT in ASE OTF?
                # When DFT is called, ASE energy, forces, and stresses should
                # get updated.
                self.initialize_train()

            # after step 1, try predicting with GP model
            else:
                # compute forces and stds with GP
                self.dft_step = False
                self.compute_properties()

                # get max uncertainty atoms
                std_in_bound, target_atoms = is_std_in_bound(
                    self.std_tolerance,
                    self.gp.force_noise,
                    self.structure,
                    max_atoms_added=self.max_atoms_added,
                    update_style=self.update_style,
                    update_threshold=self.update_threshold,
                )

                steps_since_dft = self.curr_step - self.last_dft_step
                if (not std_in_bound) and (steps_since_dft > self.min_steps_with_model):
                    # record GP forces
                    self.update_temperature()
                    self.record_state()
                    gp_energy = self.structure.potential_energy
                    gp_forces = deepcopy(self.structure.forces)
                    gp_stress = deepcopy(self.structure.stress)

                    # run DFT and record forces
                    self.dft_step = True
                    self.last_dft_step = self.curr_step
                    self.run_dft()
                    dft_frcs = deepcopy(self.structure.forces)
                    dft_stress = deepcopy(self.structure.stress)
                    dft_energy = self.structure.potential_energy

                    # run MD step & record the state
                    self.record_state()

                    # record DFT data into an .xyz file with filename self.dft_xyz.
                    # the file includes the structure, e/f/s labels and atomic
                    # indices of environments added to gp
                    self.record_dft_data(self.structure, target_atoms)

                    # compute mae and write to output
                    self.compute_mae(
                        gp_energy,
                        gp_forces,
                        gp_stress,
                        dft_energy,
                        dft_frcs,
                        dft_stress,
                    )

                    # add max uncertainty atoms to training set
                    self.update_gp(
                        target_atoms,
                        dft_frcs,
                        dft_stress=dft_stress,
                        dft_energy=dft_energy,
                    )

                    if self.write_model == 4:
                        self.checkpoint()
                        self.backup_checkpoint()

            # write gp forces
            if counter >= self.skip and not self.dft_step:
                self.update_temperature()
                self.record_state()
                counter = 0

            counter += 1

            self.md_step(forces=self.structure.forces)  # update positions by Verlet
            self.rescale_temperature(self.structure.positions)

            if self.write_model == 3:
                self.checkpoint()

        self.output.conclude_run()

        if self.write_model >= 1:
            self.write_gp()
            self.checkpoint()

    def initialize_train(self):
        # call dft and update positions
        self.run_dft()
        dft_frcs = deepcopy(self.structure.forces)
        dft_stress = deepcopy(self.structure.stress)
        dft_energy = self.structure.potential_energy

        self.update_temperature()
        self.record_state()
        self.record_dft_data(self.structure, self.init_atoms)

        # make initial gp model and predict forces
        self.update_gp(
            self.init_atoms, dft_frcs, dft_stress=dft_stress, dft_energy=dft_energy
        )
        self.structure.forces = dft_frcs

class C_ASE_OTF(ASE_OTF, C_OTF):
    """On-the-fly training module using ASE MD engine, a subclass of OTF.

    Args:
        atoms (ASE Atoms): the ASE Atoms object for the on-the-fly MD run.
        timestep (:obj:`float`): the timestep in MD. Please use ASE units, e.g. if the
            timestep is 1 fs, then set `timestep = 1 * units.fs`
        number_of_steps (:obj:`int`): the total number of steps for MD.
        dft_calc (ASE Calculator): any ASE calculator is supported,
            e.g. Espresso, VASP etc.
        md_engine (:obj:`str`): the name of MD thermostat, only `VelocityVerlet`,
            `NVTBerendsen`, `NPTBerendsen`, `NPT` and `Langevin`, `NoseHoover`
            are supported.
        md_kwargs (dict): specify the args for MD as a dictionary, the args are
            as required by the ASE MD modules consistent with the `md_engine`.
        update_settings (List[List[dict]]): array containg CHAMP simulation
            parameters to update.
        calculator (:obj:`Calculator <ase.calculators.calculator.Calculator>`): ASE calculator.
            Must have "get_uncertainties" method
          implemented.
        trajectory (ASE Trajectory): default `None`, not recommended,
            currently in experiment.

    The following arguments are for on-the-fly training, the user can also
    refer to :class:`flare.otf.OTF`

    Args:
        prev_pos_init ([type], optional): Previous positions. Defaults
            to None.
        rescale_steps (List[:obj:`int`], optional): List of frames for which the
            velocities of the atoms are rescaled. Defaults to [].
        rescale_temps (List[:obj:`int`], optional): List of rescaled temperatures.
            Defaults to [].
        calculate_energy (:obj:`bool`, optional): If True, the energy of each
            frame is calculated with the GP. Defaults to False.
        write_model (:obj:`int`, optional): If 0, write never. If 1, write at
            end of run. If 2, write after each training and end of run.
            If 3, write after each time atoms are added and end of run.
            If 4, write after each training and end of run, and back up
            after each write.
        std_tolerance_factor (float, optional): Threshold that determines
            when DFT is called. Specifies a multiple of the current noise
            hyperparameter. If the epistemic uncertainty on a force
            component exceeds this value, DFT is called. Defaults to 1.
        skip (:obj:`int`, optional): Number of frames that are skipped when
            dumping to the output file. Defaults to 0.
        init_atoms (List[:obj:`int`], optional): List of atoms from the input
            structure whose local environments and force components are
            used to train the initial GP model. If None is specified, all
            atoms are used to train the initial GP. Defaults to None.
        output_name (:obj:`str`, optional): Name of the output file. Defaults to
            'otf_run'.
        max_atoms_added (:obj:`int`, optional): Number of atoms added each time
            DFT is called. Defaults to 1.
        freeze_hyps (:obj:`int`, optional): Specifies the number of times the
            hyperparameters of the GP are optimized. After this many
            updates to the GP, the hyperparameters are frozen.
            Defaults to 10.
        n_cpus (:obj:`int`, optional): Number of cpus used during training.
            Defaults to 1.
    """

    def __init__(
        self,
        atoms,
        timestep,
        number_of_steps,
        dft_calc,
        md_engine,
        md_kwargs,
        update_settings,
        calculator=None,
        trajectory=None,
        **otf_kwargs
    ):

        self.structure = FLARE_Atoms.from_ase_atoms(atoms)
        if calculator is not None:
            self.structure.calc = calculator
        self.timestep = timestep
        self.md_engine = md_engine
        self.md_kwargs = md_kwargs
        self._kernels = None
        self.update_settings = np.array(update_settings)

        if md_engine == "CustomVerlet":
            MD = CustomVerlet
        else:
            raise NotImplementedError(md_engine + " is not implemented in ASE")

        self.md = MD(
            atoms=self.structure,
            timestep=timestep,
            trajectory=trajectory,
            **md_kwargs,
        )

        force_source = dft_source
        self.flare_calc = self.structure.calc

        # Convert ASE timestep to ps for the output file.
        flare_dt = timestep / (units.fs * 1e3)

        C_OTF.__init__(
            self,
            dt=flare_dt,
            number_of_steps=number_of_steps,
            gp=self.flare_calc.gp_model,
            force_source=force_source,
            dft_loc=dft_calc,
            dft_input=self.structure,
            **otf_kwargs,
        )

        self.flare_name = self.output_name + "_flare.json"
        self.dft_name = self.output_name + "_dft.pickle"
        self.structure_name = self.output_name + "_atoms.json"
        self.checkpt_files = [
            self.checkpt_name,
            self.flare_name,
            self.dft_name,
            self.structure_name,
            self.dft_xyz,
        ]

    def md_step(self, forces=None):
        """Perform an MD step.

        Get new position in molecular dynamics based on the forces predicted by
        FLARE_Calculator or DFT calculator

        Args:
            forces (List[float]): array containing the forces as calculated by CHAMP (or FLARE).
        """

        if self.curr_step in self.update_settings and self.dft_loc.name == "CHAMP":
            sets = self.dft_loc.parameters.settings
            ind = np.where(self.update_settings == self.curr_step)[0][0]
            changes = self.update_settings[ind][1]
            for key, val in changes.items():
                if isinstance(val, dict):
                    for key2, val2 in val.items():
                        setattr(getattr(sets, key), key2, val2)
                else:
                    setattr(sets, key, val)

        # Update previous positions.
        self.structure.prev_positions = np.copy(self.structure.positions)

        # Reset FLARE calculator.
        if self.dft_step:
            self.flare_calc.reset()
            self.structure.calc = self.flare_calc

        # Take MD step.
        # Inside the step() function, get_forces() is called
        self.md.step(forces=forces)
        self.curr_step += 1
