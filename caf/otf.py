from time import time
from copy import deepcopy
import numpy as np

from flare.otf import OTF
from flare.ase.otf import ASE_OTF
from flare.ase.atoms import FLARE_Atoms
import flare.ase.dft as dft_source
from flare.utils.learner import is_std_in_bound

from ase import units

from caf.verlet import CustomVerlet


class C_OTF(OTF):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def run(self):
        """
        Performs an on-the-fly training run.

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
            # TODO: Reinstate velocity rescaling.
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
        """
        Get new position in molecular dynamics based on the forces predicted by
        FLARE_Calculator or DFT calculator
        """

        if self.curr_step in self.update_settings:
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