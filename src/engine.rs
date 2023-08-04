/// this corresponds to the SMIRNOFF class in ForceBalance, which is a subclass
/// of OpenMM. basically it's just an interface to openmm
pub(crate) struct Engine {
    restraint_frc_index: Option<usize>,
}

impl Engine {
    /// initialization actually just seems to consist of setting up directories.
    /// not sure what options actually need to be stored at the moment
    pub(crate) fn new() -> Self {
        Self {
            restraint_frc_index: None,
        }
    }

    /// optimize the geometry and align the optimized geometry to the starting
    /// geometry. `shot` is the snapshot number to be minimized, presumably
    /// `align` handles the "and align ..." part if set to `true`. returns the
    /// energy, the rmsd in Angstrom, and a third value to be determined
    pub(crate) fn optimize(&self, shot: usize, align: bool) -> (f64, f64, f64) {
        // optional argument in Python, this is the default value
        const CRIT: f64 = 1e-4;
        // TODO if crit stays a constant, steps is also a constant 4 lol
        let steps = f64::max(1.0, -CRIT.log10()) as usize;
        self.update_simulation();
        self.set_positions(shot);

        // in-line self.set_restraint_positions(shot)
        if let Some(idx) = self.restraint_frc_index {
            // TODO none in the case I'm using, but here's the code I drafted:

            // let xyz = self.ref_mol.xyzs[shot] / 10.0; // convert to nm
            // let mut frc = self.simulation.system.get_force(idx);
            // for (i, j) in self.real_atom_idxs.iter().enumerate() {
            //     frc.set_particle_parameters(i, j, xyz[i]);
            // }
            // frc.update_parameters_in_context(&self.simulation.context);
            todo!()
        }

        let x0 = self.simulation.context.get_state().get_positions()
            [self.real_atom_idxs];
        todo!();
    }

    /// create the simulation object, or update the force field parameters in
    /// the existing simulation object. this should be run when we write a new
    /// force field XML file.
    fn update_simulation(&self) {
        todo!()
    }

    fn set_positions(&self, shot: usize) {
        todo!()
    }
}
