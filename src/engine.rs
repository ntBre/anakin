use std::{collections::HashMap, path::PathBuf};

use openff_toolkit::smirnoff::ForceField;
use openmm::{
    integrators::{Integrator, Verlet},
    Modeller, PDBFile, Simulation,
};

use crate::{forcefield::FF, objective::Target};

/// this corresponds to the SMIRNOFF class in ForceBalance, which is a subclass
/// of OpenMM. basically it's just an interface to openmm
pub(crate) struct Engine<I>
where
    I: Integrator,
{
    valkwd: Vec<String>,
    name: String,
    target: Target,
    root: PathBuf,
    srcdir: PathBuf,
    tempdir: PathBuf,
    ff: FF,
    restraint_frc_index: Option<usize>,

    /// the actual simulation. initially none, set up later
    simulation: Option<Simulation<I>>,

    /// initially empty, generated in prepare
    real_atom_idxs: Vec<usize>,
}

impl Engine<Verlet> {
    /// initialization actually just seems to consist of setting up directories.
    /// not sure what options actually need to be stored at the moment. NOTE:
    /// target is supposed to be kwargs but we obviously don't have kwargs
    pub(crate) fn new(target: Target) -> Self {
        // let pdb = PDBFile::new("testfiles/test.pdb");
        // let mut m = Modeller::new(pdb.topology, pdb.positions);
        // let ff = ForceField::load("testfiles/force-field.offxml").unwrap();
        // let topology = m.topology;
        // let system = ff.create_system(topology);
        // let integrator = Verlet::new(1.0);
        // let simulation = Simulation::new(topology, system, integrator);
        let mut ret = Self {
            restraint_frc_index: None,
            simulation: None,
            valkwd: vec![
                // these are from smirnoffio
                "ffxml".to_string(),
                "pdb".to_string(),
                "mol2".to_string(),
                "platname".to_string(),
                "precision".to_string(),
                "mmopts".to_string(),
                "vsite_bonds".to_string(),
                "implicit_solvent".to_string(),
                "restrain_k".to_string(),
                "freeze_atoms".to_string(),
                // these are from engine
                "mol".to_string(),
                "coords".to_string(),
                "name".to_string(),
                "target".to_string(),
                "pbc".to_string(),
                "FF".to_string(),
                "nonbonded_method".to_string(),
                "nonbonded_cutoff".to_string(),
            ],
            name: "openmm".to_string(),
            srcdir: target.root.join(&target.tgtdir),
            tempdir: target.root.join(&target.tempdir),
            root: target.root.clone(),
            ff: target.ff.clone(),
            target,
            real_atom_idxs: Vec::new(),
        };
        ret.set_opts();
        ret
    }

    fn set_opts(&mut self) {
        todo!()
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

        use openmm::state::DataType as D;
        let positions = self
            .simulation
            .as_ref()
            .unwrap()
            .context
            .get_state(D::Positions)
            .get_positions();

        let x0 = positions
            .into_iter()
            .enumerate()
            .filter(|(p, _)| self.real_atom_idxs.contains(p));
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
