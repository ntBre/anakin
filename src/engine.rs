use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

use log::debug;
use openff_toolkit::smirnoff::ForceField;
use openff_toolkit::topology::molecule::Molecule as OffMolecule;
use openff_toolkit::topology::Topology as OffTopology;
use openmm::pdb::modeller::Modeller;
use openmm::pdb::PDBFile;
use openmm::{
    integrators::{Integrator, Verlet},
    platform::Platform,
    system::System,
    topology::Vec3,
    Simulation,
};

use crate::{
    forcefield::FF, molecule::Molecule, objective::Target,
    objective::TargetType, Dvec,
};

/// box for periodic boundary conditions
struct PBox;

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

    platname: String,

    precision: String,

    /// initially empty, set to Some in Self::set_opts
    platform: Option<Platform>,

    ism: Option<()>,

    // these are initially None/empty and then initialized in `read_src`
    mol: Option<Molecule>,
    ref_mol: Option<Molecule>,
    mol2: Vec<String>,
    abspdb: Option<PathBuf>,
    pdb: Option<PDBFile>,
    off_topology: Option<OffTopology>,
    restrain_k: f64,
    freeze_atoms: Vec<usize>,
    pbc: bool,
    has_virtual_sites: bool,
    xyz_omms: Vec<(Vec<Vec3>, Option<PBox>)>,
    modeller: Option<Modeller>,
    system: Option<System>,
}

impl Engine<Verlet> {
    /// initialization actually just seems to consist of setting up directories.
    /// not sure what options actually need to be stored at the moment. NOTE:
    /// target is supposed to be kwargs but we obviously don't have kwargs
    pub(crate) fn new(target: Target) -> Self {
        let pdb = PDBFile::new("testfiles/test.pdb");
        let mut m = Modeller::new(pdb.topology, pdb.positions);
        let ff = ForceField::load("testfiles/force-field.offxml").unwrap();
        // let topology = m.topology;
        // TODO crashing here because "Cannot create a context for a System with
        // no particles," suggesting that ff.create_system is the important
        // place to fix this

        debug!("creating engine");
        // TODO this is supposed to be an OFF topology
        // let system = ff.create_system(topology.clone());
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
            platname: Default::default(),
            precision: Default::default(),
            platform: Default::default(),
            ism: Default::default(),
            mol: None,
            ref_mol: None,
            mol2: Vec::new(),
            abspdb: None,
            pdb: None,
            off_topology: None,
            restrain_k: 1.0,
            freeze_atoms: Vec::new(),
            pbc: false,
            has_virtual_sites: false,
            xyz_omms: Vec::new(),
            modeller: None,
            system: None,
        };
        // Step 1: set up options
        ret.set_opts();
        let TargetType::Torsion {
            pdb,
            mol2,
            coords,
            metadata,
            ndim,
            freeze_atoms,
            ns,
            writelevel,
            restrain_k,
            attenuate,
            energy_denom,
            energy_upper,
            wts,
            eqm,
            mol,
        } = &ret.target.typ
        else {
            panic!("can't handle this type yet");
        };
        // Step 2: Read data from the source directory
        ret.read_src(&pdb.clone(), mol.clone(), vec![mol2.clone()]);
        // Step 3: Prepare the temporary directory
        ret.prepare();

        debug!("finished creating engine");

        ret
    }

    /// set OpenMM-specific options
    fn set_opts(&mut self) {
        // TODO these are supposed to come from self.target, but I don't have
        // those fields on Target. just taking them as defaults
        self.platname = "Reference".to_owned();
        self.precision = "double".to_owned();

        // TODO could check that platname is a registered OpenMM Platform, but I
        // don't really care. just assume it is registered.

        self.platform = Some(Platform::by_name(&self.platname));
        // TODO check if platname is CUDA or OpenCL, but we're sticking with
        // Reference for now
        self.ism = None;
    }

    /// read files from the source directory. Provide a molecule object or a
    /// coordinate file. Add an optional PDB file for residues, atom names, etc.
    /// `pdb` is a PDB file, `mol` is a [Molecule], `mol2` is a list of .mol2
    /// files. In forcebalance these are all passed as kwargs
    fn read_src(&mut self, pdb: &str, mol: Molecule, mol2: Vec<String>) {
        let pdbfnm = self.target.root.join(&self.target.tgtdir).join(pdb);
        assert!(pdbfnm.exists(), "PDB file {:?} does not exist", pdbfnm);
        self.mol = Some(mol);
        self.mol2 = mol2;
        for fnm in &self.mol2 {
            let p = self.target.root.join(&self.target.tgtdir).join(&fnm);
            assert!(p.exists(), "file {p:?} does not exist");
        }
        self.abspdb = Some(pdbfnm.canonicalize().unwrap());
        let mpdb = Molecule::from_path(pdbfnm);

        // TODO might have to set self.mol.Data.[chain, atomname, resid,
        // resname, elem]

        self.ref_mol = self.mol.clone();
    }

    /// optimize the geometry and align the optimized geometry to the starting
    /// geometry. `shot` is the snapshot number to be minimized, presumably
    /// `align` handles the "and align ..." part if set to `true`. returns the
    /// energy, the rmsd in Angstrom, and a third value to be determined
    pub(crate) fn optimize(
        &mut self,
        shot: usize,
        align: bool,
    ) -> (f64, f64, f64) {
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
    fn update_simulation(&mut self) {
        debug!("updating simulation");
        let system = self
            .ff
            .openff_forcefield
            .create_system(self.off_topology.clone().unwrap());
        debug!("created system");
        self.system = Some(system);
        if let Some(sim) = &mut self.simulation {
            update_simulation_parameters(sim);
        } else {
            debug!("creating simulation");
            self.create_simulation()
        }
        debug!("finished updating simulation");
    }

    fn set_positions(&self, shot: usize) {
        todo!()
    }

    /// Prepare the calculation. Note that we don't create the Simulation object
    /// yet because that may depend on MD integrator parameters, thermostat,
    /// barostat, etc. This is mostly copied from OpenMM.prepare(), but we are
    /// calling ForceField from the OpenFF toolkit and ignoring AMOEBA stuff.
    fn prepare(&mut self) {
        if let Some(pdb) = &self.abspdb {
            let pdb = PDBFile::new(pdb);
            self.pdb = Some(pdb);
        }
        // load mol2 files for smirnoff topology
        let mut openff_mols = Vec::new();
        for fnm in &self.mol2 {
            let fnm = self.root.join(&self.target.tgtdir).join(fnm);
            openff_mols.push(OffMolecule::from_file(fnm));
        }
        self.off_topology = Some(OffTopology::from_openmm(
            &self.pdb.as_ref().unwrap().topology,
            openff_mols,
        ));
        self.restrain_k = 1.0;
        // TODO where did this come from?? I know it's from kwargs, but where
        // did it get determined
        self.freeze_atoms = vec![5, 6, 7, 8];

        // TODO add this to self.mmopts
        // self.ff.rigid_water;
        self.ff.make(&Dvec::zeros(self.ff.np), &self.tempdir);

        // periodic boundary conditions
        self.pbc = false;

        // Apply the FF parameters to the system. Currently this is the only way
        // to determine if the FF will apply virtual sites to the system.
        let interchange = self
            .ff
            .openff_forcefield
            .create_interchange(self.off_topology.as_ref().unwrap().clone());

        self.has_virtual_sites = !interchange.virtual_sites.is_empty();

        // handled by default
        // NOTE: no outer loop because I only have one xyz for now
        let mol = self.mol.as_ref().unwrap();
        let mut xyz_omm = Vec::new();
        for x in mol.xyzs.row_iter() {
            xyz_omm.push(Vec3::new(x[0], x[1], x[2]));
        }
        // TODO if I support pbc, obtain the box here
        self.xyz_omms.push((xyz_omm, None));
        // let openmm_topology = interchange.to_openmm_topology();
        // TODO might need some unit conversion / appending placehold positions
        // for v sites
        let openmm_positions = self.pdb.as_ref().unwrap().positions.clone();
        let openmm_topology = self.pdb.as_ref().unwrap().topology.clone();
        self.modeller = Some(Modeller::new(openmm_topology, openmm_positions));

        // build a topology and atom lists. BRW: why do we do this over and over
        // in different formats?
        let top = self.modeller.as_ref().unwrap().get_topology();
        let atoms = top.atoms();

        // TODO some weird atom_lists stuff in smirnoffio.py, skipping for now
        for f in &self.ff.fnms {
            let path = self.root.join(&self.tempdir).join(f);
            std::fs::remove_file(path)
                .unwrap_or_else(|e| panic!("failed to remove {f} with {e}"));
        }
    }

    fn create_simulation(&mut self) {
        let top = self.modeller.clone().unwrap().topology;
        let sys = self.system.clone().unwrap();
        let int = Verlet::new(1.0);
        let simulation = Simulation::new(top, sys, int);
        self.simulation = Some(simulation);
    }
}

fn update_simulation_parameters(sim: &mut Simulation<Verlet>) {
    dbg!("in sim params");
    todo!()
}
