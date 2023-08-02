use std::{collections::HashMap, path::PathBuf};

use log::{debug, error};
use serde::Deserialize;

use crate::{
    config::{self, Config},
    forcefield::FF,
    molecule::Molecule,
    work_queue::WorkQueue,
    Dmat, Dvec,
};

use self::penalty::{Penalty, PenaltyType};

pub(crate) mod penalty;

#[derive(Deserialize)]
struct Metadata {
    dihedrals: Vec<Vec<usize>>,
    grid_spacing: Vec<usize>,
    dihedral_ranges: Option<usize>,
    energy_decrease_thresh: Option<f64>,
    energy_upper_limit: f64,
    torsion_grid_ids: Vec<Vec<isize>>,
}

enum TargetType {
    // TODO some of these are likely shared options. I've seen at least
    // writelevel on both already
    Torsion {
        pdb: String,
        mol2: String,
        coords: String,
        metadata: Metadata,
        ndim: usize,
        freeze_atoms: Vec<usize>,

        /// number of snapshots
        ns: usize,

        /// how much data to write to disk
        writelevel: usize,

        /// harmonic restraint for non-torsion atoms in kcal/mol
        restrain_k: f64,

        /// attenuate the weights as a function of energy
        attenuate: bool,

        /// energy denominator for objective function
        energy_denom: f64,

        /// upper cutoff energy
        energy_upper: f64,
        // TODO might need to read_reference_data, construct engine_args, and
        // engine.
    },
    OptGeo,
}

pub(crate) struct Target {
    // options from Target base class
    /// root directory of the whole project. probably redundant with the same
    /// value on FF
    root: PathBuf,

    /// name of the target
    name: String,

    /// type of target
    typ: TargetType,

    /// relative weight of the target
    weight: f64,

    /// switch for finite difference gradients
    fdgrad: bool,

    /// switch for finite difference hessians
    fdhess: bool,

    /// switch for finite difference gradients + hessian diagonals
    fdhessdiag: bool,

    /// how many seconds to sleep
    sleepy: usize,

    /// parameter types that trigger FD gradient elements. defaults to []. NOTE:
    /// this is obviously not the right type, just picking something small as a
    /// placeholder
    fd1_pids: Vec<()>,

    /// parameter types that trigger FD hessian elements. defaults to [].
    fd2_pids: Vec<()>,

    /// finite difference step size
    pub(crate) h: f64,

    /// whether to make backup files
    backup: bool,

    /// directory to read data from
    rd: PathBuf,

    /// iteration where we turn on zero-gradient skipping
    zerograd: isize,

    /// gradient norm below which we skip
    epsgrad: f64,

    /// "dictionary" of whether to call the derivatives
    pgrad: Vec<bool>,

    /// relative directory of the target
    tgtdir: PathBuf,

    /// temporary working directory. set to temp/target_name. used for storing
    /// temporary variables that don't change through the course of the
    /// optimization
    tempbase: PathBuf,

    tempdir: PathBuf,

    /// the directory in which the simulation is running - this can be updated
    rundir: PathBuf,

    /// need the forcefield (here for now). NOTE: this seems insane. I'm going
    /// to use so much memory if I clone forcefields onto every target
    ff: FF,

    /// counts how often the objective function was computed
    xct: usize,

    /// Counts how often the gradient was computed
    gct: usize,

    /// Counts how often the Hessian was computed
    hct: usize,

    /// Whether to read indicate.log from file when restarting an aborted run.
    read_indicate: bool,

    /// Whether to write indicate.log at every iteration (true for all but remote.)
    write_indicate: bool,

    /// Whether to read objective.p from file when restarting an aborted run.
    read_objective: bool,

    /// Whether to write objective.p at every iteration (true for all but remote.)
    write_objective: bool,

    /// whether the target has been evaluated yet
    evaluated: bool,

    /// whether the previous optimization step was good
    pub(crate) good_step: bool,

    openmm_precision: String,
    openmm_platform: String,
}

impl Target {
    fn stage(&self, mvals: Dvec, order_1: bool, order_2: bool) {
        todo!()
    }

    fn get_x(&self, mvals: &Dvec) -> Extra {
        todo!()
    }

    fn get_g(&self, mvals: &Dvec) -> Extra {
        todo!()
    }

    fn get_h(&self, mvals: &Dvec) -> Extra {
        todo!()
    }
}

#[derive(Default)]
// TODO consider renaming this once I figure out what it's supposed to be. only
// one instance of it is called regularization. other times we're associating
// target names with it
pub struct Regularization {
    w: f64,
    x: f64,
}

impl Regularization {
    pub fn new(w: f64, x: f64) -> Self {
        Self { w, x }
    }
}

pub(crate) struct Extra(f64, Dvec, Dmat);

#[derive(Clone)]
pub(crate) struct ObjMap {
    pub(crate) x0: f64,
    pub(crate) g0: Dvec,
    pub(crate) h0: Dmat,

    /// objective function value
    pub(crate) x: f64,

    /// objective function first derivative
    pub(crate) g: Dvec,

    /// objective function second derivative
    pub(crate) h: Dmat,
}

impl ObjMap {
    pub(crate) fn zeros(size: usize) -> Self {
        Self {
            x0: 0.0,
            g0: Dvec::zeros(size),
            h0: Dmat::zeros(size, size),
            x: 0.0,
            g: Dvec::zeros(size),
            h: Dmat::zeros(size, size),
        }
    }
}

// TODO figure out how to obtain this. going to have to keep track of it on one
// of these structs or maybe just set a global static (...) we're definitely not
// going to examine the stack trace of the program like the python
// version..........
fn in_fd() -> bool {
    false
}

pub struct Objective {
    pub(crate) penalty: Penalty,
    normalize_weights: bool,
    forcefield: FF,
    pub(crate) targets: Vec<Target>,

    /// in Python this is an entry in the ObjDict map. TODO consider making
    /// it an Option depending on how it's used
    obj_map: HashMap<String, Regularization>,

    /// assuming this means total weight
    wtot: f64,
}

impl Objective {
    pub fn new(config: &Config, forcefield: FF) -> Self {
        let ptype = match config.penalty_type.to_lowercase().as_str() {
            "hyp" | "hyper" | "l1" | "hyperbola" | "hyperbolic" => {
                PenaltyType::Hyperbolic
            }
            "para" | "parabola" | "l2" | "quadratic" | "parabolic" => {
                PenaltyType::Parabolic
            }
            "box" => PenaltyType::Box,
            _ => {
                panic!("unrecognized penalty_type: {}", config.penalty_type)
            }
        };
        let penalty = Penalty {
            fadd: config.penalty_additive,
            fmul: config.penalty_multiplicative,
            a: config.penalty_alpha,
            b: config.penalty_hyperbolic_b,
            p: config.penalty_power,
            ptype,
        };
        let mut targets = Vec::new();
        for target in &config.targets {
            let tgtdir = PathBuf::from("targets");
            if !tgtdir.is_dir() {
                panic!(
                    "The targets directory is missing! \
                     Did you finish setting up the target data?"
                );
            }
            let tgtdir = tgtdir.join(&target.name);
            let tempbase = PathBuf::from("temp");
            let tempdir = tempbase.join(&target.name);
            let rundir = tempdir.clone();
            if !config.cont {
                // delete the temporary directory and create a new one. TODO
                // "back up if desired"
                let abstempdir = forcefield.root.join(&tempdir);
                if abstempdir.is_dir() {
                    std::fs::remove_dir_all(&abstempdir).unwrap_or_else(|e| {
                        debug!("failed to remove {abstempdir:?} with {e}")
                    });
                }
                std::fs::create_dir_all(&abstempdir).unwrap_or_else(|e| {
                    debug!("failed to create {abstempdir:?} with {e}")
                });
            }
            let root = forcefield.root.clone();
            let typ = match target.typ {
                config::TargetType::Torsion => {
                    let meta_file = tgtdir.join("metadata.json");
                    let metadata: Metadata = serde_json::from_str(
                        &std::fs::read_to_string(meta_file)
                            .expect("TorsionProfileTarget needs metadata.json"),
                    )
                    .unwrap();
                    let ndim = metadata.dihedrals.len();
                    let freeze_atoms =
                        metadata.dihedrals.iter().flatten().cloned().collect();
                    let pdb = target
                        .pdb
                        .clone()
                        .unwrap_or_else(|| String::from("conf.pdb"));
                    let coords = target
                        .coords
                        .clone()
                        .unwrap_or_else(|| String::from("scan.xyz"));
                    let mol = Molecule::new(
                        root.join(&tgtdir).join(&coords),
                        root.join(&tgtdir).join(&pdb),
                    )
                    .unwrap();
                    let ns = mol.len();
                    TargetType::Torsion {
                        pdb,
                        mol2: target.mol2.clone().unwrap(),
                        coords,
                        metadata,
                        ndim,
                        freeze_atoms,
                        ns,
                        writelevel: target.writelevel,
                        restrain_k: target.restrain_k,
                        attenuate: target.attenuate,
                        energy_denom: target.energy_denom,
                        energy_upper: target.energy_upper,
                    }
                }
                config::TargetType::OptGeo => TargetType::OptGeo,
            };
            targets.push(Target {
                good_step: false,
                h: config.finite_difference_h,
                name: target.name.clone(),
                weight: target.weight,
                evaluated: false,
                openmm_precision: target
                    .openmm_precision
                    .clone()
                    .unwrap_or_else(|| String::from("double")),
                openmm_platform: target
                    .openmm_platform
                    .clone()
                    .unwrap_or_else(|| String::from("Reference")),
                typ,
                root,
                fdgrad: target.fdgrad,
                fdhess: target.fdhess,
                fdhessdiag: target.fdhessdiag,
                sleepy: target.sleepy,
                fd1_pids: Vec::new(),
                fd2_pids: Vec::new(),
                backup: config.backup,
                rd: PathBuf::from(target.read.clone().unwrap_or_default()),
                zerograd: config.zerograd,
                epsgrad: target.epsgrad,
                pgrad: vec![true; forcefield.np],
                tgtdir,
                tempbase,
                tempdir,
                rundir,
                ff: forcefield.clone(),
                xct: 0,
                gct: 0,
                hct: 0,
                read_indicate: true,
                write_indicate: true,
                read_objective: true,
                write_objective: true,
            })
        }
        Self {
            forcefield,
            targets,
            obj_map: HashMap::new(),
            penalty,
            // TODO should be some of weights from targets
            wtot: 1.0,
            normalize_weights: config.normalize_weights,
        }
    }

    pub(crate) fn full(&mut self, vals: Dvec, order: i32) -> ObjMap {
        // TODO borrow here if it doesn't get consumed
        let mut objective = self.target_terms(vals.clone(), order);
        let extra = if self.forcefield.use_pvals {
            self.penalty
                .compute(self.forcefield.create_mvals(vals), &mut objective)
        } else {
            self.penalty.compute(vals, &mut objective)
        };
        objective.x0 = objective.x;
        objective.g0 = objective.g.clone();
        objective.h0 = objective.h.clone();

        if !in_fd() {
            self.obj_map.insert(
                "Regularization".to_owned(),
                Regularization::new(1.0, extra.0),
            );
        }

        objective.x += extra.0;
        objective.g += extra.1;
        objective.h += extra.2;

        objective
    }

    fn target_terms(&mut self, mvals: Dvec, order: i32) -> ObjMap {
        let mut objective = ObjMap::zeros(self.forcefield.np);
        for tgt in &self.targets {
            // TODO consider borrowing
            tgt.stage(mvals.clone(), order >= 1, order >= 2);
        }

        // TODO might be supposed to use an existing work queue. see
        // nifty/getWorkQueue at some point. this might even be a psqs
        // situation.
        let wq = WorkQueue::new();
        for tgt in self.targets.iter_mut() {
            let ans: Extra = match order {
                0 => tgt.get_x(&mvals),
                1 => tgt.get_g(&mvals),
                2 => tgt.get_h(&mvals),
                _ => unimplemented!(),
            };
            if !in_fd() {
                self.obj_map.insert(
                    tgt.name.clone(),
                    Regularization::new(tgt.weight / self.wtot, ans.0),
                );
            }
            objective.x += ans.0 * tgt.weight / self.wtot;
            objective.g += ans.1 * tgt.weight / self.wtot;
            objective.h += ans.2 * tgt.weight / self.wtot;
        }

        for tgt in self.targets.iter_mut() {
            tgt.evaluated = true;
        }

        // prevent exact zeros on Hessian diagonal
        for i in 0..self.forcefield.np {
            if objective.h[(i, i)] == 0.0 {
                objective.h[(i, i)] = 1.0;
            }
        }

        objective
    }
}
