use std::{
    collections::{HashMap, HashSet},
    path::{Path, PathBuf},
};

use openff_toolkit::smirnoff::{ForceField, Unit};

use crate::{config::Config, utils::orthogonalize, Dmat, Dvec};

mod utils;

use utils::*;

#[derive(Clone)]
struct Digraph;

#[derive(Clone)]
pub struct FF {
    // ============================
    // options from the config file
    // ============================
    /// the root directory of the project
    pub(crate) root: PathBuf,

    /// force field names
    fnms: Vec<String>,

    /// force field directory relative to the project directory
    ffdir: String,

    /// priors given by the user
    priors: HashMap<String, HashMap<String, f64>>,

    /// whether to constrain the charges
    constrain_charge: bool,

    /// whether to use the logarithmic map
    logarithmic_map: bool,

    /// amoeba polarization type
    amoeba_pol: String,

    /// AMOEBA mutual dipole convergence tolerance.
    amoeba_eps: f64,

    /// Switch for rigid water molecules
    rigid_water: bool,

    /// Switch for constraining bonds involving hydrogen
    constrain_h: bool,

    /// Bypass the transformation matrix and use the physical parameters
    /// directly used in: Creating the force field; advanced usage, be careful.
    /// TODO taken from the config file
    pub(crate) use_pvals: bool,

    /// Allow duplicate parameter names (internally construct unique names)
    duplicate_pnames: bool,

    // ==================
    // variables set here
    // ==================
    /// the content of all force field files stored in memory
    ffdata: Vec<ForceField>,

    /// parallel vector to `self.ffdata` saying whether each entry is in xml
    /// format
    ffdata_isxml: HashMap<String, bool>,

    /// smirnoff force field name. TODO luckily for me, forcebalance appears to
    /// restrict you to a single smirnoff forcefield, so `ffdata`,
    /// `ffdata_isxml`, `offxml`, and `openff_forcefield` should all be combined
    /// into a single field probably called `forcefield` (as I had
    /// originally...)
    offxml: String,

    openff_forcefield: ForceField,

    /// improved representation of pfields as a networkx graph. NOTE: please,
    /// please, please let this be unused.
    ptree: Digraph,

    /// list of rescaling factors
    rs: Dvec,

    /// transformation matrix for mathematical -> physical parameters
    tm: Dmat,

    /// transformation matrix for mathematical -> physical parameters
    tmi: Dmat,

    /// indices to be excluded from the Hessian update
    pub(crate) excision: Vec<usize>,

    /// The total number of parameters.
    pub(crate) np: usize,

    /// initial value of physical parameters
    pvals0: Dvec,

    /// a list of atom names
    atomnames: Vec<String>,

    /// each entry in these vecs corresponds to one parameter of the associated
    /// type in `self.openff_forcefield`. to construct a new force field with
    /// the optimized values, loop over these vectors; for a NoOpt, take the
    /// value of the parameter straight from the original force field; for an
    /// Opt, loop over its `inner` entries, taking the values from the indices
    /// into pvals and the field names and units from the Strings
    bonds_to_optimize: Vec<Param>,
    angles_to_optimize: Vec<Param>,
    propers_to_optimize: Vec<Param>,
}

#[derive(Clone, Debug)]
enum Param {
    Opt {
        inner: Vec<(usize, String, Unit)>,
    },
    /// not being optimized, pass it along straight from the original force
    /// field
    NoOpt,
}

impl FF {
    pub fn new(config: &Config) -> Self {
        // combination of addff and addff_xml in python version
        let mut forcefields = Vec::new();
        let mut offxml = None;
        let mut openff_forcefield = None;
        let mut ffdata_isxml = HashMap::new();
        let ffdir = Path::new(&config.ffdir);
        let mut bonds_to_optimize = Vec::new();
        let mut angles_to_optimize = Vec::new();
        let mut propers_to_optimize = Vec::new();
        let mut pvals0 = Vec::new();
        let fnms = config.forcefield.clone();
        for fnm in &fnms {
            offxml = Some(fnm.clone());
            let ff_file = &ffdir.join(fnm);
            let ff = match ForceField::load(ff_file) {
                Ok(ff) => ff,
                Err(e) => panic!("failed to open {ff_file:?} with {e}"),
            };
            openff_forcefield = Some(ff.clone());
            ffdata_isxml.insert(fnm.clone(), true);

            // idea here should really be to partition the parameters into those
            // we want to optimize and those to leave alone. we just need some
            // representation that lets us change the values we need to change
            // and write them back to a force field file along with the ones we
            // shouldn't change. the python code seems to be editing the strings
            // in place, but it makes a lot more sense to use the structs I
            // already have defined from SMIRNOFF
            add_bonds(&ff, &mut pvals0, &mut bonds_to_optimize);
            add_angles(&ff, &mut pvals0, &mut angles_to_optimize);
            add_propers(&ff, &mut pvals0, &mut propers_to_optimize);

            forcefields.push(ff);
        }

        let rs = rsmake(
            &bonds_to_optimize,
            &openff_forcefield,
            &angles_to_optimize,
            &propers_to_optimize,
            config,
            &pvals0,
        );

        let np = pvals0.len();

        let (transmat, excision, tmi) = mktransmat(np, &rs);

        // TODO might have to overwrite with physically-motivated values, but
        // they aren't triggered for the force field I'm testing on

        // TODO there is supposed to be more processing of the force fields here
        // to extract the useful parts. really we probably don't even need the
        // original ForceField structs
        FF {
            root: std::env::current_dir()
                .expect("somehow failed to read the working directory"),
            fnms,
            ffdir: config.ffdir.clone(),
            priors: config.priors.clone().unwrap_or_default(),
            constrain_charge: config.constrain_charge,
            logarithmic_map: config.logarithmic_map,
            amoeba_pol: config.amoeba_pol.clone().unwrap_or_default(),
            amoeba_eps: config.amoeba_eps.unwrap_or_default(),
            rigid_water: config.rigid_water,
            constrain_h: config.constrain_h,
            use_pvals: config.use_pvals,
            duplicate_pnames: config.duplicate_pnames,
            ffdata: forcefields,
            ffdata_isxml,
            offxml: offxml.unwrap(),
            openff_forcefield: openff_forcefield.unwrap(),
            ptree: Digraph,
            rs,
            tm: transmat,
            tmi,
            excision,
            np,
            pvals0: Dvec::from(pvals0),
            atomnames: Vec::new(),
            bonds_to_optimize,
            angles_to_optimize,
            propers_to_optimize,
        }
    }

    pub(crate) fn create_mvals(&self, vals: Dvec) -> Dvec {
        todo!()
    }

    pub(crate) fn create_pvals(&self, mvals: Dvec) -> Dvec {
        // TODO potentially experimental self.redirect feature

        // TODO potentially handle self.logarithmic_map option disabled by
        // default

        // python version is flat(np.dot(self.tmi, col(mvals)). I might have to
        // flatten the result of the multiplication or something here

        // redirect parameters (for the fusion penalty function)
        // TODO another redirect
        // for p in self.redirect

        // let mut pvals = (&self.tmi * mvals) + &self.pvals0;
        (&self.tmi * mvals) + &self.pvals0
    }

    pub(crate) fn make_redirect(&self, mvals: Dvec) {
        todo!()
    }

    /// this is an insane method... superficially it just makes pvals from
    /// `vals`, but it has the side effect of writing an xml file containing the
    /// FF data!
    pub(crate) fn make(&self, vals: &Dvec, dir: impl AsRef<Path>) -> Dvec {
        let pvals = if self.use_pvals {
            vals.clone()
        } else {
            // TODO ref? I think that might just move the clone deep into
            // nalgebra instead of here so probably no difference
            self.create_pvals(vals.clone())
        };

        assert_eq!(self.ffdata.len(), 1);

        let mut newffdata = self.ffdata.clone().pop().unwrap();

        // lookin pretty wet again
        for (i, bond) in self.bonds_to_optimize.iter().enumerate() {
            if let Param::Opt { inner } = bond {
                // assume the unit has not been changed
                for (p, field, _unit) in inner {
                    let mut field =
                        newffdata.bonds[i].as_hash_mut(field).unwrap();
                    field.value = pvals[*p];
                }
            }
        }

        for (i, angle) in self.angles_to_optimize.iter().enumerate() {
            if let Param::Opt { inner } = angle {
                // assume the unit has not been changed
                for (p, field, _unit) in inner {
                    let mut field =
                        newffdata.angles[i].as_hash_mut(field).unwrap();
                    field.value = pvals[*p];
                }
            }
        }

        for (i, proper) in self.propers_to_optimize.iter().enumerate() {
            if let Param::Opt { inner } = proper {
                // assume the unit has not been changed
                for (p, field, _unit) in inner {
                    let mut field = newffdata.proper_torsions[i]
                        .as_hash_mut(field)
                        .unwrap();
                    field.value = pvals[*p];
                }
            }
        }

        let fnm = self.fnms.first().unwrap();
        let path = dir.as_ref().join(fnm);
        std::fs::write(dbg!(path), newffdata.to_xml().unwrap()).unwrap();

        pvals
    }
}

fn mktransmat(np: usize, rs: &Dvec) -> (Dmat, Vec<usize>, Dmat) {
    let mut qmap: Vec<usize> = Vec::new();
    let mut qid: Vec<Vec<usize>> = Vec::new();
    let mut qid2: Vec<Vec<usize>> = Vec::new();
    let qnr = 1;
    let mut qmat2 = Dmat::identity(np, np);

    // he really loves these nested functions. TODO pull these out as
    // standalone functions
    let insert_mat = |qtrans2: Dmat, qmap: Vec<usize>| {
        let mut x = 0;
        for i in 0..np {
            if qmap.contains(&i) {
                let mut y = 0;
                for (y, &j) in qmap.iter().enumerate() {
                    qmat2[(i, j)] = qtrans2[(x, y)];
                }
                x += 1;
            }
        }
    };

    let build_qtrans2 = |tq, qid: Vec<Vec<usize>>, qmap: Vec<usize>| {
        let nq = qmap.len();
        // tq = Total number of atomic charges that are being optimized on
        // the molecule NOTE: This may be greater than the number of charge
        // parameters (nq) The reason for the "one" here is because LP
        // wanted to have multiple charge constraints at some point in the
        // future

        const ONE: usize = 1;

        let cons0 = Dmat::from_element(ONE, tq, 1.0);
        // in python, the first dimension is cons0.shape[0], but we just
        // hard-coded that to 1 right??
        let mut cons = Dmat::zeros(ONE, nq);
        let mut qtrans2 = Dmat::identity(nq, nq);

        // again, this uses cons.shape[0], but we set it to 1 !
        for i in 0..ONE {
            // cons.shape[1] but we know what that is from a few lines back
            for j in 0..nq {
                // BRW: not sure if this is still relevant, but copying from
                // Python just in case:
                //
                // Each element of qid is a list that points to atom
                // indices. LPW: This code is breaking when we're not
                // optimizing ALL the charges Replace cons0[i][k-1] with all
                // ones cons[i][j] = sum([cons0[i][k-1] for k in qid[j]])
                cons[(i, j)] = qid[j].len() as f64;
            }
            let n = cons.row(i).norm();
            let mut row = cons.row_mut(i);
            row /= n;
            for j in 0..i {
                let o = orthogonalize(
                    cons.row(i).transpose(),
                    cons.row(j).transpose(),
                )
                .transpose();
                cons.set_row(i, &o);
            }
            let mut row = qtrans2.row_mut(i);
            row.fill(0.0);
            for j in 0..nq - i - 1 {
                let o = orthogonalize(
                    qtrans2.row(i + j + 1).transpose(),
                    cons.row(i).transpose(),
                )
                .transpose();
                qtrans2.set_row(i + j + 1, &o);
            }
        }

        qtrans2
    };

    // TODO could need to build a charge constraint for each molecule. we
    // don't enter this for the case I'm looking at, maybe not SMIRNOFF in
    // general?

    let mut transmat = &qmat2 * Dmat::from_diagonal(rs);
    let mut transmat_ns = transmat.clone();
    let mut excision = HashSet::new();
    for i in 0..np {
        if transmat_ns[(i, i)].abs() < 1e-8 {
            excision.insert(i);
            transmat_ns[(i, i)] += 1.0;
        }
    }
    let excision = excision.into_iter().collect();
    for &i in &excision {
        transmat.set_column(i, &Dvec::zeros(np));
    }
    let tmi = transmat.transpose();
    (transmat, excision, tmi)
}
