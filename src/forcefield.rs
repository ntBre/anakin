use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

use openff_toolkit::smirnoff::{ForceField, Unit};

use crate::{config::Config, Dmat, Dvec};

#[derive(Clone)]
struct Digraph;

#[derive(Clone)]
pub struct FF {
    // ============================
    // options from the config file
    // ============================
    /// the root directory of the project
    root: PathBuf,

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
        // TODO clearly this is supposed to have some length
        let mut pvals0 = Vec::new();
        // TODO real value here
        let np = 0;
        let rs = Dvec::from_element(pvals0.len(), 1.0);
        let qmat2 = Dmat::identity(np, np);
        // TODO qmat2 needs to be set for real
        let transmat = qmat2 * Dmat::from_diagonal(&rs);
        let tmi = transmat.transpose();
        let fnms = config.forcefield.clone();

        // combination of addff and addff_xml in python version
        let mut forcefields = Vec::new();
        let mut offxml = None;
        let mut openff_forcefield = None;
        let mut ffdata_isxml = HashMap::new();
        let ffdir = Path::new(&config.ffdir);
        let mut bonds_to_optimize = Vec::new();
        let mut angles_to_optimize = Vec::new();
        let mut propers_to_optimize = Vec::new();
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

        dbg!(pvals0.len());

        let rs = rsmake(
            bonds_to_optimize,
            &openff_forcefield,
            angles_to_optimize,
            propers_to_optimize,
            config,
            &pvals0,
        );

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
            tm: todo!(),
            tmi,
            excision: Vec::new(),
            np,
            pvals0: Dvec::from(pvals0),
            atomnames: todo!(),
            bonds_to_optimize: todo!(),
            angles_to_optimize: todo!(),
            propers_to_optimize: todo!(),
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
        let mut pvals = (&self.tmi * mvals) + &self.pvals0;

        // redirect parameters (for the fusion penalty function)
        // TODO another redirect
        // for p in self.redirect

        pvals
    }

    pub(crate) fn make_redirect(&self, mvals: Dvec) {
        todo!()
    }
}

fn rsmake(
    bonds_to_optimize: Vec<Param>,
    openff_forcefield: &Option<ForceField>,
    angles_to_optimize: Vec<Param>,
    propers_to_optimize: Vec<Param>,
    config: &Config,
    pvals0: &Vec<f64>,
) -> Dvec {
    // extracting the maximum values of the parameters for computing the
    // rescaling factors.
    let mut max_bond = HashMap::new();
    for (i, param) in bonds_to_optimize.iter().enumerate() {
        if let Param::Opt { inner } = param {
            for (_, typ, _unit) in inner {
                let val = openff_forcefield.as_ref().unwrap().bonds[i]
                    .as_hash(typ)
                    .unwrap()
                    .value
                    .abs();
                let cur = max_bond.entry(typ.clone()).or_insert(val);
                if val > *cur {
                    *cur = val;
                }
            }
        }
    }

    let mut max_angle = HashMap::new();
    for (i, param) in angles_to_optimize.iter().enumerate() {
        if let Param::Opt { inner } = param {
            for (_, typ, _unit) in inner {
                let val = openff_forcefield.as_ref().unwrap().angles[i]
                    .as_hash(typ)
                    .unwrap()
                    .value
                    .abs();
                let cur = max_angle.entry(typ.clone()).or_insert(val);
                if val > *cur {
                    *cur = val;
                }
            }
        }
    }

    let mut max_proper = HashMap::new();
    for (i, param) in propers_to_optimize.iter().enumerate() {
        if let Param::Opt { inner } = param {
            for (_, typ, _unit) in inner {
                let val = openff_forcefield.as_ref().unwrap().proper_torsions
                    [i]
                    .as_hash(typ)
                    .unwrap()
                    .value
                    .abs();
                let cur = max_proper.entry(typ.clone()).or_insert(val);
                if val > *cur {
                    *cur = val;
                }
            }
        }
    }

    // overwrite from priors
    if let Some(map) = &config.priors {
        if let Some(bonds) = map.get("bonds") {
            for (key, value) in bonds {
                max_bond.insert(key.to_string(), *value);
            }
        }
        if let Some(angles) = map.get("angles") {
            for (key, value) in angles {
                max_angle.insert(key.to_string(), *value);
            }
        }
        if let Some(propers) = map.get("propertorsions") {
            for (key, value) in propers {
                if key == "k" {
                    for key in ["k1", "k2", "k3", "k4", "k5", "k6"] {
                        max_proper.insert(key.to_string(), *value);
                    }
                }
            }
        }
    }

    let mut rs = Dvec::from_element(pvals0.len(), 1.0);
    'outer: for i in 0..pvals0.len() {
        for bond in &bonds_to_optimize {
            if let Param::Opt { inner } = bond {
                for (x, typ, _) in inner {
                    if *x == i {
                        rs[i] = *max_bond.get(typ).unwrap();
                        continue 'outer;
                    }
                }
            }
        }
        for angle in &angles_to_optimize {
            if let Param::Opt { inner } = angle {
                for (x, typ, _) in inner {
                    if *x == i {
                        rs[i] = *max_angle.get(typ).unwrap();
                        continue 'outer;
                    }
                }
            }
        }
        for proper in &propers_to_optimize {
            if let Param::Opt { inner } = proper {
                for (x, typ, _) in inner {
                    if *x == i {
                        rs[i] = *max_proper.get(typ).unwrap();
                        continue 'outer;
                    }
                }
            }
        }
    }
    rs
}

// TODO might be able to factor these out with the use of parameter_handlers

fn add_bonds(
    ff: &ForceField,
    pvals0: &mut Vec<f64>,
    bonds_to_optimize: &mut Vec<Param>,
) {
    for (i, bond) in (&ff.bonds).into_iter().enumerate() {
        if let Some(to_optimize) = &bond.parameterize {
            let mut inner = Vec::new();
            for param in to_optimize.split(',').map(str::trim) {
                let p = bond.as_hash(param).unwrap();
                inner.push((pvals0.len(), param.to_owned(), p.unit.clone()));
                pvals0.push(p.value);
            }
            bonds_to_optimize.push(Param::Opt { inner });
        } else {
            bonds_to_optimize.push(Param::NoOpt);
        }
    }
}

fn add_angles(
    ff: &ForceField,
    pvals0: &mut Vec<f64>,
    angles_to_optimize: &mut Vec<Param>,
) {
    for (i, angle) in (&ff.angles).into_iter().enumerate() {
        if let Some(to_optimize) = &angle.parameterize {
            let mut inner = Vec::new();
            for param in to_optimize.split(',').map(str::trim) {
                let p = angle.as_hash(param).unwrap();
                inner.push((pvals0.len(), param.to_owned(), p.unit.clone()));
                pvals0.push(p.value);
            }
            angles_to_optimize.push(Param::Opt { inner });
        } else {
            angles_to_optimize.push(Param::NoOpt);
        }
    }
}

fn add_propers(
    ff: &ForceField,
    pvals0: &mut Vec<f64>,
    propers_to_optimize: &mut Vec<Param>,
) {
    for (i, proper) in (&ff.proper_torsions).into_iter().enumerate() {
        if let Some(to_optimize) = &proper.parameterize {
            let mut inner = Vec::new();
            for param in to_optimize.split(',').map(str::trim) {
                let p = proper.as_hash(param).unwrap();
                inner.push((pvals0.len(), param.to_owned(), p.unit.clone()));
                pvals0.push(p.value);
            }
            propers_to_optimize.push(Param::Opt { inner });
        } else {
            propers_to_optimize.push(Param::NoOpt);
        }
    }
}
