use std::{collections::HashMap, path::PathBuf};

use openff_toolkit::smirnoff::ForceField;

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
    priors: HashMap<String, String>,

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

    /// the mapping of interaction type -> parameter number
    map: HashMap<String, String>,

    /// the listing of parameter number -> interaction type
    plist: Vec<String>,

    /// the listing of parameter number -> atoms involved
    patoms: Vec<String>,

    /// A list where pfields[i] = [pid, 'file', line, field, mult, cmd],
    /// basically a new way to modify force field files; when we modify the
    /// force field file, we go to the specific line/field in a given file and
    /// change the number. NOTE: is this genius or insanity??
    pfields: Vec<String>,

    /// improved representation of pfields as a networkx graph. NOTE: please,
    /// please, please let this be unused.
    ptree: Digraph,

    /// unit strings that might appear in offxml file
    offxml_unit_strs: HashMap<String, String>,

    /// list of rescaling factors
    rs: Vec<f64>,

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
}

impl FF {
    pub fn new(config: &Config) -> Self {
        // TODO clearly this is supposed to have some length
        let pvals0 = Dvec::zeros(0);
        // TODO real value here
        let np = 0;
        let rs = Dvec::from_element(pvals0.len(), 1.0);
        let qmat2 = Dmat::identity(np, np);
        // TODO qmat2 needs to be set for real
        let transmat = qmat2 * Dmat::from_diagonal(&rs);
        let tmi = transmat.transpose();
        let fnms = config.forcefield.clone();

        // addff in python version
        let mut forcefields = Vec::new();
        let mut offxml = None;
        let mut openff_forcefield = None;
        let mut ffdata_isxml = HashMap::new();
        for fnm in &fnms {
            offxml = Some(fnm.clone());
            let ff = ForceField::load(fnm).unwrap();
            openff_forcefield = Some(ff.clone());
            ffdata_isxml.insert(fnm.clone(), true);
            forcefields.push(ff);
        }

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
            map: todo!(),
            plist: Vec::new(),
            patoms: todo!(),
            pfields: todo!(),
            ptree: Digraph,
            offxml_unit_strs: todo!(),
            rs: todo!(),
            tm: todo!(),
            tmi,
            excision: Vec::new(),
            np,
            pvals0,
            atomnames: todo!(),
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

        const CONCERN: [&str; 3] = ["polarizability", "epsilon", "VDWT"];

        for i in 0..self.np {
            if self.plist.iter().any(|j| CONCERN.contains(&j.as_str()))
                && pvals[i] * self.pvals0[i] < 0.0
            {
                pvals[i] = 0.0;
            }
        }

        // redirect parameters (for the fusion penalty function)
        // TODO another redirect
        // for p in self.redirect

        pvals
    }

    pub(crate) fn make_redirect(&self, mvals: Dvec) {
        todo!()
    }
}
