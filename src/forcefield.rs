use openff_toolkit::smirnoff::ForceField;

use crate::{config::Config, Dmat, Dvec};

#[derive(Clone)]
pub struct FF {
    /// The total number of parameters.
    pub(crate) np: usize,

    /// Bypass the transformation matrix and use the physical parameters
    /// directly used in: Creating the force field; advanced usage, be careful.
    /// TODO taken from the config file
    pub(crate) use_pvals: bool,

    /// transpose of the transformation matrix. initially empty but used later.
    /// TODO consider option
    tmi: Dmat,

    pvals0: Dvec,

    plist: Vec<String>,

    /// indices to be excluded from the Hessian update
    pub(crate) excision: Vec<usize>,

    /// force field names
    fnms: Vec<String>,

    /// force field directory relative to the project directory
    ffdir: String,

    /// the full force fields
    forcefields: Vec<ForceField>,
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
        let mut forcefields = Vec::new();
        for fnm in &fnms {
            forcefields.push(ForceField::load(fnm).unwrap());
        }

        // TODO there is supposed to be more processing of the force fields here
        // to extract the useful parts. really we probably don't even need the
        // original ForceField structs
        FF {
            np,
            use_pvals: config.use_pvals,
            tmi,
            pvals0,
            plist: Vec::new(),
            excision: Vec::new(),
            fnms,
            ffdir: config.ffdir.clone(),
            forcefields,
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

impl Default for FF {
    // TODO I'm not really sure this makes sense to have. see where it's being
    // used and if we can/need to use something else
    fn default() -> Self {
        Self {
            np: 0,
            use_pvals: false,
            tmi: Dmat::zeros(0, 0),
            pvals0: Dvec::zeros(0),
            plist: Vec::new(),
            excision: Vec::new(),
            fnms: Vec::new(),
            ffdir: String::from("forcefield"),
            forcefields: Vec::new(),
        }
    }
}
