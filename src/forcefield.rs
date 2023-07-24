use crate::{Dmat, Dvec};

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
}

impl FF {
    pub fn new() -> Self {
        FF {
            np: 0,
            use_pvals: false,
            tmi: Dmat::zeros(0, 0),
            pvals0: Dvec::zeros(0),
            plist: Vec::new(),
            excision: Vec::new(),
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
    fn default() -> Self {
        Self::new()
    }
}
