#[derive(Clone)]
pub struct FF {
    /// The total number of parameters.
    pub(crate) np: usize,

    /// Bypass the transformation matrix and use the physical parameters
    /// directly used in: Creating the force field; advanced usage, be careful.
    /// TODO taken from the config file
    pub(crate) use_pvals: bool,
}

impl FF {
    pub fn new() -> Self {
        FF {
            np: 0,
            use_pvals: false,
        }
    }

    pub(crate) fn create_mvals(&self, vals: Vec<f64>) -> Vec<f64> {
        todo!()
    }
}

impl Default for FF {
    fn default() -> Self {
        Self::new()
    }
}
