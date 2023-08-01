use std::path::Path;

pub(crate) struct Molecule;

impl Molecule {
    /// load a [Molecule] from `topology` and frames from `filename`
    pub(crate) fn new<P, Q>(filename: P, topology: Q) -> Self
    where
        P: AsRef<Path>,
        Q: AsRef<Path>,
    {
        todo!()
    }

    pub(crate) fn len(&self) -> usize {
        todo!()
    }
}
