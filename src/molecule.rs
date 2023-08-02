use std::{error::Error, path::Path};

use crate::{molecule::pdb::Pdb, Dmat};

use self::pdb::Record;

mod pdb;

// NOTE: skipped chain_id
//
// NOTE: skipping icode because I skipped it in the PDB parser too
//
// NOTE: also skipping `comms`, why would I need the comments lol
//
// NOTE: skipping `terminal` which appears to be a boolean vector saying whether
// or not an atom is the terminal one... I'm operating under the assumption that
// there is one molecule per PDB file right now, so in addition to being a
// terrible representation, it's also unnecessary for me.

/// Molecule.Data gets copied into from the read_pdb return value (Answer), so
/// our Molecule has the same fields as Answer
#[derive(Debug, PartialEq)]
pub(crate) struct Molecule {
    /// an N x 3 matrix containing the x, y, z coordinates for the atoms in the
    /// molecule
    xyzs: Dmat,

    /// vector of atomic symbols
    elem: Vec<String>,

    /// vector of bond connections
    bonds: Vec<(usize, usize)>,
    // TODO augment with fields from read_xyz
}

impl Molecule {
    /// load a [Molecule] from `filename` and frames from `topology`. We only
    /// know how to load frames from an xyz and the topology from a pdb
    pub(crate) fn new<P, Q>(
        filename: P,
        topology: Q,
    ) -> Result<Self, Box<dyn Error>>
    where
        P: AsRef<Path>,
        Q: AsRef<Path>,
    {
        let load_fnm = filename.as_ref();
        let fnm = topology.as_ref();
        let pdb = Pdb::load(fnm)?;

        let mut xyzs = Vec::new();
        let mut altloc = Vec::new();
        let mut atomname = Vec::new();
        let mut resid = Vec::new();
        let mut resname = Vec::new();
        let mut elem = Vec::new();
        let mut bonds = Vec::new();
        let mut rows = 0;
        for r in pdb.records {
            match r {
                Record::Atom {
                    name,
                    alt_loc,
                    res_name,
                    res_seq,
                    x,
                    y,
                    z,
                    element,
                    ..
                } => {
                    altloc.push(alt_loc);
                    atomname.push(name);
                    resid.push(res_seq);
                    resname.push(res_name);
                    elem.push(element);
                    xyzs.extend([x, y, z]);
                    rows += 1;
                }
                Record::Conect { atoms } => {
                    let a = atoms[0] - 1;
                    for b in &atoms[1..] {
                        let b = b - 1;
                        bonds.push((a.min(b), a.max(b)));
                    }
                }
                _ => {
                    continue;
                }
            };
        }

        Ok(Self {
            xyzs: Dmat::from_row_slice(rows, 3, &xyzs),
            elem,
            bonds,
        })
    }

    pub(crate) fn len(&self) -> usize {
        // the fields should be parallel arrays, so just pick an easy one
        self.elem.len()
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use nalgebra::dmatrix;

    use super::*;

    macro_rules! string {
        ($($s:expr$(,)*)*) => {
            vec![$($s.to_owned(),)*]
        }
    }

    #[test]
    fn load_molecule() {
        let got =
            Molecule::new("testfiles/test.xyz", "testfiles/test.pdb").unwrap();
        let want = Molecule {
            xyzs: dmatrix! [
                1.699,  0.017, -3.523;
                1.282,  0.133, -2.067;
                2.183, -0.676, -1.294;
                2.019, -0.731,  0.054;
                2.952, -1.526,  0.748;
                2.875, -1.651,  2.129;
                1.859, -0.979,  2.825;
                0.935, -0.181,  2.15;
                1.009, -0.056,  0.761;
                1.73 , -1.167,  4.598;
                1.047,  0.006,  5.128;
                3.017, -1.64 ,  5.116;
                0.701, -2.526,  4.818;
                -0.481, -2.457,  4.021;
                1.034,  0.621, -4.148;
                2.724,  0.373, -3.661;
                1.643, -1.022, -3.859;
                0.254, -0.224, -1.923;
                1.34 ,  1.174, -1.723;
                3.734, -2.023,  0.181;
                3.605, -2.241,  2.674;
                0.166,  0.339,  2.712;
                0.288,  0.569,  0.249;
                1.284, -3.352,  4.644;
                -1.128, -3.149,  4.395;
                -0.277, -2.728,  3.055;
            ],
            elem: string![
                "C", "C", "O", "C", "C", "C", "C", "C", "C", "S", "O", "O",
                "N", "N", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H",
                "H", "H",
            ],
            #[rustfmt::skip]
            bonds: vec![
                (0, 1), (0, 14), (0, 15), (0, 16), (1, 2), (1, 17), (1, 18),
                (2, 3), (3, 4), (3, 8), (4, 5), (4, 19), (5, 6), (5, 20),
                (6, 7), (6, 9), (7, 8), (7, 21), (8, 22), (9, 10), (9, 11),
                (9, 12), (12, 13), (12, 23), (13, 24), (13, 25),
            ],
        };
        assert_eq!(got, want);
    }
}
