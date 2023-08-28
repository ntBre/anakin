use std::{collections::HashMap, error::Error, path::Path, sync::LazyLock};

use crate::{
    molecule::pdb::Pdb,
    utils::{np_max, np_min},
    Dmat,
};

use self::{pdb::Record, xyz::Xyz};

mod pdb;
mod xyz;

static ATOMIC_NUMBERS: LazyLock<HashMap<&str, usize>> = LazyLock::new(|| {
    [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
        "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
        "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
        "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
        "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
        "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
        "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
        "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    ]
    .into_iter()
    .enumerate()
    .map(|(a, b)| (b, a + 1))
    .collect()
});

const RADII: [f64; 96] = [
    0.31, 0.28, // H and He
    1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, // First row elements
    0.00, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, // Second row elements
    2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.61, 1.52, 1.50, 1.24, 1.32, 1.22,
    1.22, 1.20, 1.19, 1.20, 1.20,
    1.16, // Third row elements, K through Kr
    2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44,
    1.42, 1.39, 1.39, 1.38, 1.39,
    1.40, // Fourth row elements, Rb through Xe
    2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92,
    1.92, 1.89, 1.90, 1.87, // Fifth row elements, s and f blocks
    1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46,
    1.48, 1.40, 1.50, 1.50, // Fifth row elements, d and p blocks
    2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69,
];

// TODO surely I need the xyzs from the .xyz file too. I expect we'll be reading
// these files over and over in the future

/// Molecule.Data gets copied into from the read_pdb return value (Answer), so
/// our Molecule has the same fields as Answer
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Molecule {
    /// an N x 3 matrix containing the x, y, z coordinates for the atoms in the
    /// molecule. comes from read_pdb, overwriting earlier values from read_xyz
    xyzs: Dmat,

    /// vector of atomic symbols
    elem: Vec<String>,

    /// vector of bond connections. TODO comes from build_topology
    bonds: Vec<(usize, usize)>,
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
        let xyz = Xyz::load(load_fnm)?;
        let fnm = topology.as_ref();
        let pdb = Pdb::load(fnm)?;

        let mut xyzs = Vec::new();
        let mut rows = 0;
        for r in pdb.records {
            if let Record::Atom { x, y, z, .. } = r {
                xyzs.extend([x, y, z]);
                rows += 1;
            };
        }

        let mut ret = Self {
            xyzs: Dmat::from_row_slice(rows, 3, &xyzs),
            elem: xyz.elem,
            bonds: Vec::new(),
        };

        ret.build_bonds();

        Ok(ret)
    }

    pub(crate) fn from_path(fnm: impl AsRef<Path>) -> Self {
        todo!();
    }

    pub(crate) fn len(&self) -> usize {
        // the fields should be parallel arrays, so just pick an easy one
        self.elem.len()
    }

    fn build_bonds(&mut self) {
        // build a vec of covalent radii
        let mut r = Vec::new();
        for elem in &self.elem {
            r.push(if let Some(idx) = ATOMIC_NUMBERS.get(elem.as_str()) {
                RADII[idx - 1]
            } else {
                0.0
            })
        }

        // minimum distance for considering two atoms bonded
        let mindist = 1.0;

        // this looks so stupid, but I'm doing it for now. this is just a big
        // list of pairs
        let na = self.elem.len();
        let mut atom_iterator = Vec::new();
        for i in 0..na {
            for j in i + 1..na {
                atom_iterator.push((i, j));
            }
        }

        // taken from self.top_settings["Fac"] but is this a default?
        let fac = 1.2;

        // oh look now we're splitting up the pairs.......
        let bt0: Vec<_> = atom_iterator.iter().map(|(a, _)| r[*a]).collect();
        let bt1: Vec<_> = atom_iterator.iter().map(|(_, a)| r[*a]).collect();

        // holy moly what is going on here
        let mut bond_thresh = Vec::new();
        for (bt0, bt1) in bt0.iter().zip(&bt1) {
            let bt = (bt0 + bt1) * fac;
            // just hilarious
            let bt = (bt > mindist) as usize as f64 * bt
                + (bt < mindist) as usize as f64 * mindist;
            // TODO
            // these conditions evaluate to 1 if true, so if bt > mindist { bt }
            // else if bt < mindist { mindist }. aka max(bt, mindist) ???????
            bond_thresh.push(bt);
        }

        // if not toppbc
        let dxij = atom_contact(&self.xyzs, &atom_iterator);

        // create a list of atoms that each atom is bonded to
        let mut atom_bonds = vec![vec![]; na];
        for (i, d) in dxij.iter().enumerate() {
            if d < &bond_thresh[i] {
                let (ii, jj) = atom_iterator[i];
                if ii == jj {
                    // omfg we didn't even filter out the same atom pairs
                    // before... we've done all this work just to filter them
                    // here lmao
                    continue;
                }
                // oh boy I can't wait to see what we do with this... it's
                // looking a whole lot like what we already had in
                // atom_iterator!!
                atom_bonds[ii].push(jj);
                atom_bonds[jj].push(ii);
            }
        }

        let mut bond_list = Vec::new();
        for (i, bi) in atom_bonds.iter().enumerate() {
            for &j in bi {
                if i == j {
                    continue;
                }
                // TODO fragment/resid check
                if i < j {
                    bond_list.push((i, j));
                } else {
                    bond_list.push((j, i));
                }
            }
        }

        bond_list.sort();
        // holy moly whyyy
        bond_list.dedup();
        self.bonds = bond_list;
    }
}

/// compute the distances between the atoms in `xyzs` requested in `pairs`
fn atom_contact(xyzs: &Dmat, pairs: &[(usize, usize)]) -> Vec<f64> {
    let mut ret = Vec::with_capacity(pairs.len());
    for (a, b) in pairs {
        ret.push((xyzs.row(*a) - xyzs.row(*b)).norm());
    }
    ret
}

#[cfg(test)]
mod tests {
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
        assert_eq!(got.xyzs, want.xyzs);
        assert_eq!(got.elem, want.elem);
        assert_eq!(got.bonds, want.bonds);
        assert_eq!(got, want);
    }
}
