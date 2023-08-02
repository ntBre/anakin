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
pub(crate) struct Molecule {
    /// an N x 3 matrix containing the x, y, z coordinates for the atoms in the
    /// molecule
    xyzs: Dmat,

    /// vector of alternate locations
    altloc: Vec<String>,

    /// vector of atom names (element + number)
    atomname: Vec<String>,

    /// vector of atom residue ids
    resid: Vec<usize>,

    /// vector of atom residue names (seems to be UNK always for us)
    resname: Vec<String>,

    /// vector of atomic symbols
    elem: Vec<String>,

    /// vector of bond connections
    bonds: Vec<(usize, usize)>,
}

impl Molecule {
    /// load a [Molecule] from `topology` and frames from `filename`. We only
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
                    let a = atoms[0];
                    for b in &atoms[1..] {
                        bonds.push((a.min(*b), a.max(*b)));
                    }
                }
                _ => {
                    continue;
                }
            };
        }

        Ok(Self {
            xyzs: Dmat::from_row_slice(rows, 3, &xyzs),
            altloc,
            atomname,
            resid,
            resname,
            elem,
            bonds,
        })
    }

    pub(crate) fn len(&self) -> usize {
        // the fields should be parallel arrays, so just pick an easy one
        self.elem.len()
    }
}
