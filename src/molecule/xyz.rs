use std::{
    error::Error,
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use crate::Dmat;

#[derive(Default)]
pub(crate) struct Xyz {
    pub(crate) elem: Vec<String>,
    xyzs: Vec<Dmat>,
    comms: Vec<String>,
}

impl Xyz {
    /// Load an XYZ file of the form:
    ///
    /// ```text
    /// 3
    /// torsion grid (-165,)
    /// C        1.6989050157    0.0167887368   -3.5226452110
    /// C        1.2823068536    0.1330978850   -2.0673677464
    /// O        2.1832213014   -0.6763254205   -1.2942688917
    /// 3
    /// torsion grid (-150,)
    /// C        1.7069317443    0.0328736165   -3.5302888699
    /// C        1.2838266044    0.1318941209   -2.0755434894
    /// O        2.2016267160   -0.6592414574   -1.3041288483
    /// ```
    ///
    /// In other words, zero or more records consisting of a number of atoms, a
    /// comment line, and a Cartesian geometry.
    pub(crate) fn load(fnm: impl AsRef<Path>) -> Result<Self, Box<dyn Error>> {
        let f = File::open(fnm)?;
        let r = BufReader::new(f);
        let mut ret = Self::default();
        let mut na = 0;
        let mut blanks = 0;
        let mut comment = false;
        let mut xyz = Vec::new(); // buffer for holding current xyz
        for (ln, line) in r.lines().flatten().enumerate() {
            if comment {
                comment = false;
                ret.comms.push(line);
                continue;
            }
            let split: Vec<_> = line.split_ascii_whitespace().collect();
            if split.is_empty() {
                blanks += 1;
                continue;
            }
            // first line or every n+2 lines - blanks. +2 for the natom and
            // comment lines themselves. I'm not sure XYZ files are supposed to
            // have blanks at all, but this is extra safe
            if ln == 0 {
                na = split[0].parse()?;
                comment = true;
            } else if (ln - blanks) % (na + 2) == 0 {
                na = split[0].parse()?;
                comment = true;
                ret.xyzs.push(Dmat::from_row_slice(na, 3, &xyz));
                xyz.clear();
                // take the atom labels from the last entry. obviously we're
                // assuming they're all the same. lee-ping takes them from the
                // first entry "if len(elem) < na"
                ret.elem.clear();
            } else {
                assert_eq!(
                    split.len(),
                    4,
                    "expected 4 entries on XYZ coordinate line, but got {}. \
                     Line {ln}:\n`{line}`",
                    split.len()
                );
                ret.elem.push(split[0].to_owned());
                // again, not an iterator extend because of ?
                for s in &split[1..] {
                    xyz.push(s.parse()?);
                }
            }
        }
        ret.xyzs.push(Dmat::from_row_slice(na, 3, &xyz));

        Ok(ret)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_xyz() {
        let got = Xyz::load("testfiles/test.xyz").unwrap();
        assert_eq!(got.xyzs.len(), 24);
        assert_eq!(got.comms.len(), 24);
        assert_eq!(got.elem.len(), 26);
    }
}
