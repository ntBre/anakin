use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use log::warn;

#[derive(Debug)]
pub(crate) enum Record {
    Remark {
        number: usize,
        contents: String,
    },
    Atom {
        serial: usize,
        name: String,
        alt_loc: String,
        res_name: String,
        // chain_id: char, // skipped
        res_seq: usize,
        // i_code: ?? looks like ours are missing these
        x: f64,
        y: f64,
        z: f64,
        occupancy: f64,
        temp_factor: f64,
        element: String,
        // charge: isize, // but we're also missing this
    },
    Ter {
        number: usize,
        res_name: String,
        chain_id: char,
        res_seq: usize,
    },
    Conect {
        atoms: Vec<usize>,
    },
}

/// we'll start with a very literal interpretation of the [PDB spec], limited to
/// the fields we actually need to parse for ForceBalance
///
/// [pdb spec]:
/// (https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)
#[derive(Debug)]
pub(crate) struct Pdb {
    pub(crate) records: Vec<Record>,
}

impl Pdb {
    pub(crate) fn load(
        filename: impl AsRef<Path>,
    ) -> Result<Self, Box<dyn Error>> {
        let f = File::open(filename)?;
        let r = BufReader::new(f);
        let mut records = Vec::new();
        for line in r.lines().flatten() {
            let split: Vec<&str> = line.split_ascii_whitespace().collect();
            if split.is_empty() {
                continue;
            }
            let record = match split[0] {
                "REMARK" => Record::Remark {
                    number: split[1].parse()?,
                    contents: split[2..].join(" "),
                },
                "ATOM" | "HETATM" => Record::Atom {
                    serial: split[1].parse()?,
                    name: split[2].to_owned(),
                    alt_loc: split[3].to_owned(),
                    res_name: split[4].to_owned(),
                    res_seq: split[5].parse()?,
                    x: split[6].parse()?,
                    y: split[7].parse()?,
                    z: split[8].parse()?,
                    occupancy: split[9].parse()?,
                    temp_factor: split[10].parse()?,
                    element: split[11].to_owned(),
                },
                "TER" => Record::Ter {
                    number: split[1].parse()?,
                    res_name: split[2].to_owned(),
                    chain_id: split[3].parse()?,
                    res_seq: split[4].parse()?,
                },
                "CONECT" => {
                    // not iterator because I want ? to return
                    let mut atoms = Vec::new();
                    for s in &split[1..] {
                        atoms.push(s.parse()?);
                    }
                    Record::Conect { atoms }
                }
                f => {
                    warn!("ignoring unknown field {f}");
                    continue;
                }
            };
            records.push(record);
        }

        Ok(Pdb { records })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_pdb() {
        let got = Pdb::load("testfiles/test.pdb").unwrap();
        assert_eq!(got.records.len(), 54);
    }
}
