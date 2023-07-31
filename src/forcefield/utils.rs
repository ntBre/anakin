use std::collections::HashMap;

use openff_toolkit::smirnoff::ForceField;

use crate::{config::Config, forcefield::Param, Dvec};

pub(super) fn rsmake(
    bonds_to_optimize: &[Param],
    openff_forcefield: &Option<ForceField>,
    angles_to_optimize: &[Param],
    propers_to_optimize: &[Param],
    config: &Config,
    pvals0: &Vec<f64>,
) -> Dvec {
    // extracting the maximum values of the parameters for computing the
    // rescaling factors.
    let mut max_bond = HashMap::new();
    for (i, param) in bonds_to_optimize.iter().enumerate() {
        if let Param::Opt { inner } = param {
            for (_, typ, _unit) in inner {
                let val = openff_forcefield.as_ref().unwrap().bonds[i]
                    .as_hash(typ)
                    .unwrap()
                    .value
                    .abs();
                let cur = max_bond.entry(typ.clone()).or_insert(val);
                if val > *cur {
                    *cur = val;
                }
            }
        }
    }

    let mut max_angle = HashMap::new();
    for (i, param) in angles_to_optimize.iter().enumerate() {
        if let Param::Opt { inner } = param {
            for (_, typ, _unit) in inner {
                let val = openff_forcefield.as_ref().unwrap().angles[i]
                    .as_hash(typ)
                    .unwrap()
                    .value
                    .abs();
                let cur = max_angle.entry(typ.clone()).or_insert(val);
                if val > *cur {
                    *cur = val;
                }
            }
        }
    }

    let mut max_proper = HashMap::new();
    for (i, param) in propers_to_optimize.iter().enumerate() {
        if let Param::Opt { inner } = param {
            for (_, typ, _unit) in inner {
                let val = openff_forcefield.as_ref().unwrap().proper_torsions
                    [i]
                    .as_hash(typ)
                    .unwrap()
                    .value
                    .abs();
                let cur = max_proper.entry(typ.clone()).or_insert(val);
                if val > *cur {
                    *cur = val;
                }
            }
        }
    }

    // overwrite from priors
    if let Some(map) = &config.priors {
        if let Some(bonds) = map.get("bonds") {
            for (key, value) in bonds {
                max_bond.insert(key.to_string(), *value);
            }
        }
        if let Some(angles) = map.get("angles") {
            for (key, value) in angles {
                max_angle.insert(key.to_string(), *value);
            }
        }
        if let Some(propers) = map.get("propertorsions") {
            for (key, value) in propers {
                if key == "k" {
                    for key in ["k1", "k2", "k3", "k4", "k5", "k6"] {
                        max_proper.insert(key.to_string(), *value);
                    }
                }
            }
        }
    }

    let mut rs = Dvec::from_element(pvals0.len(), 1.0);
    'outer: for i in 0..pvals0.len() {
        for bond in bonds_to_optimize {
            if let Param::Opt { inner } = bond {
                for (x, typ, _) in inner {
                    if *x == i {
                        rs[i] = *max_bond.get(typ).unwrap();
                        continue 'outer;
                    }
                }
            }
        }
        for angle in angles_to_optimize {
            if let Param::Opt { inner } = angle {
                for (x, typ, _) in inner {
                    if *x == i {
                        rs[i] = *max_angle.get(typ).unwrap();
                        continue 'outer;
                    }
                }
            }
        }
        for proper in propers_to_optimize {
            if let Param::Opt { inner } = proper {
                for (x, typ, _) in inner {
                    if *x == i {
                        rs[i] = *max_proper.get(typ).unwrap();
                        continue 'outer;
                    }
                }
            }
        }
    }
    rs
}

// TODO might be able to factor these out with the use of parameter_handlers

pub(super) fn add_bonds(
    ff: &ForceField,
    pvals0: &mut Vec<f64>,
    bonds_to_optimize: &mut Vec<Param>,
) {
    for (i, bond) in (&ff.bonds).into_iter().enumerate() {
        if let Some(to_optimize) = &bond.parameterize {
            let mut inner = Vec::new();
            for param in to_optimize.split(',').map(str::trim) {
                let p = bond.as_hash(param).unwrap();
                inner.push((pvals0.len(), param.to_owned(), p.unit.clone()));
                pvals0.push(p.value);
            }
            bonds_to_optimize.push(Param::Opt { inner });
        } else {
            bonds_to_optimize.push(Param::NoOpt);
        }
    }
}

pub(super) fn add_angles(
    ff: &ForceField,
    pvals0: &mut Vec<f64>,
    angles_to_optimize: &mut Vec<Param>,
) {
    for (i, angle) in (&ff.angles).into_iter().enumerate() {
        if let Some(to_optimize) = &angle.parameterize {
            let mut inner = Vec::new();
            for param in to_optimize.split(',').map(str::trim) {
                let p = angle.as_hash(param).unwrap();
                inner.push((pvals0.len(), param.to_owned(), p.unit.clone()));
                pvals0.push(p.value);
            }
            angles_to_optimize.push(Param::Opt { inner });
        } else {
            angles_to_optimize.push(Param::NoOpt);
        }
    }
}

pub(super) fn add_propers(
    ff: &ForceField,
    pvals0: &mut Vec<f64>,
    propers_to_optimize: &mut Vec<Param>,
) {
    for (i, proper) in (&ff.proper_torsions).into_iter().enumerate() {
        if let Some(to_optimize) = &proper.parameterize {
            let mut inner = Vec::new();
            for param in to_optimize.split(',').map(str::trim) {
                let p = proper.as_hash(param).unwrap();
                inner.push((pvals0.len(), param.to_owned(), p.unit.clone()));
                pvals0.push(p.value);
            }
            propers_to_optimize.push(Param::Opt { inner });
        } else {
            propers_to_optimize.push(Param::NoOpt);
        }
    }
}
