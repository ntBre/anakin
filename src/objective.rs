use std::collections::HashMap;

use crate::{
    config::Config, forcefield::FF, work_queue::WorkQueue, Dmat, Dvec,
};

use self::penalty::{Penalty, PenaltyType};

pub(crate) mod penalty;

pub(crate) struct Target {
    pub(crate) good_step: bool,

    /// current step size for finite differences
    pub(crate) h: f64,

    /// bsave in Python but I think the b means bool..
    save: bool,

    name: String,

    weight: f64,

    evaluated: bool,
}

impl Target {
    fn stage(&self, mvals: Dvec, order_1: bool, order_2: bool) {
        todo!()
    }

    fn get_x(&self, mvals: &Dvec) -> Extra {
        todo!()
    }

    fn get_g(&self, mvals: &Dvec) -> Extra {
        todo!()
    }

    fn get_h(&self, mvals: &Dvec) -> Extra {
        todo!()
    }
}

#[derive(Default)]
// TODO consider renaming this once I figure out what it's supposed to be. only
// one instance of it is called regularization. other times we're associating
// target names with it
pub struct Regularization {
    w: f64,
    x: f64,
}

impl Regularization {
    pub fn new(w: f64, x: f64) -> Self {
        Self { w, x }
    }
}

pub(crate) struct Extra(f64, Dvec, Dmat);

#[derive(Clone)]
pub(crate) struct ObjMap {
    pub(crate) x0: f64,
    pub(crate) g0: Dvec,
    pub(crate) h0: Dmat,

    /// objective function value
    pub(crate) x: f64,

    /// objective function first derivative
    pub(crate) g: Dvec,

    /// objective function second derivative
    pub(crate) h: Dmat,
}

impl ObjMap {
    pub(crate) fn zeros(size: usize) -> Self {
        Self {
            x0: 0.0,
            g0: Dvec::zeros(size),
            h0: Dmat::zeros(size, size),
            x: 0.0,
            g: Dvec::zeros(size),
            h: Dmat::zeros(size, size),
        }
    }
}

// TODO figure out how to obtain this. going to have to keep track of it on one
// of these structs or maybe just set a global static (...) we're definitely not
// going to examine the stack trace of the program like the python
// version..........
fn in_fd() -> bool {
    false
}

pub struct Objective {
    pub(crate) penalty: Penalty,
    normalize_weights: bool,
    forcefield: FF,
    pub(crate) targets: Vec<Target>,

    /// in Python this is an entry in the ObjDict map. TODO consider making
    /// it an Option depending on how it's used
    obj_map: HashMap<String, Regularization>,

    /// assuming this means total weight
    wtot: f64,
}

impl Objective {
    pub fn new(config: &Config, forcefield: FF) -> Self {
        let ptype = match config.penalty_type.to_lowercase().as_str() {
            "hyp" | "hyper" | "l1" | "hyperbola" | "hyperbolic" => {
                PenaltyType::Hyperbolic
            }
            "para" | "parabola" | "l2" | "quadratic" | "parabolic" => {
                PenaltyType::Parabolic
            }
            "box" => PenaltyType::Box,
            _ => {
                panic!("unrecognized penalty_type: {}", config.penalty_type)
            }
        };
        let penalty = Penalty {
            fadd: config.penalty_additive,
            fmul: config.penalty_multiplicative,
            a: config.penalty_alpha,
            b: config.penalty_hyperbolic_b,
            p: config.penalty_power,
            ptype,
        };
        Self {
            forcefield,
            targets: Vec::new(),
            obj_map: HashMap::new(),
            penalty,
            // TODO should be some of weights from targets
            wtot: 1.0,
            normalize_weights: config.normalize_weights,
        }
    }

    pub(crate) fn full(&mut self, vals: Dvec, order: i32) -> ObjMap {
        // TODO borrow here if it doesn't get consumed
        let mut objective = self.target_terms(vals.clone(), order);
        let extra = if self.forcefield.use_pvals {
            self.penalty
                .compute(self.forcefield.create_mvals(vals), &mut objective)
        } else {
            self.penalty.compute(vals, &mut objective)
        };
        objective.x0 = objective.x;
        objective.g0 = objective.g.clone();
        objective.h0 = objective.h.clone();

        if !in_fd() {
            self.obj_map.insert(
                "Regularization".to_owned(),
                Regularization::new(1.0, extra.0),
            );
        }

        objective.x += extra.0;
        objective.g += extra.1;
        objective.h += extra.2;

        objective
    }

    fn target_terms(&mut self, mvals: Dvec, order: i32) -> ObjMap {
        let mut objective = ObjMap::zeros(self.forcefield.np);
        for tgt in &self.targets {
            // TODO consider borrowing
            tgt.stage(mvals.clone(), order >= 1, order >= 2);
        }

        // TODO might be supposed to use an existing work queue. see
        // nifty/getWorkQueue at some point. this might even be a psqs
        // situation.
        let wq = WorkQueue::new();
        for tgt in self.targets.iter_mut() {
            tgt.save = true;
            let ans: Extra = match order {
                0 => tgt.get_x(&mvals),
                1 => tgt.get_g(&mvals),
                2 => tgt.get_h(&mvals),
                _ => unimplemented!(),
            };
            if !in_fd() {
                self.obj_map.insert(
                    tgt.name.clone(),
                    Regularization::new(tgt.weight / self.wtot, ans.0),
                );
            }
            objective.x += ans.0 * tgt.weight / self.wtot;
            objective.g += ans.1 * tgt.weight / self.wtot;
            objective.h += ans.2 * tgt.weight / self.wtot;
        }

        for tgt in self.targets.iter_mut() {
            tgt.evaluated = true;
        }

        // prevent exact zeros on Hessian diagonal
        for i in 0..self.forcefield.np {
            if objective.h[(i, i)] == 0.0 {
                objective.h[(i, i)] = 1.0;
            }
        }

        objective
    }
}
