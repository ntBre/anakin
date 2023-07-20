use crate::forcefield::FF;

pub(crate) struct Target {
    pub(crate) good_step: bool,

    /// current step size for finite differences
    pub(crate) h: f64,
}

#[derive(Default)]
pub struct Regularization {
    w: f64,
    x: f64,
}

impl Regularization {
    pub fn new(w: f64, x: f64) -> Self {
        Self { w, x }
    }
}

struct Extra(f64, Vec<f64>, Vec<f64>);

struct Penalty;

impl Penalty {
    fn compute(&self, vals: Vec<f64>, objective: &mut ObjMap) -> Extra {
        todo!()
    }
}

pub struct Objective {
    forcefield: FF,
    pub(crate) targets: Vec<Target>,

    /// in Python this is an entry in the ObjDict map. TODO consider making
    /// it an Option depending on how it's used
    regularization: Regularization,

    penalty: Penalty,
}

pub(crate) struct ObjMap {
    pub(crate) x0: f64,
    pub(crate) g0: Vec<f64>,
    pub(crate) h0: Vec<f64>,
    pub(crate) x: f64,
    pub(crate) g: Vec<f64>,
    pub(crate) h: Vec<f64>,
}

impl Objective {
    pub fn new(forcefield: FF) -> Self {
        Self {
            forcefield,
            targets: Vec::new(),
            regularization: Default::default(),
            penalty: Penalty,
        }
    }

    pub(crate) fn full(&mut self, vals: Vec<f64>, order: i32) -> ObjMap {
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

        // TODO figure out how to obtain this. going to have to keep track
        // of it on one of these structs or maybe just set a global static
        // (...) we're definitely not going to examine the stack trace of
        // the program like the python version..........
        let in_finite_difference = false;
        if !in_finite_difference {
            self.regularization = Regularization::new(1.0, extra.0);
        }

        objective.x += extra.0;
        // TODO just use + if I switch to some kind of array package
        for (i, g) in objective.g.iter_mut().enumerate() {
            *g += extra.1[i];
        }
        for (i, h) in objective.h.iter_mut().enumerate() {
            *h += extra.2[i];
        }

        objective
    }

    fn target_terms(&self, vals: Vec<f64>, order: i32) -> ObjMap {
        todo!()
    }
}
