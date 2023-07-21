use crate::{
    forcefield::FF,
    objective::{ObjMap, Objective},
};

/// a struct representing all possible convergence criteria
struct Criteria {
    step: bool,
    grad: bool,
    obj: bool,
}

impl Criteria {
    fn new() -> Self {
        Self {
            step: false,
            grad: false,
            obj: false,
        }
    }
}

pub struct Optimizer {
    objective: Objective,
    forcefield: FF,

    trust0: f64,
    adapt_fac: f64,
    mvals0: Vec<f64>,

    /// initial step size for numerical finite difference
    h0: f64,

    /// current step size for numerical finite difference
    h: f64,

    /// finite difference factor
    fdf: f64,

    /// Some target types (looks like liquid, lipid, and thermo) introduce
    /// uncertainty into the objective function, necessitating re-evaluation of
    /// the objective function when a step is rejected.
    uncert: bool,

    convergence_objective: f64,

    iteration: usize,

    good_step: bool,
}

impl Optimizer {
    pub fn new(objective: Objective, forcefield: FF) -> Self {
        Self {
            objective,
            forcefield,
            /// TODO take these from config file. I'm taking default values from
            /// the config file produced by the example script. maybe these
            /// should be on Optimizer::default
            trust0: 0.1,
            adapt_fac: 0.25,
            uncert: false,
            convergence_objective: 1e-4,
            iteration: 0,
            good_step: false,
            mvals0: Vec::new(), // TODO option? or empty vector enough
            h0: 0.001,          // taken from finite_difference_h config
            h: 0.001,
            fdf: 0.1,
        }
    }

    pub fn run(&mut self) {
        // TODO can call the appropriate optimizer here. in my case it's
        // OPTIMIZE, which calls NewtonRaphson
        self.newton_raphson();
    }

    fn newton_raphson(&mut self) {
        self.main_optimizer(false);
    }

    /// The main ForceBalance adaptive trust-radius pseudo-Newton optimizer.
    /// Tried and true in many situations. :)
    ///
    /// Usually this function is called with the BFGS or NewtonRaphson
    /// method. The NewtonRaphson method is consistently the best method I
    /// have, because I always provide at least an approximate Hessian to
    /// the objective function. The BFGS method works well, but if gradients
    /// are cheap the SciPy_BFGS method also works nicely.
    ///
    /// The method adaptively changes the step size. If the step is
    /// sufficiently good (i.e. the objective function goes down by a large
    /// fraction of the predicted decrease), then the step size is
    /// increased; if the step is bad, then it rejects the step and tries
    /// again.
    ///
    /// The optimization is terminated after either a function value or step
    /// size tolerance is reached.
    ///
    /// `bfgs` switches between BFGS (`true`) or Newton-Raphson (`false`)
    fn main_optimizer(&mut self, bfgs: bool) -> ! {
        // this is only for printing
        let detail = if self.trust0 < 0.0 {
            "hessian diagonal search"
        } else if self.adapt_fac != 0.0 {
            "adaptive trust radius"
        } else {
            "trust radius"
        };
        eprintln!("{detail}");

        // warn if optimization is unlikely to converge
        if self.uncert && self.convergence_objective < 1e-3 {
            eprintln!(
                "Condensed phase targets detected - may not converge \
			   with current choice of convergence objective {:e}. \
			   Recommended range is 1e-2 - 1e-1 for this option.",
                self.convergence_objective
            )
        }

        // initialize variables

        // order of derivatives
        let ord = if bfgs { 1 } else { 2 };
        self.set_good_step(true);
        let best_step = 1; // TODO should be zero?
                           // objective function history
        let x_hist: Vec<usize> = Vec::new(); // TODO dummy type
                                             // trust radius
        let trust = self.trust0.abs();
        // current value of the parameters
        let xk = self.mvals0.clone();
        // current optimization step
        let dx = vec![0.0; self.forcefield.np];
        // length of the current optimization step
        let ndx = 0.0;
        // ratio of actual objective function change to expected change
        let quality = 0.0;
        // threshold for "low quality step" which decreases trust radius
        let threlq = 0.25;
        // threshold for "high quality step" which increases trust radius
        let threhq = 0.75;
        // optimization steps before this one are ineligible for
        // consideration for "best step"
        let best_start = 0;
        // convergence criteria
        let criteria_satisfied = Criteria::new();

        // non-linear iterations
        loop {
            let using_checkpoint = false;
            let (x, g, h) = if using_checkpoint {
                todo!();
            } else {
                self.adjh(trust);
                // TODO see if xk can be borrowed here. it might even need to be
                // mutated?
                let ObjMap { x, g, h, .. } =
                    self.objective.full(xk.clone(), ord);
                (x, g, h)
            };
        }
    }

    pub fn set_good_step(&mut self, good_step: bool) {
        self.good_step = good_step;
        for target in self.objective.targets.iter_mut() {
            target.good_step = good_step;
        }
    }

    fn adjh(&mut self, trust: f64) {
        // the finite difference step size should be at most 1% of the trust
        // radius
        let h = f64::min(self.fdf * trust, self.h0);
        if h != self.h {
            // setting finite difference step size
            self.h = h;
            for target in self.objective.targets.iter_mut() {
                target.h = h;
            }
        }
    }
}
