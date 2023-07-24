use std::cell::RefCell;

use nalgebra::SVD;

use crate::{
    forcefield::FF,
    objective::{penalty::PenaltyType, ObjMap, Objective},
    optimizer::solvers::{hyper_solver, para_solver},
    utils::{invert_svd, std_dev},
    Dmat, Dvec,
};

mod optimize;
mod solvers;

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

    fn sum(&self) -> usize {
        self.step as usize + self.grad as usize + self.obj as usize
    }
}

pub struct Optimizer {
    objective: RefCell<Objective>,
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

    /// (float) Convergence criterion of gradient norm used in: Main Optimizer
    convergence_gradient: f64,

    /// (float) Convergence criterion of step size (just needs to fall below
    /// this threshold) used in: Main Optimizer
    convergence_step: f64,

    /// (int) The number of convergence criteria that must be met for main
    /// optimizer to converge used in: Main Optimizer
    criteria: usize,

    iteration: usize,

    /// initial iteration number, used for checkpoints
    iter_init: usize,

    /// (int) Maximum number of steps in an optimization used in: Main Optimizer
    good_step: bool,

    /// (float) Error tolerance; the optimizer will only reject steps that
    /// increase the objective function by more than this number. used in: Main
    /// Optimizer
    err_tol: f64,

    /// (float) Minimum trust radius (if the trust radius is tiny, then noisy
    /// optimizations become really gnarly) used in: Main Optimizer
    mintrust: f64,

    /// (int) Maximum number of steps in an optimization used in: Main Optimizer
    maxstep: usize,

    /// (float) Damping factor that ties down the trust radius to trust0;
    /// decrease for a more variable step size. used in: Main Optimizer
    adapt_damp: f64,

    hist: usize,

    /// (bool) Allow convergence on "low quality" steps
    converge_lowq: bool,

    /// (float) Optimization will "fail" if step falls below this size used in:
    /// Main Optimizer
    step_lowerbound: f64,

    eps: f64,

    /// levenberg-marquardt guess. lm_guess from config
    lmg: f64,

    /// (float) Search tolerance; used only when trust radius is negative,
    /// dictates convergence threshold of nonlinear search. search_tolerance
    /// from config
    search_tol: f64,
}

impl Optimizer {
    pub fn new(objective: Objective, forcefield: FF) -> Self {
        Self {
            objective: RefCell::new(objective),
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
            iter_init: 0,
            err_tol: 0.0,
            mintrust: 0.0,
            maxstep: 100,
            adapt_damp: 0.5,
            hist: 0,
            converge_lowq: false,
            convergence_gradient: 1e-3,
            convergence_step: 1e-4,
            criteria: 1,
            // eig_lowerbound. this comment got moved around. is this the right
            // place?
            step_lowerbound: 1e-6,
            eps: 1e-4,
            lmg: 1.0,
            search_tol: 1e-4,
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
    fn main_optimizer(&mut self, bfgs: bool) -> Dvec {
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
        let mut x_hist: Vec<f64> = Vec::new();
        // trust radius
        let mut trust = self.trust0.abs();
        // current value of the parameters
        let mut xk = Dvec::from(self.mvals0.clone());
        // current optimization step
        let mut dx = Dvec::zeros(self.forcefield.np);
        // length of the current optimization step
        let mut ndx = 0.0;
        // ratio of actual objective function change to expected change
        let quality = 0.0;
        // threshold for "low quality step" which decreases trust radius
        let threlq = 0.25;
        // threshold for "high quality step" which increases trust radius
        let threhq = 0.75;
        // optimization steps before this one are ineligible for
        // consideration for "best step"
        let mut best_start = 0;
        // convergence criteria
        let mut criteria_satisfied = Criteria::new();

        // non-linear iterations
        let mut x_prev = 0.0;
        let mut xk_prev = Dvec::zeros(0);
        let mut pk_prev = Dvec::zeros(0);
        let mut nxk = 0.0;
        let mut ngd = 0.0;
        let mut g_prev = Dvec::zeros(0);
        let mut h_stor = Dmat::zeros(0, 0);
        let mut dx_expect = 0.0;
        let mut datastore = ObjMap::zeros(0);
        let mut bump = false;
        let mut curr_trust = 0.0;
        let mut x = 0.0;
        let mut g = Dvec::zeros(0);
        let mut h = Dmat::zeros(0, 0);
        loop {
            let using_checkpoint = false;
            let mut data = if using_checkpoint {
                todo!();
            } else {
                self.adjh(trust);
                // TODO see if xk can be borrowed here. it might even need to be
                // mutated?
                self.objective.borrow_mut().full(xk.clone(), ord)
            };

            ObjMap { x, g, h, .. } = data.clone();

            // assess optimization step
            if self.iteration > self.iter_init {
                let dx_actual = x - x_prev;
                // true if x < best step in history
                let mut best_step = x < *x_hist[best_start..]
                    .iter()
                    .min_by(|a, b| a.total_cmp(b))
                    .unwrap_or(&x_hist[best_start]);
                let mut quality = if dx_expect == 0.0 {
                    0.0
                } else {
                    dx_actual / dx_expect
                };
                if x > (x_prev + self.err_tol.max(self.convergence_objective)) {
                    // reject step if objective function rises
                    self.set_good_step(false);
                    xk = xk_prev.clone();
                    let trust =
                        self.mintrust.max(ndx * (1.0 / (1.0 + self.adapt_fac)));
                    if self.uncert {
                        // re-evaluate the objective function and gradient at
                        // the previous parameters
                        self.iteration += 1;
                        best_start = self.iteration - self.iter_init;
                        best_step = true;
                        self.adjh(trust);
                        x_hist.push(x);

                        // TODO write a checkpoint

                        if self.iteration == self.maxstep {
                            eprintln!(
                                "Maximum number of optimization steps \
				       reached ({})",
                                self.iteration
                            );
                            break;
                        }

                        let ObjMap { x, g, h, .. } =
                            self.objective.borrow_mut().full(xk.clone(), ord);
                        self.set_good_step(true);
                        ndx = 0.0;
                        nxk = xk.norm();
                        ngd = g.norm();
                        quality = 0.0;
                    } else {
                        // go back to the start of the loop and take a reduced
                        // step
                        x = x_prev;
                        g = g_prev.clone();
                        h = h_stor.clone();
                        data = datastore.clone();
                    }
                } else {
                    self.set_good_step(true);
                    // adjust step size based on step quality
                    if quality <= threlq && self.trust0 > 0.0 {
                        trust = self
                            .mintrust
                            .max(ndx * (1.0 / (1.0 + self.adapt_fac)));
                    } else if quality >= threhq && bump && self.trust0 > 0.0 {
                        curr_trust = trust;
                        trust += (self.adapt_fac
                            * trust
                            * f64::exp(
                                -1.0 * self.adapt_damp
                                    * (trust / self.trust0 - 1.0),
                            ))
                    }

                    if best_step {
                        // TODO backups
                        // if self.backup

                        // TODO output
                        // self.ff.make(xk);
                    }

                    // Hessian update for BFGS
                    if bfgs {
                        let mut hnew = h_stor.clone();
                        dx = &xk - xk_prev;
                        let dy = &g - g_prev;
                        let mat1 =
                            &dy * &dy.transpose() / (dy.transpose().dot(&dx));
                        let mat2 = ((&hnew * &dx) * &(&hnew * &dx).transpose())
                            / (dx.transpose() * &hnew * dx)[(0, 0)];
                        hnew += mat1 - mat2;
                        h = hnew.clone();
                        data.h = h.clone();
                    }
                    // TODO possibly experimental deleting lines in parameter
                    // file
                }
            }
            x_hist.push(x);

            // Take the stdev over the previous (hist) values. Multiply by 2, so
            // when hist=2 this is simply the difference.
            let l = x_hist.len();
            let stdfront = 2.0
                * if x_hist.len() > self.hist {
                    std_dev(&x_hist[l - self.hist..])
                } else {
                    std_dev(&x_hist)
                };

            let dx2_sign = if x_hist.len() > 1 && x_hist[l - 1] < x_hist[l - 2]
            {
                -1
            } else {
                1
            };

            nxk = xk.norm();
            ngd = g.norm();

            // TODO some printing and bookkeeping on objdict_last. personally I
            // would rather handle all of the logging by dumping to a log file
            // and only printing the best step, or even requiring the user to
            // determine the best step like I did in semp. it adds a lot of
            // bookkeeping to maintain all of these separate best_* things.
            // could also pack them all into a struct at least

            // check convergence criteria
            //
            // Three conditions:
            // 1. If any of the targets contain statistical uncertainty
            //
            // 2. If the step quality is above the low quality threshold or if
            //    we allow convergence at low quality steps
            //
            // 3. If we continued from a previous run and this is the first
            //    objective function evaluation
            if self.uncert
                || (self.converge_lowq || quality > threlq)
                || (self.iteration == self.iter_init && self.iter_init > 0)
            {
                criteria_satisfied.grad = ngd < self.convergence_gradient;

                // TODO comparison assign
                criteria_satisfied.step = ndx < self.convergence_step
                    && ndx >= 0.0
                    && self.iteration > self.iter_init;

                criteria_satisfied.obj = stdfront < self.convergence_objective
                    && x_hist.len() > self.hist;
            }
            if criteria_satisfied.sum() >= self.criteria {
                break;
            }

            // update optimization variables before the next step

            // previous data from objective function call
            datastore = data.clone();
            // previous objective function and derivatives
            x_prev = x;
            g_prev = g.clone();
            h_stor = h.clone();
            // previous optimization variables
            xk_prev = xk.clone();
            // previous physical parameters
            pk_prev = self.forcefield.create_pvals(xk.clone());

            // calculate optimization step
            (dx, dx_expect, bump) = self.step(xk.clone(), data, trust);
            xk += &dx;
            ndx = dx.norm();

            // second convergence check
            if self.uncert
                || (self.converge_lowq || quality > threlq)
                || (self.iteration == self.iter_init && self.iter_init > 0)
            {
                criteria_satisfied.step = ndx < self.convergence_step
                    && ndx >= 0.0
                    && self.iteration > self.iter_init;
            }
            if criteria_satisfied.sum() >= self.criteria {
                break;
            }

            // TODO logging

            // either of these conditions could cause a convergence failure:
            // the maximum number of cycles is reached
            if self.iteration == self.maxstep {
                break;
            }

            // the step size is too small to continue. sometimes happens for
            // very rough objective functions
            if self.iteration > self.iter_init && ndx < self.step_lowerbound {
                break;
            }

            self.iteration += 1;
            if self.trust0 < 0.0 {
                trust = ndx;
            }

            // TODO self.print_vals and print them
            let print_vals = false;
            if print_vals {
                todo!();
            }

            // TODO write checkpoint
        }

        let cnvgd = criteria_satisfied.sum() >= self.criteria;
        if cnvgd {
            println!("optimization converged");
        } else {
            println!("convergence failure");
        }

        xk
    }

    pub fn set_good_step(&mut self, good_step: bool) {
        self.good_step = good_step;
        for target in self.objective.borrow_mut().targets.iter_mut() {
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
            for target in self.objective.borrow_mut().targets.iter_mut() {
                target.h = h;
            }
        }
    }

    fn step(
        &mut self,
        xk: Dvec,
        data: ObjMap,
        trust: f64,
    ) -> (Dvec, f64, bool) {
        // TODO could be a field on `self` but doesnt seem necessary
        let bhyp = !matches!(
            self.objective.borrow().penalty.ptype,
            PenaltyType::Parabolic | PenaltyType::Box
        );

        let (x, mut g, mut h) = if bhyp {
            (data.x0, data.g0, data.h0)
        } else {
            (data.x, data.g, data.h)
        };

        let mut h1 = h.clone();
        for x in &self.forcefield.excision {
            h1 = h1.remove_row(*x);
            h1 = h1.remove_column(*x);
        }

        let eig = h1.eigenvalues().unwrap();
        let emin = eig.iter().min_by(|a, b| a.total_cmp(b)).unwrap();

        if *emin < self.eps {
            // mix in steepest descent step if Hessian minimum eigenvalue is
            // negative (Experiment (???))
            let adj = self.eps.max(0.01 * emin.abs()) - emin;
            let (r, _) = h.shape();
            h += adj * Dmat::identity(r, r);
        }

        let mut xkd = xk.clone();
        if bhyp {
            for x in &self.forcefield.excision {
                g = g.remove_row(*x);
                h1 = h1.remove_row(*x);
                h1 = h1.remove_column(*x);
                xkd = xkd.remove_row(*x);
            }

            // TODO a bunch of Hyper stuff
            todo!()
        } else {
            let g0 = g.clone();
            let h0 = h.clone();
            for x in &self.forcefield.excision {
                g = g.remove_row(*x);
                h = h.remove_row(*x);
                h = h.remove_column(*x);

                let hi = invert_svd(h.clone());
                let mut dx = -1.0 * hi * &g;

                for i in &self.forcefield.excision {
                    dx = dx.insert_row(*i, 0.0);
                }
            }
        }

        let solver = if bhyp { hyper_solver } else { para_solver };

        // TODO counting micro iterations

        let trust_fun = |l| {
            let n = solver(l).0.norm();
            let d = n - trust;
            d * d
        };

        // evaluate ONLY the objective function. most useful when the objective
        // is cheap, but the derivative is expensive
        let search_fun = |l| {
            // dx is how much the step changes from the previous step
            let (dx, sol) = solver(l);
            // this is our trial step
            let xk_ = dx + &xk;
            let result = self.objective.borrow_mut().full(xk_, 0).x - data.x;
            // TODO if not self.retain_micro_outputs delete output dirs

            // TODO search_fun.micro += 1

            result
        };

        let h_fun = |l| {
            let n = solver(l).0.norm();
            let d = n - self.h;
            d * d
        };

        let mut bump = false;
        let (mut dx, mut expect);
        if self.trust0 > 0.0 {
            // trust region code
            (dx, expect) = solver(1.0);
            let dxnorm = dx.norm();
            if dxnorm > trust {
                bump = true;
                // Tried a few optimizers here, seems like Brent works well.
                // Okay, the problem with Brent is that the tolerance is
                // fractional. If the optimized value is zero, then it takes a
                // lot of meaningless steps.
                let lopt =
                    optimize::brent(trust_fun, self.lmg..4.0 * self.lmg, 1e-6)
                        .0;
                (dx, expect) = solver(lopt);
                let dxnorm = dx.norm();
            } // else we found the step
        } else {
            // search code

            // first obtain a step that is roughly the same length as the
            // provided trust radiusu
            (dx, expect) = solver(1.0);
            let mut dxnorm = dx.norm();
            let mut lopt;
            if dxnorm > trust {
                lopt =
                    optimize::brent(trust_fun, self.lmg..4.0 * self.lmg, 1e-4)
                        .0;
                (dx, expect) = solver(lopt);
                dxnorm = dx.norm();
            } else {
                lopt = 1.0;
            }
            bump = false;
            // TODO micro = 0

            // TODO this time we need the full output
            let result =
                optimize::brent(search_fun, lopt..4.0 * lopt, self.search_tol);
            if result.1 > 0.0 {
                lopt = optimize::brent(h_fun, self.lmg..4.0 * self.lmg, 1e-6).0;
                (dx, expect) = solver(lopt);
                dxnorm = dx.norm();
                // restarting search with step size dxnorm
                let result = optimize::brent(
                    search_fun,
                    lopt..4.0 * lopt,
                    self.search_tol,
                );
            }
            (dx, _) = solver(result.0);
            expect = result.1;
            // optimization step found
        }

        use PenaltyType::*;
        // TODO is_fuse() ?
        if matches!(
            self.objective.borrow().penalty.ptype,
            Fuse | FuseL0 | FuseBarrier
        ) {
            self.forcefield.make_redirect(&dx + xk);
        }

        (dx, expect, bump)
    }
}
