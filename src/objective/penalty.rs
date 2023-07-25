//! Penalty functions for regularizing the force field optimizer.
//!
//! The purpose for this module is to improve the behavior of our optimizer;
//! essentially, our problem is fraught with 'linear dependencies', a.k.a.
//! directions in the parameter space that the objective function does not
//! respond to. This would happen if a parameter is just plain useless, or if
//! there are two or more parameters that describe the same thing.
//!
//! To accomplish these objectives, a penalty function is added to the objective
//! function. Generally, the more the parameters change (i.e. the greater the
//! norm of the parameter vector), the greater the penalty. Note that this is
//! added on after all of the other contributions have been computed. This only
//! matters if the penalty 'multiplies' the objective function: Obj +
//! Obj*Penalty, but we also have the option of an additive penalty: Obj +
//! Penalty.
//!
//! Statistically, this is called regularization. If the penalty function is the
//! norm squared of the parameter vector, it is called ridge regression. There
//! is also the option of using simply the norm, and this is called lasso, but I
//! think it presents problems for the optimizer that I need to work out.
//!
//! Note that the penalty functions can be considered as part of a 'maximum
//! likelihood' framework in which we assume a PRIOR PROBABILITY of the force
//! field parameters around their initial values. The penalty function is
//! related to the prior by an exponential. Ridge regression corresponds to a
//! Gaussian prior and lasso corresponds to an exponential prior. There is also
//! 'elastic net regression' which interpolates between Gaussian and exponential
//! using a tuning parameter.
//!
//! Our priors are adjustable too - there is one parameter, which is the width
//! of the distribution. We can even use a noninformative prior for the
//! distribution widths (hyperprior!). These are all important things to
//! consider later.
//!
//! Importantly, note that here there is no code that treats the distribution
//! width. That is because the distribution width is wrapped up in the rescaling
//! factors, which is essentially a coordinate transformation on the parameter
//! space. More documentation on this will follow, perhaps in the 'rsmake'
//! method.

use nalgebra::Dyn;

use crate::{forcefield::FF, Dmat, Dvec};

use super::{Extra, ObjMap};

#[derive(Default)]
#[repr(u8)]
pub(crate) enum PenaltyType {
    // TODO serde alias this to hyp, hyper, l1, and hyperbola if I end up
    // loading straight from file. otherwise I might need to break all of those
    // (and change the repr) out if they are used differently later on
    Hyperbolic = 1,

    // aliases: para, parabola, l2, quadratic
    #[default]
    Parabolic = 2,

    Box = 3,

    Fuse = 4,

    FuseL0 = 5,

    FuseBarrier = 6,
}

#[derive(Default)]
pub(crate) struct Penalty {
    fadd: f64,
    fmul: f64,
    a: f64,
    b: f64,
    p: f64,
    pub(crate) ptype: PenaltyType,
    // spacings: HashMap?
}

impl Penalty {
    pub(crate) fn compute(&self, mvals: Dvec, objective: &mut ObjMap) -> Extra {
        // TODO borrow?
        let (k0, k1, k2) = self.pen_tab(mvals.clone());
        let (mut xadd, mut gadd, mut hadd);
        if self.fadd > 0.0 {
            xadd = k0 * self.fadd;
            gadd = &k1 * self.fadd;
            hadd = &k2 * self.fadd;
        } else {
            let np = mvals.len();
            xadd = 0.0;
            gadd = Dvec::zeros(np);
            hadd = Dmat::zeros(np, np);
        }

        if self.fmul > 0.0 {
            let ObjMap { x, g, h, .. } = objective;
            xadd += (*x * k0) * self.fmul;
            gadd += (k0 * g.clone() + *x * &k1) * self.fmul;
            let (rg, cg) = g.shape();
            let (rh, ch) = h.shape();
            // seems weird to use rg/cg on k1 but they should be the same for
            // the matmul to work
            let gk1 = g.clone().reshape_generic(Dyn(1), Dyn(rg * cg))
                * k1.clone().reshape_generic(Dyn(rg * cg), Dyn(1));
            let k1g = k1.clone().reshape_generic(Dyn(1), Dyn(rg * cg))
                * g.clone().reshape_generic(Dyn(rg * cg), Dyn(1));
            hadd += (k0 * h.clone() + gk1 + k1g + *x * k2) * self.fmul;
        }

        Extra(xadd, gadd, hadd)
    }

    // TODO needs a better name. forcebalance stores this as a constant field...
    // just dispatching based on penalty type
    fn pen_tab(&self, vals: Dvec) -> (f64, Dvec, Dmat) {
        match self.ptype {
            PenaltyType::Hyperbolic => self.hyp(vals),
            PenaltyType::Parabolic => self.l2_norm(vals),
            PenaltyType::Box => self.boxed(vals),
            PenaltyType::Fuse => self.fuse(vals),
            PenaltyType::FuseL0 => self.fuse_l0(vals),
            PenaltyType::FuseBarrier => self.fuse_barrier(vals),
        }
    }

    fn hyp(&self, vals: Dvec) -> (f64, Dvec, Dmat) {
        todo!()
    }

    /// Harmonic L2-norm constraints. `mvals` is the parameter vector. Returns
    /// the norm squared of the vector, the gradient of the norm, and the
    /// Hessian (just a constant)
    fn l2_norm(&self, mvals: Dvec) -> (f64, Dvec, Dmat) {
        let l = mvals.len();
        let mvals = Dvec::from(mvals);
        let (dc0, dc1, dc2);
        if self.p == 2.0 {
            dc0 = mvals.dot(&mvals);
            dc1 = 2.0 * mvals;
            dc2 = 2.0 * Dmat::identity(l, l);
        } else {
            let m2 = mvals.dot(&mvals);
            let p = self.p;
            dc0 = m2.powf(p / 2.0);
            dc1 = p * (m2.powf(p / 2.0 - 1.0)) * &mvals;
            dc2 = p * (m2.powf(p / 2.0 - 1.0)) * Dmat::identity(l, l)
                + p * (p - 2.0)
                    * (m2.powf(p / 2.0 - 2.0))
                    * (&mvals * mvals.transpose());
        }
        (dc0, dc1, dc2)
    }

    /// renamed from `box` since that's a Rust keyword
    fn boxed(&self, vals: Dvec) -> (f64, Dvec, Dmat) {
        todo!()
    }

    fn fuse(&self, vals: Dvec) -> (f64, Dvec, Dmat) {
        todo!()
    }

    fn fuse_l0(&self, vals: Dvec) -> (f64, Dvec, Dmat) {
        todo!()
    }

    fn fuse_barrier(&self, vals: Dvec) -> (f64, Dvec, Dmat) {
        todo!()
    }
}
