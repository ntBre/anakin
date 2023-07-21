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

use crate::forcefield::FF;

use super::{Extra, ObjMap};

#[derive(Default)]
#[repr(u8)]
enum PenaltyType {
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
    ff: FF,
    ptype: PenaltyType,
    // spacings: HashMap?
}

impl Penalty {
    pub(crate) fn compute(
        &self,
        vals: Vec<f64>,
        objective: &mut ObjMap,
    ) -> Extra {
        todo!()
    }
}
