use std::ops::Range;

/// Given a function of one variable and a bracket, return a local minimizer of
/// the function isolated to a fractional precision of tol. Adapted from
/// `scipy.optimize.brent`
struct Brent<F>
where
    F: Fn(f64) -> f64,
{
    func: F,
    brack: Range<f64>,
    tol: f64,
}

impl<F> Brent<F>
where
    F: Fn(f64) -> f64,
{
    fn new(func: F, brack: Range<f64>, tol: f64) -> Self {
        Self { func, brack, tol }
    }

    fn optimize(self) -> f64 {
        let (xa, xb, xc, fa, fb, fc, mut funcalls) = self.get_bracket_info();
        const MINTOL: f64 = 1e-11;
        const CG: f64 = 0.3819660;

        // BEGIN CORE ALGORITHM
        let mut x = xb;
        let mut w = xb;
        let mut v = xb;
        let mut fx = fb;
        let mut fw = fb;
        let mut fv = fb;
        let (mut a, mut b);
        if xa < xc {
            a = xa;
            b = xc;
        } else {
            a = xc;
            b = xa;
        }
        let mut deltax: f64 = 0.0;
        let mut iter = 0;

        const MAX_ITER: usize = 500;
        while iter < MAX_ITER {
            let tol1 = self.tol * x.abs() * MINTOL;
            let tol2 = 2.0 * tol1;
            let xmid = 0.5 * (a + b);
            if (x - xmid).abs() < (tol2 - 0.5 * (b - a)) {
                break; // convergence achieved
            }
            // apparently this is kind of a bug in scipy's own implementation,
            // but this will only be bound in the true case of the conditional
            // on the first iteration. we use an option to "handle" that, but we
            // could still crash on the unwrap obviously. Their comment is "It
            // should be set before the if (but to what?)" Based on the github
            // issue 4140, it looks like the root issue actually triggering the
            // panic was negative tolerance
            let msg = "rat not set in Brent::optimize. negative tolerance?";
            let mut rat = None;
            let mut u;
            if deltax.abs() <= tol1 {
                if x >= xmid {
                    deltax = a - x; // golden section step
                } else {
                    deltax = b - x;
                }
                rat = Some(CG * deltax);
            } else {
                // parabolic step
                let tmp1 = (x - w) * (fx - fv);
                let tmp2 = (x - v) * (fx - fw);
                let mut p = (x - v) * tmp2 - (x - w) * tmp1;
                let mut tmp2 = 2.0 * (tmp2 - tmp1);
                if (tmp2 > 0.0) {
                    p = -p
                }
                tmp2 = tmp2.abs();
                let dx_temp = deltax;
                deltax = rat.expect(msg);
                // check parabolic fit
                if ((p > tmp2 * (a - x))
                    && (p < tmp2 * (b - x))
                    && (p.abs() < f64::abs(0.5 * tmp2 * dx_temp)))
                {
                    rat = Some(p * 1.0 / tmp2); // if parabolic step is useful.
                    let u = x + rat.expect(msg);
                    if ((u - a) < tol2 || (b - u) < tol2) {
                        if xmid - x >= 0.0 {
                            rat = Some(tol1)
                        } else {
                            rat = Some(-tol1)
                        }
                    }
                } else if (x >= xmid) {
                    deltax = a - x; // if it's not do a golden section step
                } else {
                    deltax = b - x;
                    rat = Some(CG * deltax);
                }
            }
            if rat.expect(msg).abs() < tol1 {
                if rat.unwrap() >= 0.0 {
                    u = x + tol1;
                } else {
                    u = x - tol1;
                }
            } else {
                u = x + rat.unwrap();
            }
            let fu = (self.func)(u);
            funcalls += 1;

            // if it's bigger than current
            if fu > fx {
                if u < x {
                    a = u;
                } else {
                    b = u;
                }
                if fu <= fw || w == x {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                } else if fu <= fv || v == x || v == w {
                    v = u;
                    fv = fu;
                }
            } else {
                if u >= x {
                    a = x
                } else {
                    b = x
                }
                v = w;
                w = x;
                x = u;
                fv = fw;
                fw = fx;
                fx = fu;
            }

            iter += 1;
        }
        // END CORE ALGORITHM

        x
    }

    fn get_bracket_info(&self) -> (f64, f64, f64, f64, f64, f64, usize) {
        bracket(&self.func, self.brack.start, self.brack.end)
    }
}

/// Bracket the minimum of a function. Given a function `func` and distinct
/// initial points `xa` and `xb`, search in the downhill direction (as defined
/// by the initial points) and return three points that bracket the minimum of
/// the function. Adapted from `scipy/optimize`.
fn bracket<F>(
    func: F,
    mut xa: f64,
    mut xb: f64,
) -> (f64, f64, f64, f64, f64, f64, usize)
where
    F: Fn(f64) -> f64,
{
    const GOLD: f64 = 1.618034;
    const VERYSMALL: f64 = 1e-21;
    const GROW_LIMIT: f64 = 110.0;
    const MAX_ITER: usize = 500;

    let mut fa = func(xa);
    let mut fb = func(xb);
    if fa < fb {
        (xa, xb) = (xb, xa);
        (fa, fb) = (fb, fa);
    }
    let mut xc = xb + GOLD * (xb - xa);
    let mut fc = func(xc);
    let mut funcalls = 3;

    let mut iter = 0;
    while fc < fb {
        let tmp1 = (xb - xa) * (fb - fc);
        let tmp2 = (xb - xc) * (fb - fa);
        let val = tmp2 - tmp1;
        let denom = if val.abs() < VERYSMALL {
            2.0 * VERYSMALL
        } else {
            2.0 * val
        };
        let w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom;
        let wlim = xb + GROW_LIMIT * (xc - xb);
        if iter > MAX_ITER {
            // TODO this should return a result. in Python you can try/except
            // the raise here unlike a panic
            panic!(
                "No valid bracket was found before the iteration limit was \
		 reached. Consider trying different initial points or \
		 increasing MAX_ITER."
            );
        }
        iter += 1;
        let mut fw;
        if (w - xc) * (xb - w) > 0.0 {
            fw = func(w);
            funcalls += 1;
            if (fw < fc) {
                xa = xb;
                xb = w;
                fa = fb;
                fb = fw;
                break;
            } else if (fw > fb) {
                xc = w;
                fc = fw;
                break;
            }
            let w = xc + GOLD * (xc - xb);
            fw = func(w);
            funcalls += 1;
        } else if (w - wlim) * (wlim - xc) >= 0.0 {
            let w = wlim;
            fw = func(w);
            funcalls += 1;
        } else if (w - wlim) * (xc - w) > 0.0 {
            fw = func(w);
            funcalls += 1;
            if (fw < fc) {
                xb = xc;
                xc = w;
                let w = xc + GOLD * (xc - xb);
                fb = fc;
                fc = fw;
                fw = func(w);
                funcalls += 1;
            }
        } else {
            let w = xc + GOLD * (xc - xb);
            fw = func(w);
            funcalls += 1;
        }
        xa = xb;
        xb = xc;
        xc = w;
        fa = fb;
        fb = fc;
        fc = fw;
    }

    // three conditions for a valid bracket
    let cond1 = (fb < fc && fb <= fa) || (fb < fa && fb <= fc);
    let cond2 = (xa < xb && xb < xc) || (xc < xb && xb < xa);
    let cond3 = xa.is_finite() && xb.is_finite() && xc.is_finite();
    if !(cond1 && cond2 && cond3) {
        // TODO another result. This one actually returns the data in the error
        panic!(
            "The algorithm terminated without finding a valid bracket. Consider \
	     trying different initial points"
        );
    }
    (xa, xb, xc, fa, fb, fc, funcalls)
}
