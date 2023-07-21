//! Rust version of the optimization package ForceBalance

type Dvec = nalgebra::DVector<f64>;
type Dmat = nalgebra::DMatrix<f64>;

#[allow(unused)]
pub mod forcefield;

#[allow(unused)]
pub mod objective;

#[allow(unused)]
pub mod optimizer;

mod work_queue;

mod utils {
    /// compute the population standard deviation of `hist`
    pub(crate) fn std_dev(hist: &[f64]) -> f64 {
        let count = hist.len() as f64;
        let avg = hist.iter().sum::<f64>() / count;
        let variance = hist
            .iter()
            .map(|&v| {
                let diff = v - avg;
                diff * diff
            })
            .sum::<f64>()
            / count;
        variance.sqrt()
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_std_dev() {
            let got = std_dev(&[2., 4., 4., 4., 5., 5., 7., 9.]);
            let want = 2.0;
        }
    }
}
