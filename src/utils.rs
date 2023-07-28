use nalgebra::SVD;

use crate::{Dmat, Dvec};

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

/// invert `h` using singular value decomposition
pub(crate) fn invert_svd(h: Dmat) -> Dmat {
    let svd = SVD::new(h, true, true);
    let u = svd.u.unwrap();
    let s = svd.singular_values;
    let vh = svd.v_t.unwrap();
    let v = vh.transpose();
    let uh = u.transpose();

    let (r, _) = s.shape();
    let mut si = s.clone();
    for i in 0..r {
        if s[i].abs() > 1e-12 {
            si[i] = 1.0 / s[i];
        } else {
            si[i] = 0.0;
        }
    }

    let si = Dmat::from_diagonal(&si);

    v * si * uh
}

/// Given two vectors `v1` and `v2`, project out the component of `v1` that is
/// along the `v2` direction
pub(crate) fn orthogonalize(v1: Dvec, v2: Dvec) -> Dvec {
    let v2u = &v2 / v2.norm();
    &v1 - &v2u * (v1.dot(&v2u))
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use nalgebra::{dmatrix, dvector};

    use super::*;

    #[test]
    fn test_std_dev() {
        let got = std_dev(&[2., 4., 4., 4., 5., 5., 7., 9.]);
        let want = 2.0;
        assert_eq!(got, want);
    }

    #[test]
    fn test_invert_svd() {
        let a = dmatrix![
             0.36697719, -1.66159978,  0.67428787,  0.95471869, -0.56679188, -0.01711798;
            -0.74538269,  1.51043451,  1.86606766,  1.54650471,  1.10091666, -1.22139565;
            -1.74799674, -0.27737085, -0.81011218,  0.75476982,  0.48963132, -0.75768452;
             0.38177729,  0.95066395, -0.95018952, -1.82094117, -1.22058972, -0.76688585;
             0.34602198, -0.36669003, -0.97698248, -0.92653817, -0.74624743, -1.00819961;
            -1.07075401, -0.67128857,  0.80680508,  0.91201023, -1.14635459, 0.52916454;
            -0.29440139, -0.06310497, -0.95573141,  1.01467853,  0.28952005, -0.28664023;
            -0.26080475, -0.18001284, -0.77877666,  1.43770459, -0.72468674, 1.91552729;
             0.22023536,  0.66911524,  0.7576564 ,  1.6363398 ,  2.18313475, -0.87904127;
        ];
        let got = invert_svd(a);
        let want = dmatrix![
            0.21306264, -0.00498933, -0.38024463,  0.10292814,  0.14231749, -0.17653875,  0.11594121,  0.11249802,  0.17857546;
            -0.22821694,  0.23634011, -0.12831318,  0.23931617, -0.03767579, -0.0177474 ,  0.05129079,  0.19568922,  0.02010601;
            0.06142217,  0.13086315, -0.11865457, -0.0845286 , -0.12488938, 0.14349874, -0.22904299, -0.16509635, -0.06246162;
            0.18096369,  0.15882583, -0.07600512,  0.08034847,  0.11870443, 0.042973  ,  0.24201302,  0.25871529,  0.14578433;
            -0.15791544, -0.18541514,  0.13327825, -0.27297548, -0.16557467, -0.20425462, -0.08803939, -0.16444261,  0.12128662;
            -0.19306354, -0.14504148, -0.06447834, -0.16048321, -0.25867453, -0.01536834, -0.1388115 ,  0.14817648, -0.08114219;
        ];
        assert_abs_diff_eq!(got, want, epsilon = 1e-8);
    }

    #[test]
    fn test_orthogonalize() {
        let v1 = dvector![1.0, 2.0, 3.0];
        let v2 = dvector![4.0, 5.0, 6.0];
        let got = orthogonalize(v1, v2);
        let want = dvector![-0.66233766, -0.07792208, 0.50649351];
        assert_abs_diff_eq!(got, want, epsilon = 1e-8);
    }
}
