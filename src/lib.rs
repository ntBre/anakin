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
