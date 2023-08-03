//! Rust version of the optimization package ForceBalance

#![feature(lazy_cell)]
#![allow(unused)]

type Dvec = nalgebra::DVector<f64>;
type Dmat = nalgebra::DMatrix<f64>;

pub mod config;
pub mod forcefield;
mod molecule;
pub mod objective;
pub mod optimizer;
mod utils;
mod work_queue;
