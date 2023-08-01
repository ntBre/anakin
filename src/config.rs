use std::{collections::HashMap, error::Error, fs::read_to_string, path::Path};

use serde::Deserialize;

mod default_fns;

use default_fns::*;

/// generated partially by `scripts/convert_config.py`
#[derive(Deserialize)]
pub struct Config {
    /// Path for GROMACS executables (if not the default)
    gmxpath: Option<String>,

    /// The suffix of GROMACS executables
    gmxsuffix: Option<String>,

    /// Path for TINKER executables (if not the default)
    tinkerpath: Option<String>,

    /// Type of the penalty: L2, L1 or Box
    #[serde(default = "default_penalty_type")]
    pub(crate) penalty_type: String,

    /// Values to scan in the parameter space, given like this: -0.1:0.1:11
    scan_vals: Option<String>,

    /// Name of the restart file we read from
    readchk: Option<String>,

    /// Name of the restart file we write to (can be same as readchk)
    writechk: Option<String>,

    /// Directory containing force fields, relative to project directory
    #[serde(default = "default_ffdir")]
    pub(crate) ffdir: String,

    /// The AMOEBA polarization type, either direct, mutual, or nonpolarizable.
    /// TODO if we use this it should actually be an enum
    pub(crate) amoeba_pol: Option<String>,

    /// Path to AMBER installation directory (leave blank to use AMBERHOME
    /// environment variable.
    amberhome: Option<String>,

    /// Maximum number of steps in an optimization
    #[serde(default = "default_maxstep")]
    maxstep: isize,

    /// Number of good optimization steps to average over when checking the
    /// objective convergence criterion
    #[serde(default = "default_objective_history")]
    objective_history: isize,

    /// The number of convergence criteria that must be met for main optimizer
    /// to converge
    #[serde(default = "default_criteria")]
    criteria: isize,

    /// Number of beads in ring polymer MD (zero to disable)
    rpmd_beads: Option<isize>,

    /// Set to a nonnegative number to turn on zero gradient skipping at that
    /// optimization step.
    #[serde(default = "default_zerograd")]
    zerograd: isize,

    /// Write temp directories to backup before wiping them
    #[serde(default = "default_backup")]
    backup: bool,

    /// Write the checkpoint file at every optimization step
    #[serde(default = "default_writechk_step")]
    writechk_step: bool,

    /// Whether to retain the output files of completed micro iterations
    #[serde(default = "default_retain_micro_outputs")]
    retain_micro_outputs: bool,

    /// Allow convergence on "low quality" steps
    #[serde(default)]
    converge_lowq: bool,

    /// Specify whether there are virtual sites in the simulation (being fitted
    /// or not). Enforces calculation of vsite positions.
    #[serde(default)]
    have_vsite: bool,

    /// Specify whether to constrain the charges on the molecules.
    pub(crate) constrain_charge: bool,

    /// Print the objective function gradient at every step
    #[serde(default = "default_print_gradient")]
    print_gradient: bool,

    /// Optimize in the space of log-variables
    #[serde(default)]
    pub(crate) logarithmic_map: bool,

    /// Print the objective function Hessian at every step
    #[serde(default)]
    print_hessian: bool,

    /// Print the mathematical and physical parameters at every step
    #[serde(default = "default_print_parameters")]
    print_parameters: bool,

    /// Normalize the weights for the fitting targets
    #[serde(default = "default_normalize_weights")]
    pub(crate) normalize_weights: bool,

    /// Set to false to suppress printing options that are equal to their
    /// defaults
    #[serde(default)]
    verbose_options: bool,

    /// Perform calculations using rigid water molecules.
    #[serde(default)]
    pub(crate) rigid_water: bool,

    /// Perform calculations with contrained hydrogen bond lengths.
    #[serde(default)]
    pub(crate) constrain_h: bool,

    /// Generate bonds from virtual sites to host atom bonded atoms.
    #[serde(default)]
    vsite_bonds: bool,

    /// Bypass the transformation matrix and use the physical parameters
    /// directly
    #[serde(default)]
    pub(crate) use_pvals: bool,

    /// Re-evaluate the objective function and gradients when the step is
    /// rejected (for noisy objective functions).
    #[serde(default)]
    reevaluate: bool,

    /// Continue the current run from where we left off (supports mid-iteration
    /// recovery).
    #[serde(rename = "continue")]
    #[serde(default)]
    cont: bool,

    /// Allow duplicate parameter names (only if you know what you are doing!
    #[serde(default)]
    pub(crate) duplicate_pnames: bool,

    /// Levenberg-Marquardt trust radius; set to negative for nonlinear search
    #[serde(default = "default_trust0")]
    trust0: f64,

    /// Minimum trust radius (if the trust radius is tiny, then noisy
    /// optimizations become really gnarly)
    mintrust: f64,

    /// Convergence criterion of objective function (in MainOptimizer this is
    /// the stdev of X2 over [objective_history] steps)
    #[serde(default = "default_convergence_objective")]
    convergence_objective: f64,

    /// Convergence criterion of gradient norm
    #[serde(default = "default_convergence_gradient")]
    convergence_gradient: f64,

    /// Convergence criterion of step size (just needs to fall below this
    /// threshold)
    #[serde(default = "default_convergence_step")]
    convergence_step: f64,

    /// Minimum eigenvalue for applying steepest descent correction
    #[serde(default = "default_eig_lowerbound")]
    eig_lowerbound: f64,

    /// Optimization will "fail" if step falls below this size
    #[serde(default = "default_step_lowerbound")]
    step_lowerbound: f64,

    /// Guess value for bracketing line search in trust radius algorithm
    #[serde(default = "default_lm_guess")]
    lm_guess: f64,

    /// Step size for finite difference derivatives in many functions
    #[serde(default = "default_finite_difference_h")]
    finite_difference_h: f64,

    /// Make sure that the finite difference step size does not exceed this
    /// multiple of the trust radius.
    #[serde(default = "default_finite_difference_factor")]
    finite_difference_factor: f64,

    /// Factor for additive penalty function in objective function
    pub(crate) penalty_additive: f64,

    /// Factor for multiplicative penalty function in objective function
    #[serde(default)]
    pub(crate) penalty_multiplicative: f64,

    /// Extra parameter for fusion penalty function. Dictates position of log
    /// barrier or L1-L0 switch distance
    #[serde(default = "default_penalty_alpha")]
    pub(crate) penalty_alpha: f64,

    /// Cusp region for hyperbolic constraint; for x=0, the Hessian is a/2b
    #[serde(default = "default_penalty_hyperbolic_b")]
    pub(crate) penalty_hyperbolic_b: f64,

    /// Power of the Euclidean norm of the parameter vector (default 2.0 is
    /// normal L2 penalty)
    #[serde(default = "default_penalty_power")]
    pub(crate) penalty_power: f64,

    /// The step size is increased / decreased by up to this much in the event
    /// of a good / bad step; increase for a more variable step size.
    #[serde(default = "default_adaptive_factor")]
    adaptive_factor: f64,

    /// Damping factor that ties down the trust radius to trust0; decrease for a
    /// more variable step size.
    #[serde(default = "default_adaptive_damping")]
    adaptive_damping: f64,

    /// Error tolerance; the optimizer will only reject steps that increase the
    /// objective function by more than this number.
    error_tolerance: f64,

    /// Search tolerance; used only when trust radius is negative, dictates
    /// convergence threshold of nonlinear search.
    #[serde(default = "default_search_tolerance")]
    search_tolerance: f64,

    /// The AMOEBA mutual polarization criterion.
    pub(crate) amoeba_eps: Option<f64>,

    /// The names of force fields, corresponding to files
    /// `forcefields/filename.ext`.
    pub(crate) forcefield: Vec<String>,

    pub(crate) priors: Option<HashMap<String, HashMap<String, f64>>>,
}

impl Config {
    pub fn load<P>(path: P) -> Result<Self, Box<dyn Error>>
    where
        P: AsRef<Path>,
    {
        let contents = read_to_string(path)?;
        let config = toml::from_str(&contents)?;
        Ok(config)
    }
}
