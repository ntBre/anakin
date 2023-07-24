use std::{error::Error, fs::read_to_string, path::Path};

use serde::Deserialize;

mod default_fns;

use default_fns::*;

/// generated partially by `scripts/convert_config.py`
#[derive(Deserialize)]
#[allow(unused)]
pub struct Config {
    /// Path for GROMACS executables (if not the default)
    gmxpath: Option<String>,

    /// The suffix of GROMACS executables
    gmxsuffix: Option<String>,

    /// Path for TINKER executables (if not the default)
    tinkerpath: Option<String>,

    /// Type of the penalty: L2, L1 or Box
    #[serde(default = "default_penalty_type")]
    penalty_type: String,

    /// Values to scan in the parameter space, given like this: -0.1:0.1:11
    scan_vals: Option<String>,

    /// Name of the restart file we read from
    readchk: Option<String>,

    /// Name of the restart file we write to (can be same as readchk)
    writechk: Option<String>,

    /// Directory containing force fields, relative to project directory
    #[serde(default = "default_ffdir")]
    ffdir: String,

    /// The AMOEBA polarization type, either direct, mutual, or nonpolarizable.
    amoeba_pol: Option<String>,

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

    /// The port number to use for Work Queue
    wq_port: isize,

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
    constrain_charge: bool,

    /// Print the objective function gradient at every step
    #[serde(default = "default_print_gradient")]
    print_gradient: bool,

    /// Optimize in the space of log-variables
    #[serde(default)]
    logarithmic_map: bool,

    /// Print the objective function Hessian at every step
    #[serde(default)]
    print_hessian: bool,

    /// Print the mathematical and physical parameters at every step
    #[serde(default = "default_print_parameters")]
    print_parameters: bool,

    /// Normalize the weights for the fitting targets
    #[serde(default = "default_normalize_weights")]
    normalize_weights: bool,

    /// Set to false to suppress printing options that are equal to their
    /// defaults
    #[serde(default)]
    verbose_options: bool,

    /// Perform calculations using rigid water molecules.
    #[serde(default)]
    rigid_water: bool,

    /// Perform calculations with contrained hydrogen bond lengths.
    #[serde(default)]
    constrain_h: bool,

    /// Generate bonds from virtual sites to host atom bonded atoms.
    #[serde(default)]
    vsite_bonds: bool,

    /// Bypass the transformation matrix and use the physical parameters
    /// directly
    #[serde(default)]
    pub(crate) use_pvals: bool,

    /// Execute Work Queue tasks and local calculations asynchronously for
    /// improved speed
    asynchronous: bool,

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
    duplicate_pnames: bool,
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
