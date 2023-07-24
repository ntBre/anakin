use serde::Deserialize;

/// generated partially by `scripts/convert_config.py`
#[derive(Deserialize)]
struct Config {
    /// Path for GROMACS executables (if not the default)
    gmxpath: String,

    /// The suffix of GROMACS executables
    gmxsuffix: String,

    /// Path for TINKER executables (if not the default)
    tinkerpath: String,

    /// Type of the penalty: L2, L1 or Box
    #[serde(default = "default_penalty_type")]
    penalty_type: String,

    /// Values to scan in the parameter space, given like this: -0.1:0.1:11
    scan_vals: String,

    /// Name of the restart file we read from
    readchk: String,

    /// Name of the restart file we write to (can be same as readchk)
    writechk: String,

    /// Directory containing force fields, relative to project directory
    #[serde(default = "default_ffdir")]
    ffdir: String,

    /// The AMOEBA polarization type, either direct, mutual, or nonpolarizable.
    amoeba_pol: String,

    /// Path to AMBER installation directory (leave blank to use AMBERHOME
    /// environment variable.
    amberhome: String,

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
    rpmd_beads: isize,

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
    converge_lowq: bool,

    /// Specify whether there are virtual sites in the simulation (being fitted
    /// or not). Enforces calculation of vsite positions.
    have_vsite: bool,

    /// Specify whether to constrain the charges on the molecules.
    constrain_charge: bool,

    /// Print the objective function gradient at every step
    #[serde(default = "default_print_gradient")]
    print_gradient: bool,

    /// Optimize in the space of log-variables
    logarithmic_map: bool,

    /// Print the objective function Hessian at every step
    print_hessian: bool,

    /// Print the mathematical and physical parameters at every step
    #[serde(default = "default_print_parameters")]
    print_parameters: bool,

    /// Normalize the weights for the fitting targets
    #[serde(default = "default_normalize_weights")]
    normalize_weights: bool,

    /// Set to false to suppress printing options that are equal to their
    /// defaults
    verbose_options: bool,

    /// Perform calculations using rigid water molecules.
    rigid_water: bool,

    /// Perform calculations with contrained hydrogen bond lengths.
    constrain_h: bool,

    /// Generate bonds from virtual sites to host atom bonded atoms.
    vsite_bonds: bool,

    /// Bypass the transformation matrix and use the physical parameters
    /// directly
    use_pvals: bool,

    /// Execute Work Queue tasks and local calculations asynchronously for
    /// improved speed
    asynchronous: bool,

    /// Re-evaluate the objective function and gradients when the step is
    /// rejected (for noisy objective functions).
    reevaluate: bool,

    /// Continue the current run from where we left off (supports mid-iteration
    /// recovery).
    #[serde(rename = "continue")]
    cont: bool,

    /// Allow duplicate parameter names (only if you know what you are doing!
    duplicate_pnames: bool,
}

fn default_penalty_type() -> String {
    String::from("L2")
}

fn default_ffdir() -> String {
    String::from("forcefield")
}

fn default_maxstep() -> isize {
    100
}

fn default_objective_history() -> isize {
    2
}

fn default_criteria() -> isize {
    1
}

fn default_zerograd() -> isize {
    -1
}

fn default_backup() -> bool {
    true
}

fn default_writechk_step() -> bool {
    true
}

fn default_retain_micro_outputs() -> bool {
    true
}

fn default_print_gradient() -> bool {
    true
}

fn default_print_parameters() -> bool {
    true
}

fn default_normalize_weights() -> bool {
    true
}
