use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub(crate) enum TargetType {
    #[serde(alias = "TorsionProfile_SMIRNOFF")]
    Torsion,
    #[serde(alias = "OptGeoTarget_SMIRNOFF")]
    OptGeo,
}

/// struct corresponding to tgt_opts in ForceBalance
#[derive(Debug, Deserialize)]
pub(crate) struct Target {
    pub(crate) name: String,

    #[serde(rename = "type")]
    pub(crate) typ: TargetType,

    /// The resolution of mapping interactions to net forces and torques for
    /// groups of atoms. In order of resolution: molecule > residue >
    /// charge-group
    #[serde(default = "default_force_map")]
    force_map: String,

    /// Interaction fragment 1: a selection of atoms specified using atoms and
    /// dashes, e.g. 1-6 to select the first through sixth atom (i.e. list
    /// numbering starts from 1)
    #[serde(default)]
    fragment1: String,

    /// Interaction fragment 2: a selection of atoms specified using atoms and
    /// dashes, e.g. 7-11 to select atoms 7 through 11.
    #[serde(default)]
    fragment2: String,

    /// Precision of OpenMM calculation if using CUDA or OpenCL platform. Choose
    /// either single, double or mixed ; defaults to the OpenMM default.
    #[serde(default)]
    pub(crate) openmm_precision: Option<String>,

    /// OpenMM platform. Choose either Reference, CUDA or OpenCL. AMOEBA is on
    /// Reference or CUDA only.
    #[serde(default)]
    pub(crate) openmm_platform: Option<String>,

    /// Text file containing quantum data. If not provided, will search for a
    /// default (qdata.txt).
    #[serde(default)]
    qdata_txt: Option<String>,

    /// Text file containing interacting systems. If not provided, will search
    /// for a default.
    #[serde(default = "default_inter_txt")]
    inter_txt: String,

    /// Text file containing extra options for Optimized Geometry target. If not
    /// provided, will search for a default.
    #[serde(default = "default_optgeo_options_txt")]
    optgeo_options_txt: String,

    /// Reassign modes before fitting frequencies, using either linear
    /// assignment "permute" or maximum overlap "overlap".
    #[serde(default)]
    reassign_modes: Option<String>,

    /// Provide file name for condensed phase coordinates.
    #[serde(default)]
    liquid_coords: Option<String>,

    /// Provide file name for gas phase coordinates.
    #[serde(default)]
    gas_coords: Option<String>,

    /// Provide file name for condensed phase NVT coordinates.
    #[serde(default)]
    nvt_coords: Option<String>,

    /// Provide file name for lipid coordinates.
    #[serde(default)]
    lipid_coords: Option<String>,

    /// Coordinates for single point evaluation; if not provided, will search
    /// for a default.
    #[serde(default)]
    pub(crate) coords: Option<String>,

    /// PDB file mainly used for building OpenMM and AMBER systems.
    #[serde(default)]
    pub(crate) pdb: Option<String>,

    /// MOL2 file needed to set up the system (in addition to any specified
    /// under forcefield). NOTE: this is a list in Python, but we seem to use
    /// only single ones so far
    #[serde(default)]
    pub(crate) mol2: Option<String>,

    /// Gromacs .mdp files. If not provided, will search for default.
    #[serde(default)]
    gmx_mdp: Option<String>,

    /// Gromacs .top files. If not provided, will search for default.
    #[serde(default)]
    gmx_top: Option<String>,

    /// Gromacs .ndx files. If not provided, will search for default.
    #[serde(default)]
    gmx_ndx: Option<String>,

    /// File containing commands for "tleap" when setting up AMBER simulations.
    #[serde(default)]
    amber_leapcmd: Option<String>,

    /// TINKER .key files. If not provided, will search for default.
    #[serde(default)]
    tinker_key: Option<String>,

    /// Text file containing experimental data.
    #[serde(default = "default_expdata_txt")]
    expdata_txt: String,

    /// Text file containing experimental data.
    #[serde(default = "default_hfedata_txt")]
    hfedata_txt: String,

    /// Method for calculating hydration energies (single point, FEP, TI).
    #[serde(default = "default_hfemode")]
    hfemode: String,

    /// Provide a temporary directory ".tmp" to read data from a previous
    /// calculation on the initial iteration (for instance, to restart an
    /// aborted run).
    #[serde(default)]
    pub(crate) read: Option<String>,

    /// Specify an optional prefix script to run in front of rtarget.py, for
    /// loading environment variables
    #[serde(default)]
    remote_prefix: String,

    /// Number of fitting atoms; defaults to all of them. Use a comma and dash
    /// style list (1,2-5), atoms numbered from one, inclusive
    #[serde(default)]
    fitatoms: Option<String>,

    /// Specify a subset of molecules to fit. The rest are used for
    /// cross-validation.
    #[serde(default)]
    subset: Option<String>,

    /// Name of the barostat to use for equilibration.
    #[serde(default = "default_gmx_eq_barostat")]
    gmx_eq_barostat: String,

    /// JSON file containing options for the OpenFF Evaluator target. If not
    /// provided, will search for a default.
    #[serde(default = "default_evaluator_input")]
    evaluator_input: String,

    /// A SQLite database containing the electrostatic property data to train
    /// against.
    #[serde(default = "default_recharge_esp_store")]
    recharge_esp_store: String,

    /// The type of electrostatic property to train against [esp,
    /// electric-field].
    #[serde(default = "default_recharge_property")]
    recharge_property: String,

    /// How to treat relative (MM-QM) energies. "average": Subtract out the mean
    /// gap (default). "qm_minimum": Reference all MM and QM energies to the
    /// structure with minimum QM energy. "absolute": Use absolute energies in
    /// fitting, do not subtract out any energy gap.
    #[serde(default = "default_energy_mode")]
    energy_mode: String,

    /// Number of snapshots; defaults to all of the snapshots
    #[serde(default)]
    shots: Option<usize>,

    /// Wait a number of seconds every time this target is visited (gives me a
    /// chance to ctrl+C)
    #[serde(default = "default_sleepy")]
    pub(crate) sleepy: usize,

    /// Number of time steps for the liquid production run.
    #[serde(default = "default_liquid_md_steps")]
    liquid_md_steps: usize,

    /// Number of time steps for the liquid equilibration run.
    #[serde(default = "default_liquid_eq_steps")]
    liquid_eq_steps: usize,

    /// Number of time steps for the lipid production run.
    #[serde(default = "default_lipid_md_steps")]
    lipid_md_steps: usize,

    /// Number of time steps for the lipid equilibration run.
    #[serde(default = "default_lipid_eq_steps")]
    lipid_eq_steps: usize,

    /// Number of steps in the liquid simulation between MC barostat volume
    /// adjustments.
    #[serde(default = "default_n_mcbarostat")]
    n_mcbarostat: usize,

    /// Number of time steps for the gas production run, if different from
    /// default.
    #[serde(default = "default_gas_md_steps")]
    gas_md_steps: usize,

    /// Number of time steps for the gas equilibration run, if different from
    /// default.
    #[serde(default = "default_gas_eq_steps")]
    gas_eq_steps: usize,

    /// Number of time steps for the liquid NVT production run.
    #[serde(default = "default_nvt_md_steps")]
    nvt_md_steps: usize,

    /// Number of time steps for the liquid NVT equilibration run.
    #[serde(default = "default_nvt_eq_steps")]
    nvt_eq_steps: usize,

    /// Affects the amount of data being printed to the temp directory.
    #[serde(default = "default_writelevel")]
    pub(crate) writelevel: usize,

    /// Set the number of threads used by Gromacs or TINKER processes in MD
    /// simulations
    #[serde(default = "default_md_threads")]
    md_threads: usize,

    /// Whether to save trajectories. 0 = Never save; 1 = Delete if optimization
    /// step is good; 2 = Always save
    #[serde(default = "default_save_traj")]
    save_traj: usize,

    /// Number of time steps for the equilibration run.
    #[serde(default = "default_eq_steps")]
    eq_steps: usize,

    /// Number of time steps for the production run.
    #[serde(default = "default_md_steps")]
    md_steps: usize,

    /// Number of simulations required to calculate quantities.
    #[serde(default = "default_n_sim_chain")]
    n_sim_chain: usize,

    /// Provide the number of molecules in the structure (defaults to
    /// auto-detect).
    #[serde(default)]
    n_molecules: Option<usize>,

    /// Specify an hessian target objective function normalization method.
    #[serde(default = "default_hess_normalize_type")]
    hess_normalize_type: usize,

    /// Weight of the target (determines its importance vs. other targets)
    #[serde(default = "default_weight")]
    pub(crate) weight: f64,

    /// Energy denominator in kcal/mol for objective function and lower energy
    /// limit for attenuating weights where applicable
    #[serde(default = "default_energy_denom")]
    pub(crate) energy_denom: f64,

    /// Upper energy cutoff in kcal/mol for setting weights to zero, used to
    /// exclude super-repulsive interactions
    #[serde(default = "default_energy_upper")]
    pub(crate) energy_upper: f64,

    /// Weight of experimental density
    #[serde(default = "default_w_rho")]
    w_rho: f64,

    /// Weight of enthalpy of vaporization
    #[serde(default = "default_w_hvap")]
    w_hvap: f64,

    /// Weight of thermal expansion coefficient
    #[serde(default = "default_w_alpha")]
    w_alpha: f64,

    /// Weight of isothermal compressibility
    #[serde(default = "default_w_kappa")]
    w_kappa: f64,

    /// Weight of isobaric heat capacity
    #[serde(default = "default_w_cp")]
    w_cp: f64,

    /// Weight of dielectric constant
    #[serde(default = "default_w_eps0")]
    w_eps0: f64,

    /// Weight of average area per lipid
    #[serde(default = "default_w_al")]
    w_al: f64,

    /// Weight of deuterium order parameter
    #[serde(default = "default_w_scd")]
    w_scd: f64,

    /// Weight of energy
    #[serde(default = "default_w_energy")]
    w_energy: f64,

    /// Weight of atomistic forces
    #[serde(default = "default_w_force")]
    w_force: f64,

    /// Weight of surface tension
    #[serde(default = "default_w_surf_ten")]
    w_surf_ten: f64,

    /// Weight of net forces (condensed to molecules, residues, or charge
    /// groups)
    #[serde(default = "default_w_netforce")]
    w_netforce: f64,

    /// Weight of torques (condensed to molecules, residues, or charge groups)
    #[serde(default = "default_w_torque")]
    w_torque: f64,

    /// Weight of RESP
    #[serde(default = "default_w_resp")]
    w_resp: f64,

    /// a
    #[serde(default = "default_resp_a")]
    resp_a: f64,

    /// b
    #[serde(default = "default_resp_b")]
    resp_b: f64,

    /// Simulation temperature for hydration free energies (Kelvin)
    #[serde(default = "default_hfe_temperature")]
    hfe_temperature: f64,

    /// Simulation temperature for hydration free energies (atm)
    #[serde(default = "default_hfe_pressure")]
    hfe_pressure: f64,

    /// If nonzero, override the Energy RMS used to normalize the energy part of
    /// the objective function term
    #[serde(default = "default_energy_rms_override")]
    energy_rms_override: f64,

    /// If nonzero, override the Force RMS used to normalize the energy part of
    /// the objective function term
    #[serde(default = "default_force_rms_override")]
    force_rms_override: f64,

    /// RMSD normalization for optimized geometries in Angstrom
    #[serde(default = "default_rmsd_denom")]
    rmsd_denom: f64,

    /// Frequency normalization (in wavenumber) for vibrational frequencies
    #[serde(default = "default_wavenumber_tol")]
    wavenumber_tol: f64,

    /// Dipole normalization (Debye) ; set to 0 if a zero weight is desired
    #[serde(default = "default_dipole_denom")]
    dipole_denom: f64,

    /// Quadrupole normalization (Buckingham) ; set to 0 if a zero weight is
    /// desired
    #[serde(default = "default_quadrupole_denom")]
    quadrupole_denom: f64,

    /// Dipole polarizability tensor normalization (cubic Angstrom) ; set to 0
    /// if a zero weight is desired
    #[serde(default = "default_polarizability_denom")]
    polarizability_denom: f64,

    /// Time step size for the liquid simulation.
    #[serde(default = "default_liquid_timestep")]
    liquid_timestep: f64,

    /// Time interval for saving coordinates for the liquid production run.
    #[serde(default = "default_liquid_interval")]
    liquid_interval: f64,

    /// Time step size for the gas simulation (if zero, use default in external
    /// script.).
    #[serde(default = "default_gas_timestep")]
    gas_timestep: f64,

    /// Time interval for saving coordinates for the gas production run (if
    /// zero, use default in external script.)
    #[serde(default = "default_gas_interval")]
    gas_interval: f64,

    /// Time step size for the lipid simulation.
    #[serde(default = "default_lipid_timestep")]
    lipid_timestep: f64,

    /// Time interval for saving coordinates for the lipid production run.
    #[serde(default = "default_lipid_interval")]
    lipid_interval: f64,

    /// Time step size for the NVT simulation.
    #[serde(default = "default_nvt_timestep")]
    nvt_timestep: f64,

    /// Time interval for saving coordinates for the NVT simulation production
    /// run.
    #[serde(default = "default_nvt_interval")]
    nvt_interval: f64,

    /// Gas-phase dipole parameter for self-polarization correction (in debye).
    #[serde(default = "default_self_pol_mu0")]
    self_pol_mu0: f64,

    /// Polarizability parameter for self-polarization correction (in debye).
    #[serde(default = "default_self_pol_alpha")]
    self_pol_alpha: f64,

    /// Gradient below this threshold will be set to zero.
    #[serde(default = "default_epsgrad")]
    pub(crate) epsgrad: f64,

    /// Snapshots with (E_MM - E_QM) < 0.0 will have their weights increased by
    /// this factor. Only valid if energy_mode is set to "qm_minimum".
    #[serde(default = "default_energy_asymmetry")]
    energy_asymmetry: f64,

    /// Cutoff for nonbonded interactions (passed to engines).
    #[serde(default)]
    nonbonded_cutoff: Option<f64>,

    /// Cutoff for vdW interactions if different from other nonbonded
    /// interactions
    #[serde(default)]
    vdw_cutoff: Option<f64>,

    /// Step size for finite difference derivatives for liquid targets in
    /// pure_num_grad
    #[serde(default = "default_liquid_fdiff_h")]
    liquid_fdiff_h: f64,

    /// Force constant for harmonic positional energy restraints
    #[serde(default = "default_restrain_k")]
    pub(crate) restrain_k: f64,

    /// Finite difference gradient of objective function w/r.t. specified
    /// parameters
    #[serde(default)]
    pub(crate) fdgrad: bool,

    /// Finite difference Hessian of objective function w/r.t. specified
    /// parameters
    #[serde(default)]
    pub(crate) fdhess: bool,

    /// Finite difference Hessian diagonals w/r.t. specified parameters (costs
    /// 2np times a objective calculation)
    #[serde(default)]
    pub(crate) fdhessdiag: bool,

    /// Compute all energies and forces in one fell swoop where possible(as
    /// opposed to calling the simulation code once per snapshot)
    #[serde(default = "default_all_at_once")]
    all_at_once: bool,

    /// For OpenMM or other codes with Python interface: Compute energies and
    /// forces internally
    #[serde(default = "default_run_internal")]
    run_internal: bool,

    /// Enable the energy objective function
    #[serde(default = "default_energy")]
    energy: bool,

    /// Enable the force objective function
    #[serde(default = "default_force")]
    force: bool,

    /// Enable the RESP objective function
    #[serde(default)]
    resp: bool,

    /// Call Q-Chem to do MM COSMO on MM snapshots.
    #[serde(default)]
    do_cosmo: bool,

    /// Perform a geometry optimization before computing properties
    #[serde(default = "default_optimize_geometry")]
    optimize_geometry: bool,

    /// Normalize interaction energies using 1/sqrt(denom**2 + (E(qm)-denom)**2)
    /// for energies more positive than denom.
    #[serde(default)]
    pub(crate) attenuate: bool,

    /// Divide objective function by the number of snapshots / vibrations
    #[serde(default)]
    normalize: bool,

    /// Normalize the condensed phase property contributions to the liquid /
    /// lipid property target
    #[serde(default)]
    w_normalize: bool,

    /// Give the user a chance to fill in condensed phase stuff on the zeroth
    /// step
    #[serde(default)]
    manual: bool,

    /// Don't target the average enthalpy of vaporization and allow it to freely
    /// float (experimental)
    #[serde(default)]
    hvap_subaverage: bool,

    /// Force the external npt.py script to crash if CUDA Platform not available
    #[serde(default)]
    force_cuda: bool,

    /// Enable anisotropic box scaling (e.g. for crystals or two-phase
    /// simulations) in external npt.py script
    #[serde(default)]
    anisotropic_box: bool,

    /// Enable multiple-timestep integrator in external npt.py script
    #[serde(default)]
    mts_integrator: bool,

    /// Minimize the energy of the system prior to running dynamics
    #[serde(default = "default_minimize_energy")]
    minimize_energy: bool,

    /// Evaluate target as a remote work_queue task
    #[serde(default)]
    remote: bool,

    /// Adapt to simulation uncertainty by combining property estimations and
    /// adjusting simulation length.
    #[serde(default)]
    adapt_errors: bool,

    /// When running remote target, back up files at the remote location.
    #[serde(default)]
    remote_backup: bool,

    /// Pure numerical gradients -- launch two additional simulations for each
    /// perturbed forcefield parameter, and compute derivatives using 3-point
    /// formula. (This is very expensive and should only serve as a sanity
    /// check)
    #[serde(default)]
    pure_num_grad: bool,
}
fn default_weight() -> f64 {
    1.0
}

fn default_energy_denom() -> f64 {
    1.0
}

fn default_energy_upper() -> f64 {
    30.0
}

fn default_w_rho() -> f64 {
    1.0
}

fn default_w_hvap() -> f64 {
    1.0
}

fn default_w_alpha() -> f64 {
    1.0
}

fn default_w_kappa() -> f64 {
    1.0
}

fn default_w_cp() -> f64 {
    1.0
}

fn default_w_eps0() -> f64 {
    1.0
}

fn default_w_al() -> f64 {
    1.0
}

fn default_w_scd() -> f64 {
    1.0
}

fn default_w_energy() -> f64 {
    1.0
}

fn default_w_force() -> f64 {
    1.0
}

fn default_w_surf_ten() -> f64 {
    0.0
}

fn default_w_netforce() -> f64 {
    0.0
}

fn default_w_torque() -> f64 {
    0.0
}

fn default_w_resp() -> f64 {
    0.0
}

fn default_resp_a() -> f64 {
    0.001
}

fn default_resp_b() -> f64 {
    0.1
}

fn default_hfe_temperature() -> f64 {
    298.15
}

fn default_hfe_pressure() -> f64 {
    1.0
}

fn default_energy_rms_override() -> f64 {
    0.0
}

fn default_force_rms_override() -> f64 {
    0.0
}

fn default_rmsd_denom() -> f64 {
    0.1
}

fn default_wavenumber_tol() -> f64 {
    10.0
}

fn default_dipole_denom() -> f64 {
    1.0
}

fn default_quadrupole_denom() -> f64 {
    1.0
}

fn default_polarizability_denom() -> f64 {
    1.0
}

fn default_liquid_timestep() -> f64 {
    1.0
}

fn default_liquid_interval() -> f64 {
    0.1
}

fn default_gas_timestep() -> f64 {
    1.0
}

fn default_gas_interval() -> f64 {
    0.1
}

fn default_lipid_timestep() -> f64 {
    1.0
}

fn default_lipid_interval() -> f64 {
    0.1
}

fn default_nvt_timestep() -> f64 {
    1.0
}

fn default_nvt_interval() -> f64 {
    0.1
}

fn default_self_pol_mu0() -> f64 {
    0.0
}

fn default_self_pol_alpha() -> f64 {
    0.0
}

fn default_epsgrad() -> f64 {
    0.0
}

fn default_energy_asymmetry() -> f64 {
    1.0
}

fn default_liquid_fdiff_h() -> f64 {
    1e-2
}

fn default_restrain_k() -> f64 {
    1.0
}

fn default_force_map() -> String {
    String::from("residue")
}

fn default_inter_txt() -> String {
    String::from("interactions.txt")
}

fn default_optgeo_options_txt() -> String {
    String::from("optgeo_options.txt")
}

fn default_gmx_eq_barostat() -> String {
    String::from("berendsen")
}

fn default_evaluator_input() -> String {
    String::from("evaluator_input.json")
}

fn default_recharge_esp_store() -> String {
    String::from("esp-store.sqlite")
}

fn default_recharge_property() -> String {
    String::from("esp")
}

fn default_energy_mode() -> String {
    String::from("average")
}

fn default_expdata_txt() -> String {
    String::from("expset.txt")
}

fn default_hfedata_txt() -> String {
    String::from("hfedata.txt")
}

fn default_hfemode() -> String {
    String::from("single")
}

fn default_sleepy() -> usize {
    0
}

fn default_liquid_md_steps() -> usize {
    10000
}

fn default_liquid_eq_steps() -> usize {
    1000
}

fn default_lipid_md_steps() -> usize {
    10000
}

fn default_lipid_eq_steps() -> usize {
    1000
}

fn default_n_mcbarostat() -> usize {
    25
}

fn default_gas_md_steps() -> usize {
    100000
}

fn default_gas_eq_steps() -> usize {
    10000
}

fn default_nvt_md_steps() -> usize {
    100000
}

fn default_nvt_eq_steps() -> usize {
    10000
}

fn default_writelevel() -> usize {
    0
}

fn default_md_threads() -> usize {
    1
}

fn default_save_traj() -> usize {
    0
}

fn default_eq_steps() -> usize {
    20000
}

fn default_md_steps() -> usize {
    50000
}

fn default_n_sim_chain() -> usize {
    1
}

fn default_hess_normalize_type() -> usize {
    0
}

fn default_all_at_once() -> bool {
    true
}

fn default_run_internal() -> bool {
    true
}

fn default_energy() -> bool {
    true
}

fn default_force() -> bool {
    true
}

fn default_optimize_geometry() -> bool {
    true
}

fn default_minimize_energy() -> bool {
    true
}
