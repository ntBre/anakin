use serde::Deserialize;

/// struct corresponding to tgt_opts in ForceBalance
#[derive(Debug, Deserialize)]
pub(crate) struct Target {
    name: String,

    /// Weight of the target (determines its importance vs. other targets)
    #[serde(default = "default_weight")]
    weight: f64,

    /// Energy denominator in kcal/mol for objective function and lower energy
    /// limit for attenuating weights where applicable
    #[serde(default = "default_energy_denom")]
    energy_denom: f64,

    /// Upper energy cutoff in kcal/mol for setting weights to zero, used to
    /// exclude super-repulsive interactions
    #[serde(default = "default_energy_upper")]
    energy_upper: f64,

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
    epsgrad: f64,

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
    restrain_k: f64,
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
