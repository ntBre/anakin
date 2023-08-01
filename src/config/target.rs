use serde::Deserialize;

/// struct corresponding to tgt_opts in ForceBalance
#[derive(Debug, Deserialize)]
pub(crate) struct Target {
    name: String,

    /// Weight of the target (determines its importance vs. other targets)
    weight: f64, // 1.0

    /// Energy denominator in kcal/mol for objective function and lower energy
    /// limit for attenuating weights where applicable
    energy_denom: f64, // 1.0

    /// Upper energy cutoff in kcal/mol for setting weights to zero, used to
    /// exclude super-repulsive interactions
    energy_upper: f64, // 30.0

    /// Weight of experimental density
    w_rho: f64, // 1.0

    /// Weight of enthalpy of vaporization
    w_hvap: f64, // 1.0

    /// Weight of thermal expansion coefficient
    w_alpha: f64, // 1.0

    /// Weight of isothermal compressibility
    w_kappa: f64, // 1.0

    /// Weight of isobaric heat capacity
    w_cp: f64, // 1.0

    /// Weight of dielectric constant
    w_eps0: f64, // 1.0

    /// Weight of average area per lipid
    w_al: f64, // 1.0

    /// Weight of deuterium order parameter
    w_scd: f64, // 1.0

    /// Weight of energy
    w_energy: f64, // 1.0

    /// Weight of atomistic forces
    w_force: f64, // 1.0

    /// Weight of surface tension
    w_surf_ten: f64, // 0.0

    /// Weight of net forces (condensed to molecules, residues, or charge
    /// groups)
    w_netforce: f64, // 0.0

    /// Weight of torques (condensed to molecules, residues, or charge groups)
    w_torque: f64, // 0.0

    /// Weight of RESP
    w_resp: f64, // 0.0

    /// a
    resp_a: f64, // 0.001

    /// b
    resp_b: f64, // 0.1

    /// Simulation temperature for hydration free energies (Kelvin)
    hfe_temperature: f64, // 298.15

    /// Simulation temperature for hydration free energies (atm)
    hfe_pressure: f64, // 1.0

    /// If nonzero, override the Energy RMS used to normalize the energy part of
    /// the objective function term
    energy_rms_override: f64, // 0.0

    /// If nonzero, override the Force RMS used to normalize the energy part of
    /// the objective function term
    force_rms_override: f64, // 0.0

    /// RMSD normalization for optimized geometries in Angstrom
    rmsd_denom: f64, // 0.1

    /// Frequency normalization (in wavenumber) for vibrational frequencies
    wavenumber_tol: f64, // 10.0

    /// Dipole normalization (Debye) ; set to 0 if a zero weight is desired
    dipole_denom: f64, // 1.0

    /// Quadrupole normalization (Buckingham) ; set to 0 if a zero weight is
    /// desired
    quadrupole_denom: f64, // 1.0

    /// Dipole polarizability tensor normalization (cubic Angstrom) ; set to 0
    /// if a zero weight is desired
    polarizability_denom: f64, // 1.0

    /// Time step size for the liquid simulation.
    liquid_timestep: f64, // 1.0

    /// Time interval for saving coordinates for the liquid production run.
    liquid_interval: f64, // 0.1

    /// Time step size for the gas simulation (if zero, use default in external
    /// script.).
    gas_timestep: f64, // 1.0

    /// Time interval for saving coordinates for the gas production run (if
    /// zero, use default in external script.)
    gas_interval: f64, // 0.1

    /// Time step size for the lipid simulation.
    lipid_timestep: f64, // 1.0

    /// Time interval for saving coordinates for the lipid production run.
    lipid_interval: f64, // 0.1

    /// Time step size for the NVT simulation.
    nvt_timestep: f64, // 1.0

    /// Time interval for saving coordinates for the NVT simulation production
    /// run.
    nvt_interval: f64, // 0.1

    /// Gas-phase dipole parameter for self-polarization correction (in debye).
    self_pol_mu0: f64, // 0.0

    /// Polarizability parameter for self-polarization correction (in debye).
    self_pol_alpha: f64, // 0.0

    /// Gradient below this threshold will be set to zero.
    epsgrad: f64, // 0.0

    /// Snapshots with (E_MM - E_QM) < 0.0 will have their weights increased by
    /// this factor. Only valid if energy_mode is set to "qm_minimum".
    energy_asymmetry: f64, // 1.0

    /// Cutoff for nonbonded interactions (passed to engines).
    nonbonded_cutoff: f64, // None

    /// Cutoff for vdW interactions if different from other nonbonded
    /// interactions
    vdw_cutoff: f64, // None

    /// Step size for finite difference derivatives for liquid targets in
    /// pure_num_grad
    liquid_fdiff_h: f64, // 1e-2

    /// Force constant for harmonic positional energy restraints
    restrain_k: f64, // 1.0
}
