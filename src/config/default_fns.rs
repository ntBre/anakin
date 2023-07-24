pub(crate) fn default_penalty_type() -> String {
    String::from("L2")
}

pub(crate) fn default_ffdir() -> String {
    String::from("forcefield")
}

pub(crate) fn default_maxstep() -> isize {
    100
}

pub(crate) fn default_objective_history() -> isize {
    2
}

pub(crate) fn default_criteria() -> isize {
    1
}

pub(crate) fn default_zerograd() -> isize {
    -1
}

pub(crate) fn default_backup() -> bool {
    true
}

pub(crate) fn default_writechk_step() -> bool {
    true
}

pub(crate) fn default_retain_micro_outputs() -> bool {
    true
}

pub(crate) fn default_print_gradient() -> bool {
    true
}

pub(crate) fn default_print_parameters() -> bool {
    true
}

pub(crate) fn default_normalize_weights() -> bool {
    true
}

pub(crate) fn default_trust0() -> f64 {
    0.1
}

pub(crate) fn default_convergence_objective() -> f64 {
    0.0001
}

pub(crate) fn default_convergence_gradient() -> f64 {
    0.001
}

pub(crate) fn default_convergence_step() -> f64 {
    0.0001
}

pub(crate) fn default_eig_lowerbound() -> f64 {
    0.0001
}

pub(crate) fn default_step_lowerbound() -> f64 {
    1e-06
}

pub(crate) fn default_lm_guess() -> f64 {
    1.0
}

pub(crate) fn default_finite_difference_h() -> f64 {
    0.001
}

pub(crate) fn default_finite_difference_factor() -> f64 {
    0.1
}

pub(crate) fn default_penalty_alpha() -> f64 {
    0.001
}

pub(crate) fn default_penalty_hyperbolic_b() -> f64 {
    1e-06
}

pub(crate) fn default_penalty_power() -> f64 {
    2.0
}

pub(crate) fn default_adaptive_factor() -> f64 {
    0.25
}

pub(crate) fn default_adaptive_damping() -> f64 {
    0.5
}

pub(crate) fn default_search_tolerance() -> f64 {
    0.0001
}
