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
