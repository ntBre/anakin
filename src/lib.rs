//! Rust version of the optimization package ForceBalance

pub mod forcefield {
    #[derive(Clone)]
    pub struct FF;

    impl FF {
        pub fn new() -> Self {
            FF
        }
    }
}

pub mod objective {
    use crate::forcefield::FF;

    pub struct Objective {
        forcefield: FF,
    }

    impl Objective {
        pub fn new(forcefield: FF) -> Self {
            Self { forcefield }
        }
    }
}

pub mod optimizer {
    use crate::{forcefield::FF, objective::Objective};

    pub struct Optimizer {
        objective: Objective,
        forcefield: FF,
    }

    impl Optimizer {
        pub fn new(objective: Objective, forcefield: FF) -> Self {
            Self {
                objective,
                forcefield,
            }
        }

        pub fn run(&self) {}
    }
}
