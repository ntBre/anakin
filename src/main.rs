use anakin::{forcefield::FF, objective::Objective, optimizer::Optimizer};

fn main() {
    let forcefield = FF::new();
    // TODO probably this forcefield is supposed to be a reference and not
    // separate clones, but we'll see
    let objective = Objective::new(forcefield.clone());
    let optimizer = Optimizer::new(objective, forcefield);

    optimizer.run();
}
