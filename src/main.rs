use std::error::Error;

use anakin::{
    config::Config, forcefield::FF, objective::Objective, optimizer::Optimizer,
};

fn main() -> Result<(), Box<dyn Error>> {
    let config = match Config::load("testfiles/optimize.in") {
        Ok(c) => c,
        Err(e) => {
            panic!("{e:#?}");
        }
    };
    let forcefield = FF::new(&config);
    // TODO probably this forcefield is supposed to be a reference and not
    // separate clones, but we'll see
    let objective = Objective::new(forcefield.clone());
    let mut optimizer = Optimizer::new(objective, forcefield);

    optimizer.run();

    Ok(())
}
