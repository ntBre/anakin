use std::error::Error;

use anakin::{
    config::Config, forcefield::FF, objective::Objective, optimizer::Optimizer,
};
use log::debug;

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::init();
    let config = match Config::load("testfiles/optimize.in") {
        Ok(c) => c,
        Err(e) => {
            panic!("{e:#?}");
        }
    };
    debug!("finished loading config");
    let forcefield = FF::new(&config);
    debug!("finished initializing forcefield");
    // TODO probably this forcefield is supposed to be a reference and not
    // separate clones, but we'll see
    let objective = Objective::new(forcefield.clone());
    debug!("finished initializing objective");
    let mut optimizer = Optimizer::new(objective, forcefield);
    debug!("finished initializing optimizer");

    optimizer.run();

    Ok(())
}
