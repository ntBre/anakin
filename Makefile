run:
	LD_LIBRARY_PATH=/home/brent/omsf/clone/openmm/build \
	RUST_LOG=debug RUST_BACKTRACE=1 cargo run

clippy:
	cargo clippy

test:
	LD_LIBRARY_PATH=/home/brent/omsf/clone/openmm/build \
	cargo test $(ARGS)
