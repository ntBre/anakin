run:
	RUST_LOG=debug RUST_BACKTRACE=1 cargo run

clippy:
	cargo clippy

test:
	cargo test $(ARGS)
