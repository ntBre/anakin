run:
	LD_LIBRARY_PATH=$(shell conda info --base)/lib RUST_BACKTRACE=1 \
	cargo run
