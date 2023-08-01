#!/usr/bin/awk -f

# convert tgt_opt_types from parser.py to Rust types

/bool,/ {
    gsub(/:/, "", $1)
    name = $1
    def = $NF
    if (def != "0") {
	printf "#[serde(default = \"default_%s\")]\n", name
	fns[i++] = sprintf("fn default_%s() -> bool { true }\n", name)
    } else {
	printf "#[serde(default)]\n", name
    }
    printf "%s: %s\n", $1, $2
    next
}

{
    print
}

END {
    for (f in fns) print fns[f]
}
