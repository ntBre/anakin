from parser import gen_opts_types

default_funcs = []


def build_strings():
    for name, value in gen_opts_types["strings"].items():
        (default, _priority, docs, *rest) = value
        if default:
            default_fn = f"default_{name}"
            func = f"""fn {default_fn}() -> String {{
        String::from(\"{default}\")
            }}
        """
            default_funcs.append(func)
            print(
                f"""/// {docs}
        #[serde(default = \"{default_fn}\")]
        {name}: String,
        """
            )
        else:
            print(
                f"""/// {docs}
        {name}: String,
        """
            )


def build_ints():
    for name, value in gen_opts_types["ints"].items():
        (default, _priority, docs, *rest) = value
        if default:
            default_fn = f"default_{name}"
            func = f"""fn {default_fn}() -> isize {{
        {default}
            }}
        """
            default_funcs.append(func)
            print(
                f"""/// {docs}
        #[serde(default = \"{default_fn}\")]
        {name}: isize,
        """
            )
        else:
            print(
                f"""/// {docs}
        {name}: isize,
        """
            )


def build_bools():
    for name, value in gen_opts_types["bools"].items():
        (default, _priority, docs, *rest) = value
        if default:
            default = bool(default)  # cast from ints
            default_fn = f"default_{name}"
            func = f"""fn {default_fn}() -> bool {{
        {str(default).lower()}
            }}
        """
            default_funcs.append(func)
            print(
                f"""/// {docs}
        #[serde(default = \"{default_fn}\")]
        {name}: bool,
        """
            )
        else:
            print(
                f"""/// {docs}
        {name}: bool,
        """
            )


build_strings()
build_ints()
build_bools()

for func in default_funcs:
    print(func)
