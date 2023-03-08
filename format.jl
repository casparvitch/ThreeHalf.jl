using JuliaFormatter

format(
    ["./src/", "./test/", "format.jl"],
    verbose = true,
    format_markdown = true,
    conditional_to_if = true,
    format_docstrings = true,
    always_use_return = true,
    import_to_using = true,
    remove_extra_newlines = true,
    whitespace_ops_in_indices = true,
    whitespace_typedefs = true,
    margin = 80,
    always_for_in = true,
)
