# Writes llms.txt at the root of the repository from docs/src/llms.md: the same text with
# each documentation link `[`name`](@ref)` reduced to `name`. test/docrefs.jl checks that
# the file matches the page, and the page's links are checked against the exported names
# there, so a name that stops existing fails the guard.
#
#   julia dev/llms.jl

llms_text(md::AbstractString) = replace(md, r"\[(`[^`\n]+`)\]\(@ref[^)\n]*\)" => s"\1")

if abspath(PROGRAM_FILE) == @__FILE__
    root = dirname(@__DIR__)
    write(joinpath(root, "llms.txt"), llms_text(read(joinpath(root, "docs", "src", "llms.md"), String)))
    println("wrote llms.txt")
end
