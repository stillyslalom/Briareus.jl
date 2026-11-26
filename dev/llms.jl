# Writes llms.txt at the root of the repository from docs/src/llms.md: the same text with
# each documentation link `[`name`](@ref)` reduced to `name`. It also writes the Claude
# skill skills/briareus/SKILL.md, the same text under the header that a skill needs.
# test/docrefs.jl checks that both files match the page, and the page's links are
# checked against the exported names there, so a name that stops existing fails the
# guard.
#
#   julia dev/llms.jl

llms_text(md::AbstractString) = replace(md, r"\[(`[^`\n]+`)\]\(@ref[^)\n]*\)" => s"\1")

# The header of the skill: its name, and the description from which an agent decides to
# load it.
const SKILL_HEADER = """
---
name: briareus
description: Write, run and review Briareus.jl scripts and briareus commands for laser pointing, view-factor x-ray drive, TaLIS and 1ω clearance, beam rules of engagement and VisRad workspaces at NIF and OMEGA. Use for hohlraum, halfraum, backlighter and direct-drive designs, .vrw files and .briar sessions.
---

"""

skill_text(md::AbstractString) = SKILL_HEADER * llms_text(md)

if abspath(PROGRAM_FILE) == @__FILE__
    root = dirname(@__DIR__)
    page = read(joinpath(root, "docs", "src", "llms.md"), String)
    write(joinpath(root, "llms.txt"), llms_text(page))
    mkpath(joinpath(root, "skills", "briareus"))
    write(joinpath(root, "skills", "briareus", "SKILL.md"), skill_text(page))
    println("wrote llms.txt and skills/briareus/SKILL.md")
end
