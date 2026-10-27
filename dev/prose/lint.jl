# Prose lint for Briareus: American spelling everywhere, and the docs register (no "we"/"you",
# no software as the subject of verbs for minds or bodies) in README, docs pages, tutorial
# comments and src docstrings (src comments get the spelling check only). Base Julia only, so
# it runs in a hook without a project environment (`--compile=min`: under a second).
# The mind-verb check is a lexicon of inflected forms, a heuristic: it reads no syntax.
#
#   julia dev/prose/lint.jl [--fix] [paths...]   report every hit (default: all tracked files);
#                                                --fix applies the spelling corrections
#   julia dev/prose/lint.jl --staged             pre-commit: hits on lines added since HEAD;
#                                                spelling and we/you fail, mind verbs warn
#   julia dev/prose/lint.jl --hook               Claude Code PostToolUse: tool JSON on stdin;
#                                                hits in the edited text go to stderr, exit 2
#
# Only prose is read: Markdown outside code, Julia comments and docstrings (outside code
# blocks and inline code), `#` lines of TOML and CSV files. Quoted text is exempt from the
# register checks. ALLOW holds lines that are never flagged.

const FIXED = Dict(
    "centre" => "center", "centres" => "centers", "centred" => "centered",
    "centring" => "centering", "centreline" => "centerline", "centrelines" => "centerlines",
    "colour" => "color", "colours" => "colors", "coloured" => "colored", "colouring" => "coloring",
    "behaviour" => "behavior", "behaviours" => "behaviors", "favour" => "favor",
    "favours" => "favors", "favoured" => "favored", "favourable" => "favorable",
    "neighbour" => "neighbor", "neighbours" => "neighbors", "neighbouring" => "neighboring",
    "flavour" => "flavor", "vapour" => "vapor", "labour" => "labor", "honour" => "honor",
    "rumour" => "rumor", "metre" => "meter", "metres" => "meters", "millimetre" => "millimeter",
    "millimetres" => "millimeters", "micrometre" => "micrometer", "micrometres" => "micrometers",
    "centimetre" => "centimeter", "centimetres" => "centimeters", "kilometre" => "kilometer",
    "kilometres" => "kilometers", "nanometre" => "nanometer", "nanometres" => "nanometers",
    "litre" => "liter", "litres" => "liters", "fibre" => "fiber", "fibres" => "fibers",
    "labelled" => "labeled", "labelling" => "labeling", "modelled" => "modeled",
    "modelling" => "modeling", "cancelled" => "canceled", "cancelling" => "canceling",
    "travelled" => "traveled", "travelling" => "traveling", "signalled" => "signaled",
    "signalling" => "signaling", "totalled" => "totaled", "totalling" => "totaling",
    "fuelled" => "fueled", "levelled" => "leveled", "levelling" => "leveling",
    "tunnelling" => "tunneling", "channelled" => "channeled", "funnelled" => "funneled",
    "analyse" => "analyze", "analysed" => "analyzed", "analysing" => "analyzing",
    "catalyse" => "catalyze", "paralyse" => "paralyze",
    "analogue" => "analog", "analogues" => "analogs", "catalogue" => "catalog",
    "catalogues" => "catalogs", "programme" => "program", "programmes" => "programs",
    "licence" => "license", "defence" => "defense", "grey" => "gray", "greys" => "grays",
    "artefact" => "artifact", "artefacts" => "artifacts", "towards" => "toward",
    "whilst" => "while", "amongst" => "among", "aluminium" => "aluminum", "sulphur" => "sulfur",
    "judgement" => "judgment", "ageing" => "aging", "focussed" => "focused",
    "focussing" => "focusing", "manoeuvre" => "maneuver", "tyre" => "tire", "mould" => "mold",
    "sceptical" => "skeptical", "learnt" => "learned", "spelt" => "spelled", "fulfil" => "fulfill",
)
# -ise/-isation become -ize/-ization, except words spelled -ise in American English too
const ISE = r"^([a-z]{3,}?)is(e|ed|es|ing|ation|ations|er|ers)$"
const ISE_OK = Set(split("""
    rise rises risen rising arise arises noise noises raise raises raised raising premise
    premises imprecise precise concise disguise revise revised revises exercise exercises
    promise promises comprise comprises comprised compromise advise advised surprise surprised
    surprising surprisingly expertise supervise devise devised enterprise franchise improvise
    chastise despise praise poise cruise treatise paradise merchandise demise excise incise
    apprise reprise valise anise mortise televise precisely concisely"""))

# verbs for minds and bodies, third person singular only: base forms are mostly imperatives
# or plurals with human subjects, and past forms mostly participles ("decided per block",
# "the pairs tried"); "agrees", "aims", "assumes" and "knows" are left out: academic usage
# for quantities and models ("the two agree", "the model assumes")
const MIND = Set(split("""
    asks says tells explains claims argues decides thinks wants tries seeks intends expects
    believes prefers swallows cares refuses insists promises remembers forgets learns chooses
    worries likes hates hopes wishes feels complains"""))
# base forms, flagged after a plural noun ("the pieces remember", "the priors expect")
const MIND_BASE = Set(split("""
    ask say tell explain claim argue decide think want try seek intend expect believe prefer
    swallow care refuse insist promise remember forget learn choose worry hate hope wish
    feel complain"""))
# words in -s that are not plural nouns
const NOT_PLURAL = Set(split("""
    this its his is was has does as us thus plus less unless across always perhaps whereas yes
    sometimes ones others axis basis analysis gas lens mass bias radius focus status"""))
# words before a mind verb that make it a person's (or a citation's)
const HUMAN = Set(split("""
    i he she they who user users author authors reader readers one someone people al. al
    physicists experimenters designers"""))
const PERSON = Set(["we", "us", "our", "ours", "ourselves", "you", "your", "yours", "yourself"])

# (file, substring): lines never flagged
const ALLOW = [
    ("reference/research-notes.md", "moves the lens towards the target"),  # a quotation
]
const EXTS = (".jl", ".md", ".toml", ".txt", ".csv")
const REGISTER = r"^(README\.md|docs/src/.*\.md|docs/literate/.*\.jl|src/.*\.jl)$"

# --- prose masks: same-length Char vectors, non-prose replaced by spaces ---

blank!(cs, r) = (cs[r] .= ' '; cs)

"Mask inline code, math, link targets and HTML tags in one Markdown line."
function mask_inline!(cs::Vector{Char})
    s = String(cs)
    for re in (r"``.*?``", r"`[^`]*`", r"\$[^$]*\$", r"\]\([^)]*\)", r"<[^>]*>", r"https?://\S+")
        for m in eachmatch(re, s)
            link = m.match[1] == ']'              # keep the link text's closing bracket
            a = length(s, 1, m.offset) - 1 + link
            blank!(cs, a+1:a+length(m.match)-link)
        end
        s = String(cs)
    end
    cs
end

"""Per line: a Char vector with everything but prose blanked, and whether it is a Julia
comment. Lines outside `only` (all lines if `nothing`) come back empty."""
function prose(rel::String, lines::Vector{String}, only = nothing)
    out = fill(Char[], length(lines))      # never mutated: masks below are built from collect
    comment = falses(length(lines))
    ext = lowercase(splitext(rel)[2])
    fence = false        # inside ``` in Markdown (or in a docstring / comment block)
    doc = false          # inside a Julia triple-quoted string
    # the lines asked for and those that can change the state (the hook runs interpreted,
    # where a pass over every line of a long file costs a second)
    marks = findall(l -> occursin("```", l) || occursin("\"\"\"", l), lines)
    visit = only === nothing ? eachindex(lines) : sort!(collect(union(only, marks)))
    for i in visit
        l = lines[i]
        cs = only === nothing || i in only ? collect(l) : Char[]
        st = lstrip(l)
        if ext == ".md" || ext == ".txt"
            if startswith(st, "```")
                fence = !fence
                out[i] = blank!(cs, 1:length(cs))
            else
                # (the docs use fenced code only: indented lines are list continuations)
                out[i] = fence ? blank!(cs, 1:length(cs)) : mask_inline!(cs)
            end
        elseif ext == ".jl"
            if startswith(st, "\"\"\"") || endswith(rstrip(l), "\"\"\"") && doc
                single = length(st) > 6 && endswith(rstrip(st), "\"\"\"") && startswith(st, "\"\"\"")
                if single && !isempty(cs)
                    out[i] = mask_inline!(cs)
                    k = something(findfirst("\"\"\"", l), 1:3)
                    blank!(cs, 1:length(l, 1, last(k)))
                    blank!(cs, max(1, length(cs) - 2):length(cs))
                elseif !single
                    doc = !doc
                    fence = false
                    out[i] = blank!(cs, 1:length(cs))
                end
            elseif doc
                if startswith(st, "```")
                    fence = !fence
                    out[i] = blank!(cs, 1:length(cs))
                else
                    out[i] = fence || startswith(l, "        ") ? blank!(cs, 1:length(cs)) :
                             mask_inline!(cs)
                end
            elseif occursin(r"^\s*\"[^\"]+\"\s*$", l)      # one-line docstring
                out[i] = mask_inline!(cs)
            else
                out[i] = mask_comment!(cs, l)
                comment[i] = true
            end
        elseif ext == ".toml" || ext == ".csv"
            out[i] = startswith(st, "#") ? mask_inline!(cs) : blank!(cs, 1:length(cs))
        else
            out[i] = blank!(cs, 1:length(cs))
        end
    end
    out, comment
end

"Keep a Julia line's comment (outside string literals, approximately), blank the code."
function mask_comment!(cs::Vector{Char}, l::String)
    quotes = 0
    for (j, c) in enumerate(cs)
        c == '"' && (quotes += 1)
        if c == '#' && iseven(quotes) && get(cs, j + 1, ' ') != '='
            blank!(cs, 1:j)
            return mask_inline!(cs)
        end
    end
    blank!(cs, 1:length(cs))
end

# --- checks ---

struct Hit
    line::Int
    col::Int          # Char index of the word
    kind::Symbol      # :spelling, :person, :mind
    word::String
    fix::String       # American form (spelling only)
end

function american(w::AbstractString)
    any(isuppercase, w[nextind(w, 1):end]) && !all(isuppercase, w) && return nothing  # CamelCase: a name
    lw = lowercase(w)
    a = get(FIXED, lw, nothing)
    if a === nothing
        m = match(ISE, lw)
        (m === nothing || lw in ISE_OK || endswith(lw, "wise")) && return nothing
        a = m[1] * "iz" * m[2]
    end
    all(isuppercase, w) ? uppercase(a) : isuppercase(w[1]) ? uppercasefirst(a) : a
end

"Blank quoted text (straight and curly double quotes) for the register checks."
function unquote(cs::Vector{Char})
    cs = copy(cs)
    s = String(cs)
    for m in eachmatch(r"\"[^\"]*\"|“[^”]*”", s)
        a = length(s, 1, m.offset) - 1
        blank!(cs, a+1:a+length(m.match))
    end
    cs
end

function words(cs::Vector{Char})
    out = Tuple{Int,String}[]          # (start column, word)
    j = 1
    while j <= length(cs)
        if isletter(cs[j])
            k = j
            while k < length(cs) && isletter(cs[k+1])     # "colour's" -> "colour", "s"
                k += 1
            end
            push!(out, (j, String(cs[j:k])))
            j = k + 1
        else
            j += 1
        end
    end
    out
end

"Hits in `text` (file `rel`), on the lines in `only` (all if `nothing`)."
function lint(rel::String, text::String, only = nothing)
    lines = split(text, '\n') .|> String
    masks, comment = prose(rel, lines, only)
    register = occursin(REGISTER, rel)
    allow = [s for (f, s) in ALLOW if f == rel]
    hits = Hit[]
    for i in (only === nothing ? eachindex(lines) : sort!(collect(only)))
        cs = masks[i]
        any(s -> occursin(s, lines[i]), allow) && continue
        for (j, w) in words(cs)
            a = american(w)
            a === nothing || push!(hits, Hit(i, j, :spelling, w, a))
        end
        # comments in src/ are notes for developers; the register is for docs and docstrings
        (register && !(comment[i] && startswith(rel, "src/"))) || continue
        ws = words(unquote(cs))
        for (n, (j, w)) in enumerate(ws)
            lw = lowercase(w)
            if lw in PERSON
                push!(hits, Hit(i, j, :person, w, ""))
            else
                prev = n > 1 ? lowercase(ws[n-1][2]) : ""
                next = n < length(ws) ? lowercase(ws[n+1][2]) : ""
                before = rstrip(String(cs[1:j-1]))
                # "states" is a verb before a complementizer ("the document states that")
                mind = lw in MIND || lw == "states" && next in ("that", "which", "what", "how",
                                                                 "whether", "only")
                # a base form after a plural noun ("the pieces remember")
                mind |= lw in MIND_BASE && length(prev) > 3 && endswith(prev, 's') &&
                        prev ∉ NOT_PLURAL && endswith(before, prev)
                # a person or citation as subject ("Barlow et al.", "(2024)", "Ba24"); a
                # capitalized word begins a sentence or heading (no subject on this line)
                skip = prev in HUMAN || isuppercase(w[1]) || prev in ("the", "a", "an") ||
                       occursin(r"\d{4}[a-z]?\)$|[A-Za-z]\d{2,4}$", before)
                mind && !skip && push!(hits, Hit(i, j, :mind, w, ""))
            end
        end
        # "cannot tell", "can't tell": what a quantity does not separate
        s = lowercase(String(unquote(cs)))
        for m in eachmatch(r"\b(cannot|can't|can not) tell\b", s)
            push!(hits, Hit(i, length(s, 1, m.offset), :mind, m.match, ""))
        end
    end
    hits
end

const LABEL = Dict(:spelling => "British spelling", :person => "we/you",
                   :mind => "verb for a mind or body")

function describe(rel, lines, h::Hit)
    what = h.kind == :spelling ? "\"$(h.word)\" -> \"$(h.fix)\"" : "\"$(h.word)\""
    "$rel:$(h.line): $(LABEL[h.kind]) $what: $(strip(lines[h.line]))"
end

# --- git ---

gitout(dir, args...) = try
    String(readchomp(pipeline(`git -C $dir $args`; stderr = devnull)))
catch
    nothing
end

head_text(root, rel) = something(gitout(root, "show", "HEAD:$rel"), "")
tracked(root) = split(something(gitout(root, "ls-files"), ""), '\n') .|> String

lintable(rel) = any(e -> endswith(lowercase(rel), e), EXTS) && !endswith(rel, "Manifest.toml") &&
                rel != "dev/prose/lint.jl"

"Hits on lines of `new` that are not lines of `old`."
function added_hits(rel, old, new)
    seen = Set(strip.(split(old, '\n')))
    lines = split(new, '\n')
    lint(rel, new, Set(i for (i, l) in enumerate(lines) if strip(l) ∉ seen))
end

# --- modes ---

function report(root, paths, fix)
    rels = isempty(paths) ? tracked(root) :
           [relpath(abspath(p), root) |> x -> replace(x, '\\' => '/') for p in paths]
    n = 0
    for rel in filter(lintable, rels)
        path = joinpath(root, rel)
        isfile(path) || continue
        text = read(path, String)
        hits = lint(rel, text)
        lines = split(text, '\n')
        foreach(h -> println(describe(rel, lines, h)), hits)
        n += length(hits)
        sp = filter(h -> h.kind == :spelling, hits)
        if fix && !isempty(sp)
            ls = [collect(l) for l in lines]
            for h in sort(sp; by = h -> (h.line, -h.col))
                cs = ls[h.line]
                ls[h.line] = vcat(cs[1:h.col-1], collect(h.fix), cs[h.col+length(h.word):end])
            end
            write(path, join(String.(ls), '\n'))
        end
    end
    println(n, " hits", fix ? " (spelling fixed)" : "")
    n
end

function staged(root)
    files = split(something(gitout(root, "diff", "--cached", "--name-only", "--diff-filter=AM"), ""), '\n')
    fail = false
    for rel in filter(r -> !isempty(r) && lintable(r), String.(files))
        new = something(gitout(root, "show", ":$rel"), "")
        lines = split(new, '\n')
        for h in added_hits(rel, head_text(root, rel), new)
            println(stderr, h.kind == :mind ? "warning: " : "error: ", describe(rel, lines, h))
            fail |= h.kind != :mind
        end
    end
    fail && println(stderr, "prose lint: fix the errors above (julia dev/prose/lint.jl --fix " *
                            "<file> corrects spelling); quoted text is exempt")
    fail ? 1 : 0
end

"A JSON string value for `key` in `json` (first occurrence as an object key)."
function json_string(json::String, key::String)
    m = match(Regex("[{,]\\s*\"$key\"\\s*:\\s*\""), json)
    m === nothing && return nothing
    io = IOBuffer()
    i = m.offset + ncodeunits(m.match)
    while i <= ncodeunits(json)
        c = json[i]
        if c == '"'
            return String(take!(io))
        elseif c == '\\'
            e = json[i+1]
            if e == 'u'
                u = parse(UInt16, json[i+2:i+5]; base = 16)
                i += 6
                if 0xd800 <= u < 0xdc00 && json[i] == '\\'   # surrogate pair
                    lo = parse(UInt16, json[i+2:i+5]; base = 16)
                    print(io, Char(0x10000 + (UInt32(u - 0xd800) << 10) + (lo - 0xdc00)))
                    i += 6
                else
                    print(io, Char(u))
                end
                continue
            end
            print(io, get(Dict('n' => '\n', 't' => '\t', 'r' => '\r', 'b' => '\b', 'f' => '\f'),
                          e, e))
            i += 2
            continue
        end
        print(io, c)
        i = nextind(json, i)
    end
    nothing
end

function hook()
    json = read(stdin, String)
    file = json_string(json, "file_path")
    file === nothing && return 0
    file = abspath(file)
    isfile(file) || return 0
    root = gitout(dirname(file), "rev-parse", "--show-toplevel")
    (root === nothing || !isfile(joinpath(root, "dev", "prose", "lint.jl"))) && return 0
    rel = replace(relpath(file, root), '\\' => '/')
    lintable(rel) || return 0
    new = read(file, String)
    lines = split(new, '\n')
    edit = json_string(json, "tool_name") == "Edit" && json_string(json, "new_string") !== nothing
    hits = if edit                                                # Edit: the edited lines
        ns, os = json_string(json, "new_string"), something(json_string(json, "old_string"), "")
        edited = Set{Int}()
        for r in findall(ns, new)
            l0 = length(split(SubString(new, 1, first(r)), '\n'))  # count() runs interpreted
            union!(edited, l0:l0+count(==('\n'), ns))
        end
        # only what this edit added: hits whose word occurs more often in new than old text
        count_in(t, w) = count(Regex("\\b" * w * "\\b", "i"), t)
        filter(h -> count_in(ns, h.word) > count_in(os, h.word), lint(rel, new, edited))
    else                                              # Write, MultiEdit: lines new since HEAD
        added_hits(rel, head_text(root, rel), new)
    end
    isempty(hits) && return 0
    msg = ["Prose lint (dev/prose/lint.jl) on text just written:";
           [describe(rel, lines, h) for h in hits];
           "Docs register: neutral academic prose; no we/you (use the passive or name the " *
           "thing); software takes literal verbs (computes, checks, selects), not verbs for " *
           "minds or bodies (asks, says, decides, wants). Mind-verb hits are heuristic: leave " *
           "literal uses. American spelling throughout."]
    println(stderr, join(msg, '\n'))
    2
end

function main(args)
    root = something(gitout(pwd(), "rev-parse", "--show-toplevel"), pwd())
    # a bug here must not get in the way of editing: report it and pass
    "--hook" in args && return try
        hook()
    catch e
        println(stderr, "dev/prose/lint.jl --hook failed: ", sprint(showerror, e))
        0
    end
    "--staged" in args && return staged(root)
    report(root, filter(a -> !startswith(a, "--"), args), "--fix" in args)
    0
end

exit(main(ARGS))
