# This file reads CSV files without dependencies: `read_csv`, the fields of the writers and
# the digitized tables of data/.

_csvfield(x::AbstractFloat) = @sprintf("%.10g", x)
_csvfield(x::Symbol) = _csvfield(string(x))
_csvfield(x::AbstractString) = occursin(r"[,\"]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : x
_csvfield(x) = string(x)

"""
    read_csv(path; comment = '#') -> (header::Vector{String}, rows::Vector{Vector{String}})

CSV reader without dependencies (RFC 4180 quoting). A field in double quotes may contain
the delimiter, line breaks and the comment character, with `""` for a quote. Its text
is kept as quoted, while unquoted fields are stripped of surrounding blanks. Lines that
are blank or start with `comment` (after blanks, outside a quoted field) are skipped.
Tab- or semicolon-separated files are detected from the header line. A row with fewer
fields than the header, a row with more (other than empty trailing fields), and a quoted
field left open at the end of the file throw an `ArgumentError` with the line number.
"""
function read_csv(path::AbstractString; comment = '#')
    header, rows, _ = _read_csv(path; comment)
    return header, rows
end

# Reads a digitized table in data/, with numeric columns as Float64 vectors and others
# as strings.
function _data_table(file)
    h, rows = read_csv(joinpath(@__DIR__, "..", "data", file))
    col(i) = (v = [r[i] for r in rows];
              all(x -> tryparse(Float64, x) !== nothing, v) ? parse.(Float64, v) : String.(v))
    return NamedTuple{Tuple(Symbol.(h))}(Tuple(col(i) for i in eachindex(h)))
end

# `read_csv`, which also returns the file line on which each row starts.
function _read_csv(path::AbstractString; comment = '#')
    text = read(path, String)
    startswith(text, '﻿') && (text = text[nextind(text, 1):end])     # byte-order mark
    skipline(l) = (s = lstrip(l); isempty(s) || startswith(s, comment))
    first_line = findfirst(!skipline, split(text, '\n'))
    first_line === nothing && throw(ArgumentError("$path has no rows"))
    l1 = split(text, '\n')[first_line]
    delim = occursin('\t', l1) ? '\t' : (count(==(';'), l1) > count(==(','), l1) ? ';' : ',')
    records, starts = Vector{String}[], Int[]
    fields, buf = String[], IOBuffer()
    inq = quoted = false                 # inside quotes, and whether the field has a quoted part
    line, start, i, n = 1, 1, firstindex(text), lastindex(text)
    atstart = true                       # at the start of a record
    endfield!() = (s = String(take!(buf)); push!(fields, quoted ? s : String(strip(s))); quoted = false)
    function endrecord!()
        endfield!()
        push!(records, fields); push!(starts, start)
        fields = String[]
    end
    while i <= n
        if atstart
            # a blank or comment line before a record is skipped to its end
            j = something(findnext('\n', text, i), n + 1)
            if skipline(text[i:prevind(text, j)])
                i, line = nextind(text, min(j, n)), line + 1
                j > n && break
                continue
            end
            atstart, start = false, line
        end
        c = text[i]
        if inq
            if c == '"'
                if i < n && text[nextind(text, i)] == '"'
                    write(buf, '"'); i = nextind(text, i)
                else
                    inq = false
                end
            elseif !(c == '\r' && i < n && text[nextind(text, i)] == '\n')
                c == '\n' && (line += 1)
                write(buf, c)
            end
        elseif c == '"'
            if !quoted                                 # blanks before the opening quote
                s = String(take!(buf))
                isempty(strip(s)) ? (quoted = true) : write(buf, s)
            end
            inq = true
        elseif c == delim
            endfield!()
        elseif c == '\n' || c == '\r'
            c == '\r' && i < n && text[nextind(text, i)] == '\n' && (i = nextind(text, i))
            endrecord!()
            line += 1
            atstart = true
        elseif !(quoted && isspace(c))                 # blanks after the closing quote
            write(buf, c)
        end
        i = nextind(text, i)
    end
    inq && throw(ArgumentError("$path, line $start: a quoted field is not closed"))
    atstart || endrecord!()
    header = popfirst!(records)
    popfirst!(starts)
    while length(header) > 1 && isempty(header[end])
        pop!(header)
    end
    m = length(header)
    for (r, l) in zip(records, starts)
        while length(r) > m && isempty(r[end])
            pop!(r)
        end
        length(r) == m ||
            throw(ArgumentError("$path, line $l: $(length(r)) fields, the header has $m"))
    end
    return header, records, starts
end

# Parses a number from a CSV field, or throws an error that names the file, line and
# column.
function _csv_number(s::AbstractString, path, line, column)
    x = tryparse(Float64, s)
    x === nothing && throw(ArgumentError("$path, line $line, column $column: \"$s\" is not a number"))
    return x
end
