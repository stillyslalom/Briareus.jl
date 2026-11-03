# The CSV reader (RFC 4180 quoting, comments, row checks) and the beam table's errors.

@testset "CSV reader" begin
    csv(text) = (f = tempname() * ".csv"; write(f, text); f)
    # quoted fields hold the delimiter, line breaks, the comment character and "" quotes;
    # their text is kept as quoted (blanks inside), unquoted fields are stripped
    f = csv("""
        # a comment line, then a blank one

        name, note ,value
        a,"x, y",1
        "b","two
        lines",2
          # an indented comment
        "c","# not a comment",3
        d," pad ""q"" ",4
        """)
    h, rows = read_csv(f)
    @test h == ["name", "note", "value"]
    @test rows == [["a", "x, y", "1"], ["b", "two\nlines", "2"], ["c", "# not a comment", "3"],
                   ["d", " pad \"q\" ", "4"]]
    # line numbers of the rows (a multi-line field counts its lines)
    @test Briareus._read_csv(f)[3] == [4, 5, 8, 9]
    # CRLF line ends, a byte-order mark, a tab-separated file
    @test read_csv(csv("﻿a,b\r\n1,\"2\r\n3\"\r\n")) == (["a", "b"], [["1", "2\n3"]])
    @test read_csv(csv("a\tb\n1\t\n")) == (["a", "b"], [["1", ""]])
    # empty trailing fields beyond the header are dropped; a short row or a long one is
    # refused with its line number, and so is an unclosed quote
    @test read_csv(csv("a,b\n1,2,,\n")) == (["a", "b"], [["1", "2"]])
    err(f) = try read_csv(f); "" catch e; e isa ArgumentError ? e.msg : rethrow() end
    @test occursin("line 3: 1 fields, the header has 2", err(csv("a,b\n1,2\n3\n4,5\n")))
    @test occursin("line 2: 3 fields", err(csv("a,b\n1,2,3\n")))
    @test occursin("line 2: a quoted field is not closed", err(csv("a,b\n1,\"2\n3,4\n")))
    @test_throws ArgumentError read_csv(csv("# only a comment\n\n"))
end

@testset "beam table errors name the row and column" begin
    q = nif_quads()[1:3]
    path = write_beam_table(tempname() * ".csv", q)
    lines = readlines(path)
    cols = split(lines[1], ',')
    k = findfirst(==("aim_z"), cols)
    row = split(lines[3], ',')
    row[k] = "abc"
    lines[3] = join(row, ',')
    write(path, join(lines, '\n'))
    msg = try read_beam_table(path); "" catch e; e isa ArgumentError ? e.msg : rethrow() end
    @test occursin("line 3, column aim_z: \"abc\" is not a number", msg)
    write(path, lines[1] * "\n" * lines[2] * "\n" * join(split(lines[3], ',')[1:5], ',') * "\n")
    msg = try read_beam_table(path); "" catch e; e isa ArgumentError ? e.msg : rethrow() end
    @test occursin("line 3: 5 fields", msg)
end
