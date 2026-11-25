# A small writer for self-contained HTML pages: one file with its style and a few lines
# of script for sorting tables, with no reference to any other file or to the network.
# The review packets (review_packets.jl) are built from it.

# A fragment of HTML that is written as it stands. Every other cell or text is escaped.
struct _HTML
    s::String
end

_html_text(x::_HTML) = x.s
_html_text(x) = _xml_escape(x)

const _HTML_STYLE = """
body{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;font-size:14px;line-height:1.45;color:#1a1a1a;margin:24px auto;max-width:1180px;padding:0 16px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:28px 0 8px;border-bottom:1px solid #ccc;padding-bottom:3px}h3{font-size:14px;margin:18px 0 6px}
p.lead{color:#444;margin:4px 0 12px}p.note{color:#555;font-size:13px}
table{border-collapse:collapse;margin:6px 0 12px;font-size:13px}
th,td{border:1px solid #d0d0d0;padding:3px 8px;text-align:left;vertical-align:top}
th{background:#f0f0f0}table.sortable th{cursor:pointer}td.num{text-align:right;font-variant-numeric:tabular-nums}
.chip{display:inline-block;padding:0 7px;border-radius:9px;color:#fff;font-size:12px;white-space:nowrap}
.pass,.green,.clear{background:#009E73}.review,.amber,.tight{background:#E69F00;color:#1a1a1a}.fail,.red{background:#D55E00}
.curved_safe,.specular_clear{background:#0072B2}.unchecked,.none{background:#888}
figure{margin:10px 0 18px;padding:10px;border:1px solid #d0d0d0;display:flex;flex-wrap:wrap;gap:12px;align-items:flex-start}
figcaption{font-weight:600;flex-basis:100%}figure table{max-width:560px;margin:0}div.scroll{overflow-x:auto}
details{margin:6px 0}summary{cursor:pointer;color:#333}
"""

# Sorts the body of a `table.sortable` by the clicked column: by the `data-v` number of
# each cell where given, by its text otherwise, and in the opposite order on a second click.
const _HTML_SCRIPT = """
document.querySelectorAll('table.sortable th').forEach(function (th) {
  th.addEventListener('click', function () {
    var body = th.closest('table').tBodies[0], i = th.cellIndex, d = th.dataset.d === '1' ? -1 : 1;
    th.dataset.d = d === 1 ? '1' : '0';
    Array.from(body.rows).sort(function (x, y) {
      var a = x.cells[i], b = y.cells[i];
      if (a.dataset.v !== undefined) { return d * (parseFloat(a.dataset.v) - parseFloat(b.dataset.v)); }
      return d * a.textContent.localeCompare(b.textContent);
    }).forEach(function (r) { body.appendChild(r); });
  });
});
"""

# A number cell: right-aligned text with the value as the sort key.
_html_num(x::Real, text::AbstractString) =
    _HTML(string("<td class=\"num\" data-v=\"", isfinite(x) ? @sprintf("%.6g", x) : (x > 0 ? "1e300" : "-1e300"), "\">",
                 _xml_escape(text), "</td>"))
_html_num(x::Real; digits::Integer = 3) = _html_num(x, isfinite(x) ? string(round(x; digits)) : string(x))

# A status or tier as a colored chip, with a sort rank.
function _html_chip(class, text = class; rank::Integer = 0)
    return _HTML(string("<td data-v=\"", rank, "\"><span class=\"chip ", _xml_escape(class), "\">", _xml_escape(text),
                        "</span></td>"))
end

# A table. A cell that is an `_HTML` starting with `<td` is a whole cell, and any other
# value is the text of a cell.
function _html_table(io::IO, head, rows; sortable::Bool = true, id::AbstractString = "")
    print(io, "<div class=\"scroll\"><table", isempty(id) ? "" : string(" id=\"", _xml_escape(id), "\""),
          sortable ? " class=\"sortable\"" : "", "><thead><tr>")
    for h in head
        print(io, "<th>", _html_text(h), "</th>")
    end
    println(io, "</tr></thead><tbody>")
    for row in rows
        print(io, "<tr>")
        for c in row
            if c isa _HTML && startswith(c.s, "<td")
                print(io, c.s)
            else
                print(io, "<td>", _html_text(c), "</td>")
            end
        end
        println(io, "</tr>")
    end
    println(io, "</tbody></table></div>")
    return io
end

# The whole page around `body` (HTML). `generated` is a line of text for the foot of the
# page, such as a date: nothing is written when it is `nothing`, so the page of the same
# input is the same text.
function _html_document(title::AbstractString, body::AbstractString; lead::AbstractString = "",
                        generated = nothing)
    io = IOBuffer()
    println(io, "<!DOCTYPE html>")
    println(io, "<html lang=\"en\"><head><meta charset=\"utf-8\"/>",
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>")
    println(io, "<title>", _xml_escape(title), "</title>")
    println(io, "<style>", _HTML_STYLE, "</style></head><body>")
    println(io, "<h1>", _xml_escape(title), "</h1>")
    isempty(lead) || println(io, "<p class=\"lead\">", _xml_escape(lead), "</p>")
    print(io, body)
    generated === nothing || println(io, "<p class=\"note\">", _xml_escape(generated), "</p>")
    println(io, "<script>", _HTML_SCRIPT, "</script>")
    println(io, "</body></html>")
    return String(take!(io))
end

# An SVG drawing for the inside of a page: as it stands, without a trailing newline.
_html_svg(svg::AbstractString) = _HTML(String(rstrip(svg)))
