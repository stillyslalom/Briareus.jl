// The TaLIS report shows the notes (short lines, wrapped to the panel), the per-rule
// summary and every check colored by status, with CSV export.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 6
    function statusColor(s) {
        if (s === "fail") return "#c42b1c"
        if (s === "review") return "#d08a00"
        return "#2e9e44"
    }
    RowLayout {
        visible: !app.tlHasReport
        spacing: 6
        Note {
            Layout.fillWidth: false
            text: app.tlRunning ? "checking… " + app.tlElapsed : "No report yet: run TaLIS rules"
        }
        Info { topic: "talis.intro" }
    }
    RowLayout {
        visible: app.tlHasReport
        spacing: 8
        Cell { text: app.tlBrief; font.weight: Font.DemiBold; Layout.fillWidth: true; Layout.maximumWidth: implicitWidth }
        Info { topic: "talis.report" }
        Item { Layout.fillWidth: true }
        Label { text: "inputs changed since"; visible: app.tlStale; color: "#d08a00" }
        Button { text: "Pinch"; checkable: true; checked: app.tlPinch; onToggled: Julia.br_tl_pinch(checked)
                 tip: "With a 1ω row selected, also draw the nearest side of every other beam whose 1ω margin to the same part is below the tight threshold" }
        Button { text: "Daylight view"; enabled: app.tlFocusApproach; onClicked: Julia.br_tl_view("daylight")
                 tip: "Look down the 1ω beamlet of the selected row, in orthographic projection: the gap seen between the hardware and the contour is the transverse distance of the rule" }
        Button { text: "Grazing view"; enabled: app.tlFocusApproach; onClicked: Julia.br_tl_view("grazing")
                 tip: "Look along the nearest side of the 1ω beamlet of the selected row, with the beam running up the view" }
        Button { text: "Export CSV…"; enabled: !app.tlRunning; onClicked: win.openTalisExport()
                 tip: "Write every check to a CSV file" }
    }
    // notes and the rule summary beside every check, split by a draggable handle
    SplitView {
        Layout.fillWidth: true
        Layout.fillHeight: true
        orientation: Qt.Horizontal
        ScrollView {
            id: notesView
            SplitView.preferredWidth: 330
            SplitView.minimumWidth: 200
            SplitView.fillWidth: !app.tlHasReport
            contentWidth: availableWidth
            contentHeight: notesColumn.implicitHeight
            clip: true
            ColumnLayout {
                id: notesColumn
                width: notesView.availableWidth - 12          // clear of the split handle
                spacing: 4
                // the tolerance note first, then the rule summary, then the other notes
                Repeater {
                    model: noteModel
                    delegate: RowLayout {
                        visible: model.kind === "tolerance"
                        Layout.fillWidth: true
                        spacing: 6
                        // (the ⓘ follows the text; a long note wraps to the panel)
                        Label { text: model.text; opacity: 0.72; wrapMode: Text.WordWrap
                                Layout.fillWidth: true; Layout.maximumWidth: implicitWidth }
                        Info { topic: "talis.tolerance"; Layout.alignment: Qt.AlignTop }
                        Item { Layout.fillWidth: true }
                    }
                }
                SectionLabel { text: "Rules"; visible: app.tlHasReport }
                // One rule takes two lines, its counts and then its worst check.
                Repeater {
                    model: ruleModel
                    delegate: ColumnLayout {
                        spacing: 0
                        Layout.fillWidth: true
                        Tip {
                            visible: ruleHover.hovered
                            text: model.requirement + "\nworst " + model.worst + " · " + model.where_
                        }
                        HoverHandler { id: ruleHover }
                        Label {
                            text: model.label + "  <font color='" + palette.windowText + "'>" +
                                  (model.fail > 0 ? "<font color='#c42b1c'>" + model.fail + " fail</font> · " : "") +
                                  (model.review > 0 ? "<font color='#d08a00'>" + model.review + " review</font> · " : "") +
                                  model.pass + " pass</font>"
                            textFormat: Text.StyledText
                            color: statusColor(model.status)
                            font.weight: Font.DemiBold
                            wrapMode: Text.WordWrap
                            Layout.fillWidth: true
                        }
                        Label {
                            text: "worst " + model.worst + " · " + model.where_
                            opacity: 0.7
                            wrapMode: Text.WordWrap
                            Layout.fillWidth: true
                        }
                    }
                }
                SectionLabel { text: "Notes"; visible: app.tlHasReport }
                Repeater {
                    model: noteModel
                    delegate: RowLayout {
                        visible: model.kind !== "tolerance"
                        Layout.fillWidth: true
                        spacing: 6
                        Label {
                            text: "· " + model.text
                            opacity: model.kind === "scope" || model.kind === "source" ? 0.7 : 1
                            wrapMode: Text.WordWrap
                            Layout.fillWidth: true
                            Layout.maximumWidth: implicitWidth
                        }
                        Info { topic: "talis.scope"; visible: model.kind === "scope"; Layout.alignment: Qt.AlignTop }
                        Item { Layout.fillWidth: true }
                    }
                }
            }
        }
        // every check
        ColumnLayout {
            visible: app.tlHasReport
            SplitView.fillWidth: app.tlHasReport
            SplitView.minimumWidth: 300
            spacing: 4
            HScroll {
                Layout.fillWidth: true
                Layout.fillHeight: true
                minimumContentWidth: head.totalWidth + 10
                TableHeader {
                    id: head
                    Layout.leftMargin: 10
                    columns: [{ t: "Beam", w: 60 }, { t: "Rule", w: 110 }, { t: "Component", w: 96 },
                              { t: "Status", w: 50 }, { t: "Margin", w: 80 }, { t: "Worst", w: 76 },
                              { t: "Required", w: 70 }]
                }
                ListView {
                    Layout.leftMargin: 10
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    clip: true
                    model: talisModel
                    boundsBehavior: Flickable.StopAtBounds
                    ScrollBar.vertical: ScrollBar {}
                    delegate: RowLayout {
                        spacing: 4
                        opacity: app.tlFocus < 0 || app.tlFocus === index ? 1.0 : 0.55
                        Tip { visible: rowHover.hovered && (model.note !== "" || (app.tlFocus === index && app.tlFocusText !== ""))
                              text: app.tlFocus === index && app.tlFocusText !== "" ? app.tlFocusText : model.note }
                        HoverHandler { id: rowHover }
                        TapHandler { onTapped: Julia.br_tl_focus(index) }
                        Cell { text: model.beam; Layout.preferredWidth: head.widths[0] }
                        Cell { text: model.rule; Layout.preferredWidth: head.widths[1] }
                        Cell { text: model.component; Layout.preferredWidth: head.widths[2] }
                        Cell { text: model.status; Layout.preferredWidth: head.widths[3]; color: statusColor(model.status) }
                        Cell { text: model.margin; Layout.preferredWidth: head.widths[4] }
                        Cell { text: model.worst; Layout.preferredWidth: head.widths[5]; opacity: 0.7 }
                        Cell { text: model.required; Layout.preferredWidth: head.widths[6]; opacity: 0.7 }
                    }
                }
            }
        }
        // the plane across the 1ω beamlet of the selected row, through the hardware point:
        // the 10% contour, the worst-case square, the hardware and the caliper
        ColumnLayout {
            visible: app.tlHasReport && app.tlSlice !== ""
            SplitView.preferredWidth: 280
            SplitView.minimumWidth: 160
            spacing: 2
            Label {
                text: app.tlSliceText
                wrapMode: Text.WordWrap
                Layout.fillWidth: true
                Layout.leftMargin: 8
            }
            Canvas {
                id: sliceCanvas
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.leftMargin: 8
                property string slice: app.tlSlice
                property color ink: palette.windowText
                onSliceChanged: requestPaint()
                onInkChanged: requestPaint()
                onWidthChanged: requestPaint()
                onHeightChanged: requestPaint()
                onPaint: {
                    var ctx = getContext("2d")
                    ctx.reset()
                    if (slice === "") return
                    var s = JSON.parse(slice)
                    var side = Math.min(width, height)
                    if (side < 20) return
                    var p0 = s.window[0], q1 = s.window[3]
                    var sc = side / (s.window[1] - p0)
                    function X(p) { return (p - p0) * sc }
                    function Y(q) { return (q1 - q) * sc }
                    function square(sq) {
                        ctx.beginPath()
                        ctx.rect(X(sq[0] - sq[2]), Y(sq[1] + sq[2]), 2 * sq[2] * sc, 2 * sq[2] * sc)
                    }
                    ctx.save()
                    ctx.beginPath(); ctx.rect(0, 0, side, side); ctx.clip()
                    // the worst-case square, dashed over a tint
                    square(s.worst)
                    ctx.fillStyle = Qt.rgba(0.84, 0.37, 0.0, 0.12); ctx.fill()
                    ctx.strokeStyle = "#D55E00"; ctx.lineWidth = 1
                    if (ctx.setLineDash) ctx.setLineDash([6, 4])
                    ctx.stroke()
                    if (ctx.setLineDash) ctx.setLineDash([])
                    // the 10% contour
                    square(s.nominal)
                    ctx.lineWidth = 2; ctx.stroke()
                    // the hardware
                    ctx.strokeStyle = ink; ctx.lineWidth = 1.5
                    ctx.beginPath()
                    for (var i = 0; i + 3 < s.hardware.length; i += 4) {
                        ctx.moveTo(X(s.hardware[i]), Y(s.hardware[i + 1]))
                        ctx.lineTo(X(s.hardware[i + 2]), Y(s.hardware[i + 3]))
                    }
                    ctx.stroke()
                    // the caliper and its length
                    var c = s.caliper
                    ctx.strokeStyle = "#0072B2"; ctx.fillStyle = "#0072B2"
                    ctx.beginPath(); ctx.moveTo(X(c[0]), Y(c[1])); ctx.lineTo(X(c[2]), Y(c[3])); ctx.stroke()
                    for (var j = 0; j < 4; j += 2) {
                        ctx.beginPath(); ctx.arc(X(c[j]), Y(c[j + 1]), 3, 0, 2 * Math.PI); ctx.fill()
                    }
                    ctx.font = "12px sans-serif"
                    ctx.fillText(s.distance.toFixed(3) + " mm", X((c[0] + c[2]) / 2) + 6, Y((c[1] + c[3]) / 2) - 6)
                    ctx.restore()
                    // the scale bar: 1, 2 or 5 times a power of ten within a quarter of the view
                    var span = (s.window[1] - p0) / 4
                    var bar = Math.pow(10, Math.floor(Math.log(span) / Math.LN10))
                    bar = 5 * bar <= span ? 5 * bar : 2 * bar <= span ? 2 * bar : bar
                    ctx.strokeStyle = ink; ctx.fillStyle = ink; ctx.lineWidth = 2
                    ctx.beginPath(); ctx.moveTo(4, side - 8); ctx.lineTo(4 + bar * sc, side - 8); ctx.stroke()
                    ctx.fillText(bar + " mm", 4, side - 14)
                    ctx.strokeStyle = ink; ctx.lineWidth = 0.5
                    ctx.strokeRect(0.5, 0.5, side - 1, side - 1)
                }
            }
        }
    }
}
