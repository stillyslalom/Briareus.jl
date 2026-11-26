// The Outages tab shows the chamber port map of the last rules-of-engagement scan, each
// quad (or beam) in the color of its tier, and a quad with single beams of another tier
// with that tier in its center. A click on a mark adds it to a hypothetical
// outage, or takes it out, and the table beside the map scores the units picked so far
// together against every criterion of the scan.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 6
    function tierColor(t) {
        if (t === "red") return "#D55E00"
        if (t === "amber") return "#E69F00"
        if (t === "green") return "#009E73"
        return "#BBBBBB"
    }
    RowLayout {
        spacing: 8
        Button { text: "Scan outages"; onClicked: Julia.br_roe_start()
                 tip: "Score the loss of each beam and quad against the default criteria of the design (engineering defaults)" }
        Button { text: "Clear the outage"; enabled: app.roeOutageText !== ""; onClicked: Julia.br_roe_clear()
                 tip: "Take every picked unit out of the hypothetical outage" }
        Cell { text: app.roeText === "" ? "No scan yet" : app.roeText; Layout.fillWidth: true; wrapMode: Text.Wrap }
    }
    RowLayout {
        Layout.fillWidth: true
        Layout.fillHeight: true
        spacing: 10
        visible: app.roePorts !== ""
        ColumnLayout {
            Layout.fillHeight: true
            Layout.preferredWidth: 420
            Layout.minimumWidth: 220
            spacing: 2
            Canvas {
                id: portCanvas
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.minimumHeight: 120
                property string ports: app.roePorts
                property color ink: palette.windowText
                property var marks: ports === "" ? [] : JSON.parse(ports)
                // the disks: the size of each and the centers, shared with the clicks
                readonly property real size: Math.max(0, Math.min(width / 2 - 8, height - 18))
                function center(upper) { return upper ? [size / 2, size / 2] : [size * 1.5 + 16, size / 2] }
                function at(m) { const c = center(m.upper), R = size / 2 - 8; return [c[0] + R * m.x, c[1] - R * m.y] }
                onMarksChanged: requestPaint()
                onInkChanged: requestPaint()
                onWidthChanged: requestPaint()
                onHeightChanged: requestPaint()
                onPaint: {
                    var ctx = getContext("2d")
                    ctx.reset()
                    if (size < 40) return
                    var R = size / 2 - 8
                    ctx.strokeStyle = ink; ctx.lineWidth = 0.5
                    for (var h = 0; h < 2; h++) {
                        var c = center(h === 0)
                        ctx.globalAlpha = 0.5
                        for (var g = 1; g <= 3; g++) { ctx.beginPath(); ctx.ellipse(c[0] - R * g / 3, c[1] - R * g / 3, 2 * R * g / 3, 2 * R * g / 3); ctx.stroke() }
                        ctx.globalAlpha = 1
                        ctx.fillStyle = ink; ctx.font = "11px sans-serif"
                        ctx.fillText(h === 0 ? "upper hemisphere" : "lower hemisphere, seen from above", c[0] - R, size + 12)
                    }
                    var r = Math.max(4, Math.min(9, size / 40))
                    for (var i = 0; i < marks.length; i++) {
                        var m = marks[i], p = at(m)
                        ctx.beginPath(); ctx.ellipse(p[0] - r, p[1] - r, 2 * r, 2 * r)
                        ctx.fillStyle = tierColor(m.tier); ctx.fill()
                        ctx.lineWidth = m.picked ? 3 : 0.5; ctx.strokeStyle = m.picked ? ink : "#000000"; ctx.stroke()
                        if (m.beam !== "none" && m.beam !== m.tier) {
                            ctx.beginPath(); ctx.ellipse(p[0] - 0.55 * r, p[1] - 0.55 * r, 1.1 * r, 1.1 * r)
                            ctx.fillStyle = tierColor(m.beam); ctx.fill()
                        }
                        if (m.badge !== "") {
                            // the badge of its role, at the upper right
                            ctx.font = "bold 8px sans-serif"
                            ctx.lineWidth = 2.5; ctx.strokeStyle = "#FFFFFF"
                            ctx.strokeText(m.badge, p[0] + 0.7 * r, p[1] - 0.7 * r)
                            ctx.fillStyle = ink; ctx.fillText(m.badge, p[0] + 0.7 * r, p[1] - 0.7 * r)
                        }
                    }
                }
                MouseArea {
                    anchors.fill: parent
                    hoverEnabled: true
                    // the mark under the pointer, within twice its radius
                    function hit(px, py) {
                        var best = null, bd = 1e9
                        for (var i = 0; i < portCanvas.marks.length; i++) {
                            var p = portCanvas.at(portCanvas.marks[i])
                            var d = Math.hypot(p[0] - px, p[1] - py)
                            if (d < bd) { bd = d; best = portCanvas.marks[i] }
                        }
                        return bd <= 2 * Math.max(4, Math.min(9, portCanvas.size / 40)) ? best : null
                    }
                    onClicked: (e) => { const m = hit(e.x, e.y); if (m) Julia.br_roe_pick(m.name) }
                    onPositionChanged: (e) => { const m = hit(e.x, e.y)
                        hoverName.text = !m ? "" : (m.beam === "none" || m.beam === m.tier ? m.name + ": " + m.tier
                                         : m.name + ": " + m.tier + ", one of its beams " + m.beam) + ", role " + m.role }
                }
            }
            Label { id: hoverName; text: ""; opacity: 0.8 }
            Label {
                // the key of the role badges
                readonly property var pairs: {
                    var seen = {}, out = []
                    for (var i = 0; i < portCanvas.marks.length; i++) {
                        var m = portCanvas.marks[i]
                        if (m.badge !== "" && !seen[m.badge]) { seen[m.badge] = true; out.push(m.badge + " " + m.role) }
                    }
                    return out.sort()
                }
                visible: pairs.length > 0
                text: "Roles: " + pairs.join(" · ")
                opacity: 0.8; font.pixelSize: 11
            }
            Label { text: "Center of a mark: the loss of one beam of the quad"; opacity: 0.8; font.pixelSize: 11
                    visible: portCanvas.marks.some(m => m.beam !== "none" && m.beam !== m.tier) }
        }
        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 4
            Cell { text: app.roeOutageText === "" ? "Click marks on the map to drop several units together" : app.roeOutageText
                   wrapMode: Text.Wrap; Layout.fillWidth: true }
            HScroll {
                visible: app.roeOutageRows !== "[]" && app.roeOutageRows !== ""
                Layout.fillWidth: true
                Layout.fillHeight: true
                minimumContentWidth: head.totalWidth
                TableHeader {
                    id: head
                    columns: [{ t: "Criterion", w: 150 }, { t: "Score", w: 110 }, { t: "Green ≤", w: 70 },
                              { t: "Amber ≤", w: 70 }, { t: "Tier", w: 60 }]
                }
                ListView {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    clip: true
                    model: app.roeOutageRows === "" ? [] : JSON.parse(app.roeOutageRows)
                    boundsBehavior: Flickable.StopAtBounds
                    ScrollBar.vertical: ScrollBar {}
                    delegate: RowLayout {
                        spacing: 4
                        Cell { text: modelData.criterion; Layout.preferredWidth: head.widths[0]; font.weight: modelData.binding ? Font.DemiBold : Font.Normal }
                        Cell { text: modelData.score; Layout.preferredWidth: head.widths[1] }
                        Cell { text: modelData.green; Layout.preferredWidth: head.widths[2]; opacity: 0.7 }
                        Cell { text: modelData.amber; Layout.preferredWidth: head.widths[3]; opacity: 0.7 }
                        Cell { text: modelData.tier; Layout.preferredWidth: head.widths[4]; color: tierColor(modelData.tier) }
                    }
                }
            }
        }
    }
    Item { visible: app.roePorts === ""; Layout.fillHeight: true }
}
