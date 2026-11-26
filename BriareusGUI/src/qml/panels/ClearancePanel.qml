// The clearance report shows the closest beam per component and every beam–component
// row, colored by status, with CSV export. Beside it, the keep-out map is drawn flat: the
// two hemispheres of the sphere about TCC as disks, or the plane as a square.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

RowLayout {
  spacing: 10
  function statusColor(s) {
      if (s === "fail") return "#c42b1c"
      if (s === "review") return "#d08a00"
      return "#2e9e44"
  }
  ColumnLayout {
    spacing: 6
    Layout.fillWidth: true
    Layout.fillHeight: true
    RowLayout {
        visible: !app.rvHasReport
        spacing: 6
        Note {
            Layout.fillWidth: false
            text: app.rvRunning ? "checked " + app.rvDone + " of " + app.rvTotal + " beams" : "No report yet: run Beam clearance"
            opacity: 0.72
        }
        Info { topic: "clearance.results" }
    }
    RowLayout {
        visible: app.rvHasReport
        spacing: 8
        Cell { text: app.rvBrief; font.weight: Font.DemiBold; Layout.fillWidth: true; Layout.maximumWidth: implicitWidth }
        Info { topic: "clearance.results" }
        Item { Layout.fillWidth: true }
        Label { text: "inputs changed since"; visible: app.rvStale; color: "#d08a00" }
        Button { text: "Export CSV…"; enabled: !app.rvRunning; onClicked: win.openClearanceExport()
                 tip: "Write every beam–component row to a CSV file" }
    }
    Repeater {
        model: worstModel
        delegate: Cell {
            text: model.component + ": " + model.beam + " · " +
                  (model.worst === "–" ? "no hardware point before the landing"
                                       : "clearance " + model.worst + " mm (nominal " + model.nominal + " mm)") +
                  " · " + model.power + " TW · " + model.status
            color: statusColor(model.status)
            Layout.fillWidth: true
        }
    }
    // Without a report nothing else takes the height, so the note stays at the top.
    Item { visible: !app.rvHasReport; Layout.fillHeight: true }
    HScroll {
        visible: app.rvHasReport
        Layout.fillWidth: true
        Layout.fillHeight: true
        Layout.minimumHeight: 80
        minimumContentWidth: head.totalWidth
        TableHeader {
            id: head
            columns: [{ t: "Beam", w: 64 }, { t: "Component", w: 110 }, { t: "Status", w: 56 },
                      { t: "Worst (mm)", w: 80 }, { t: "Nominal", w: 72 },
                      { t: "Power (TW)", w: 80 }, { t: "Peak (W/cm²)", w: 96 }]
        }
        ListView {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Layout.minimumHeight: 60
            clip: true
            model: clearModel
            boundsBehavior: Flickable.StopAtBounds
            ScrollBar.vertical: ScrollBar {}
            delegate: RowLayout {
                spacing: 4
                Cell { text: model.beam; Layout.preferredWidth: head.widths[0] }
                Cell { text: model.component; Layout.preferredWidth: head.widths[1] }
                Cell { text: model.status; Layout.preferredWidth: head.widths[2]; color: statusColor(model.status) }
                Cell { text: model.worst; Layout.preferredWidth: head.widths[3] }
                Cell { text: model.nominal; Layout.preferredWidth: head.widths[4]; opacity: 0.7 }
                Cell { text: model.power; Layout.preferredWidth: head.widths[5] }
                Cell { text: model.intensity; Layout.preferredWidth: head.widths[6] }
            }
        }
    }
  }
  // the keep-out map, flat: the colors of Briareus.corridor_map_svg
  ColumnLayout {
    visible: app.cmShown && app.cmMap !== ""
    Layout.preferredWidth: 380
    Layout.minimumWidth: 200
    Layout.fillHeight: true
    spacing: 4
    RowLayout {
        spacing: 8
        Cell { text: "Keep-out map"; font.weight: Font.DemiBold; Layout.fillWidth: true }
        Button { text: "Export SVG…"; onClicked: win.openCorridorExport()
                 tip: "Write the map as a vector drawing, with its legend" }
    }
    Canvas {
        id: mapCanvas
        Layout.fillWidth: true
        Layout.fillHeight: true
        Layout.minimumHeight: 120
        property string drawing: app.cmMap
        property color ink: palette.windowText
        readonly property var fills: ["#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#B7E3D4", "#BBBBBB"]
        onDrawingChanged: requestPaint()
        onInkChanged: requestPaint()
        onWidthChanged: requestPaint()
        onHeightChanged: requestPaint()
        onPaint: {
            var ctx = getContext("2d")
            ctx.reset()
            if (drawing === "") return
            var m = JSON.parse(drawing)
            var nu = m.nu, nv = m.nv, f = m.fills
            function cls(i, j) { return f.charCodeAt(i + nu * j) - 48 }
            ctx.font = "11px sans-serif"
            if (m.surface === "sphere") {
                var size = Math.min(width / 2 - 6, height - 18)
                if (size < 20) return
                var R = size / 2 - 2, dθ = 180 / nu, dφ = 2 * Math.PI / nv
                for (var h = 0; h < 2; h++) {
                    var cx = h * (size + 12) + size / 2, cy = size / 2
                    for (var i = 0; i < nu; i++) {
                        var θ0 = i * dθ, θ1 = (i + 1) * dθ
                        if (h === 1) { var t = θ0; θ0 = 180 - θ1; θ1 = 180 - t }
                        if (θ0 >= 90 - 1e-9) continue
                        var r0 = R * θ0 / 90, r1 = R * Math.min(θ1, 90) / 90
                        var j = 0
                        while (j < nv) {
                            var c = cls(i, j), k = j
                            while (k + 1 < nv && cls(i, k + 1) === c) k++
                            // the run as a polygon, one corner per cell on each arc
                            ctx.beginPath()
                            for (var a = j; a <= k + 1; a++)
                                ctx.lineTo(cx + r1 * Math.cos(a * dφ), cy - r1 * Math.sin(a * dφ))
                            for (var b = k + 1; b >= j; b--)
                                ctx.lineTo(cx + r0 * Math.cos(b * dφ), cy - r0 * Math.sin(b * dφ))
                            ctx.closePath()
                            ctx.fillStyle = fills[c - 1]; ctx.fill()
                            ctx.strokeStyle = fills[c - 1]; ctx.lineWidth = 0.5; ctx.stroke()
                            j = k + 1
                        }
                    }
                    ctx.strokeStyle = ink; ctx.globalAlpha = 0.5; ctx.lineWidth = 0.5
                    for (var g = 1; g <= 3; g++) { ctx.beginPath(); ctx.arc(cx, cy, R * g / 3, 0, 2 * Math.PI); ctx.stroke() }
                    ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx + R, cy); ctx.stroke()
                    ctx.globalAlpha = 1; ctx.fillStyle = ink
                    ctx.fillText(h === 0 ? "θ < 90° (φ = 0 to the right)" : "θ > 90°", cx - R, size + 12)
                }
            } else {
                var side = Math.min(width, height - 18)
                if (side < 20) return
                var w = side / nu
                for (var jj = 0; jj < nv; jj++) {
                    for (var ii = 0; ii < nu; ii++) {
                        ctx.fillStyle = fills[cls(ii, jj) - 1]
                        ctx.fillRect(ii * w, side - (jj + 1) * w, w + 0.5, w + 0.5)
                    }
                }
                ctx.strokeStyle = ink; ctx.lineWidth = 0.5; ctx.strokeRect(0.5, 0.5, side - 1, side - 1)
                ctx.fillStyle = ink
                ctx.fillText("±" + m.extent + " mm about the point", 2, side + 12)
            }
        }
    }
    Flow {
        Layout.fillWidth: true
        spacing: 10
        Repeater {
            model: [["#0072B2", "3ω light"], ["#D55E00", "1ω light"], ["#CC79A7", "line of sight"],
                    ["#F0E442", "open, tight"], ["#B7E3D4", "open"]].concat(
                       app.cmMap !== "" && JSON.parse(app.cmMap).fills.indexOf("6") >= 0 ? [["#BBBBBB", "shadowed by hardware"]] : [])
            delegate: Row {
                spacing: 4
                Rectangle { width: 12; height: 12; color: modelData[0]; border.width: 0.5; border.color: "#000000" }
                Label { text: modelData[1] }
            }
        }
    }
  }
}
