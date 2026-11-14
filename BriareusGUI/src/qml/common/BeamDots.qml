// A beam set's 2 × 3 matrix of dots, as VisRad's: the top row is whether its 1ω (red),
// 2ω (green) and 3ω (blue) light takes part in the calculations, the bottom row whether
// it is drawn. A filled dot is on, a ring off, a pale dot on for some of the beams.
// A click switches the dot for the node's beams (the whole selection when the node is
// in it). A 1ω or 2ω calculation dot switches once that color has a power ratio (the
// Drive's properties); without one it is shown but not switchable.
import QtQuick
import QtQuick.Controls
import jlqml

Grid {
    id: dots
    property string matrix: "000000"
    property string nodeId: ""
    property int dot: 8
    columns: 3
    rowSpacing: 2
    columnSpacing: 3
    readonly property var colors: ["#d62728", "#2ca02c", "#1f5fd6"]
    readonly property var names: ["1ω", "2ω", "3ω"]
    Repeater {
        model: 6
        delegate: Rectangle {
            id: d
            required property int index
            readonly property int h: index % 3 + 1
            readonly property bool use: index < 3
            readonly property string st: dots.matrix.charAt(index)
            readonly property bool usable: st !== "-" && !(use && h < 3 && !app.ucLive[h - 1])
            width: dots.dot
            height: dots.dot
            radius: dots.dot / 2
            color: st === "1" ? dots.colors[h - 1]
                 : st === "2" ? Qt.lighter(dots.colors[h - 1], 1.7) : "transparent"
            border.color: dots.colors[h - 1]
            border.width: 1.2
            opacity: st === "-" ? 0.15 : usable ? 1 : 0.35
            HoverHandler { id: hover }
            Tip {
                visible: hover.hovered
                text: {
                    const n = dots.names[d.h - 1]
                    if (d.st === "-") return "No " + n + " footprint for these beams"
                    if (d.use && d.h < 3 && !app.ucLive[d.h - 1]) return "Set the " + n + " power ratio in Drive first"
                    return d.use ? n + " light in the calculations" : n + " light drawn"
                }
            }
            MouseArea {
                anchors.fill: parent
                anchors.margins: -1
                enabled: d.usable
                onClicked: Julia.br_beam_dot(dots.nodeId, d.use ? "use" : "draw", d.h, d.st !== "1")
            }
        }
    }
}
