// A table's column headings, one per `columns` entry ({ t: heading, w: least width }).
// Each column is as wide as its least width or its heading, whichever is wider; the
// rows below take `widths` (and `totalWidth`) so that cells line up with the headings.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

RowLayout {
    id: header
    property var columns: []
    readonly property var widths: {
        const ws = []
        const n = heads.count                     // re-evaluated as the headings appear
        for (let i = 0; i < columns.length; ++i) {
            const it = i < n ? heads.itemAt(i) : null
            ws.push(Math.max(columns[i].w, it ? Math.ceil(it.implicitWidth) : 0))
        }
        return ws
    }
    readonly property real totalWidth: widths.reduce((a, b) => a + b, 0) +
                                       spacing * Math.max(0, widths.length - 1)
    spacing: 4
    Repeater {
        id: heads
        model: header.columns
        Label {
            text: modelData.t
            font.weight: Font.DemiBold
            Layout.preferredWidth: header.widths[index] || modelData.w
        }
    }
}
