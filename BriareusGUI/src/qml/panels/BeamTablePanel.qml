// The beam table: pointing, defocus, power and on/off per cone, quad or beam, edited in
// place. Selecting a row selects it in the tree and the 3D view.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    id: page
    spacing: 6
    readonly property var columns: [
        { key: "z", label: "aim z" }, { key: "r", label: "aim r" },
        { key: "phi", label: "tang." }, { key: "defocus", label: "defocus" },
        { key: "power", label: "TW/beam" }]
    readonly property string rowKind: app.beamsBy === "group" ? "group" : app.beamsBy

    RowLayout {
        spacing: 8
        Label { text: "Rows" }
        ComboBox {
            Layout.preferredWidth: 140
            textRole: "text"; valueRole: "value"
            model: app.hasQuads ? [{ text: "per cone", value: "group" }, { text: "per quad", value: "quad" },
                                   { text: "per beam", value: "beam" }]
                                : [{ text: "per cone", value: "group" }, { text: "per beam", value: "beam" }]
            currentIndex: Math.max(0, (app.hasQuads ? ["group", "quad", "beam"] : ["group", "beam"]).indexOf(app.beamsBy))
            onActivated: Julia.br_set_grouping(currentValue)
        }
        Cell { text: app.beamsSummary; opacity: 0.7; Layout.fillWidth: true }
    }
    ErrorText { text: app.beamsError }
    Note { text: app.beamsWarnings; color: "#d08a00"; opacity: 1 }
    // the columns at their widths, scrolled sideways on a narrow window
    HScroll {
        Layout.fillWidth: true
        Layout.fillHeight: true
        Layout.minimumHeight: 100
        minimumContentWidth: 36 + 72 + 24 + 104 + 5 * 72 + 8 * 4
        spacing: 6
        RowLayout {
            spacing: 4
            Label { text: "On"; Layout.preferredWidth: 36 }
            Label { text: "Name"; Layout.preferredWidth: 72 }
            Label { text: "n"; Layout.preferredWidth: 24 }
            Label { text: "θ / φ (°)"; Layout.preferredWidth: 104 }
            Repeater {
                model: page.columns
                Label { text: modelData.label; Layout.preferredWidth: 72; font.pixelSize: 12; opacity: 0.8 }
            }
        }
        ListView {
            id: table
            Layout.fillWidth: true
            Layout.fillHeight: true
            Layout.minimumHeight: 80
            clip: true
            model: beamModel
            boundsBehavior: Flickable.StopAtBounds
            ScrollBar.vertical: ScrollBar {}
            delegate: Rectangle {
                id: beamRow
                width: table.width
                height: 28
                readonly property int rowNumber: model.number
                readonly property string nodeId: (app.beamsBy === "group" ? "group:" : app.beamsBy === "quad" ? "quad:" : "beam:") + model.key
                readonly property var cells: [
                    { key: "z", value: model.z }, { key: "r", value: model.r },
                    { key: "phi", value: model.phiOffset }, { key: "defocus", value: model.defocus },
                    { key: "power", value: model.power }]
                readonly property color hl: palette.highlight
                color: app.selection === nodeId ? Qt.rgba(hl.r, hl.g, hl.b, 0.22) : "transparent"
                TapHandler { onTapped: Julia.br_select(beamRow.nodeId) }
                RowLayout {
                    anchors.fill: parent
                    spacing: 4
                    CheckBox {
                        Layout.preferredWidth: 36
                        checked: model.on !== 0
                        opacity: model.on === 2 ? 0.6 : 1
                        onToggled: Julia.br_set_beam_on(beamRow.rowNumber, checked)
                    }
                    Label { text: model.key; Layout.preferredWidth: 72; elide: Text.ElideRight }
                    Label { text: model.count; Layout.preferredWidth: 24; opacity: 0.7 }
                    Label { text: model.theta + " / " + model.phi; Layout.preferredWidth: 104; opacity: 0.7 }
                    Repeater {
                        model: beamRow.cells
                        TextField {
                            Layout.preferredWidth: 72
                            text: modelData.value
                            selectByMouse: true
                            onEditingFinished: if (text !== modelData.value)
                                                   Julia.br_edit_beam(beamRow.rowNumber, modelData.key, text)
                        }
                    }
                    Item { Layout.fillWidth: true }
                }
            }
        }
    }
}
