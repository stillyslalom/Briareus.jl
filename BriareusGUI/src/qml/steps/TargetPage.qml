import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

StepPage {
    title: "Target"
    guidance: "Choose a hohlraum or a backlighter foil and its dimensions (mm), and cut " +
              "diagnostic windows into a hohlraum wall. The viewer colors the mesh by " +
              "component; drag to rotate it. The coarse mesh keeps the view factors on the " +
              "Drive step to a few seconds; use the standard mesh for design numbers."

    GridLayout {
        columns: 2
        columnSpacing: 16
        rowSpacing: 8
        Label { text: "Target" }
        ComboBox {
            Layout.preferredWidth: 260
            textRole: "text"; valueRole: "value"
            model: [{ text: "Hohlraum with capsule", value: "hohlraum" },
                    { text: "Backlighter foil", value: "foil" }]
            currentIndex: app.targetKind === "foil" ? 1 : 0
            onActivated: Julia.br_set_target_kind(currentValue)
        }
        Label { text: "Mesh" }
        ComboBox {
            Layout.preferredWidth: 260
            textRole: "text"; valueRole: "value"
            model: [{ text: "Coarse (fast view factors)", value: "coarse" },
                    { text: "Standard", value: "standard" },
                    { text: "Fine", value: "fine" }]
            currentIndex: Math.max(0, ["coarse", "standard", "fine"].indexOf(app.targetResolution))
            onActivated: Julia.br_set_resolution(currentValue)
        }
    }

    Label { text: "Dimensions"; font.weight: Font.DemiBold; Layout.topMargin: 8 }
    GridLayout {
        columns: 3
        columnSpacing: 12
        rowSpacing: 6
        Repeater {
            model: app.targetKind === "foil"
                ? [{ key: "distance", label: "Distance from TCC", unit: "mm" },
                   { key: "width", label: "Foil width", unit: "mm" },
                   { key: "conversion", label: "Laser → x-ray conversion", unit: "fraction" },
                   { key: "main_radius", label: "Main-target sphere radius", unit: "mm (0: none)" }]
                : [{ key: "length", label: "Length", unit: "mm" },
                   { key: "radius", label: "Radius", unit: "mm" },
                   { key: "leh_radius", label: "LEH radius", unit: "mm" },
                   { key: "capsule_radius", label: "Capsule radius", unit: "mm (0: none)" },
                   { key: "wall_albedo", label: "Wall albedo", unit: "fraction" },
                   { key: "wall_conversion", label: "Wall laser → x-ray conversion", unit: "fraction" }]
            delegate: RowLayout {
                Layout.columnSpan: 3
                spacing: 12
                Label { text: modelData.label; Layout.preferredWidth: 220 }
                TextField {
                    Layout.preferredWidth: 110
                    text: app["t_" + modelData.key]
                    onEditingFinished: if (text !== app["t_" + modelData.key])
                                           Julia.br_edit_target(modelData.key, text)
                }
                Label { text: modelData.unit; opacity: 0.75 }
            }
        }
    }
    // diagnostic windows in the hohlraum wall
    ColumnLayout {
        id: windows
        visible: app.targetKind === "hohlraum"
        spacing: 4
        Layout.topMargin: 8
        Layout.fillWidth: true
        readonly property var fields: [
            { key: "phi", label: "φ (°)" }, { key: "z", label: "z (mm)" },
            { key: "width", label: "width" }, { key: "height", label: "height" },
            { key: "albedo", label: "albedo" }, { key: "conversion", label: "conv." }]
        RowLayout {
            spacing: 12
            Label { text: "Diagnostic windows"; font.weight: Font.DemiBold }
            Button {
                text: "Add window"
                onClicked: Julia.br_add_window()
                ToolTip.visible: hovered
                ToolTip.text: "A closed 1 × 1 mm window at the midplane, where it fits"
            }
        }
        Label {
            text: "Rectangular windows cut into the wall: centered at azimuth φ and height z, " +
                  "width in mm of arc and height in mm. A closed window re-emits with its own " +
                  "albedo and laser → x-ray conversion; an open one is a hole that laser light " +
                  "and x-rays leave through. Each window is its own component (window_1, …), " +
                  "outlined in the viewer."
            visible: app.windowCount > 0
            opacity: 0.75
            wrapMode: Text.WordWrap
            Layout.fillWidth: true
        }
        RowLayout {
            visible: app.windowCount > 0
            spacing: 4
            Label { text: ""; Layout.preferredWidth: 64 }
            Repeater {
                model: windows.fields
                Label { text: modelData.label; Layout.preferredWidth: 52; font.pixelSize: 12 }
            }
            Label { text: "open"; font.pixelSize: 12 }
        }
        Repeater {
            model: windowModel
            delegate: RowLayout {
                id: windowRow
                readonly property int rowNumber: model.number
                readonly property bool isOpen: model.open
                readonly property var cells: [
                    { key: "phi", value: model.phi }, { key: "z", value: model.z },
                    { key: "width", value: model.width }, { key: "height", value: model.height },
                    { key: "albedo", value: model.albedo }, { key: "conversion", value: model.conversion }]
                spacing: 4
                Label { text: model.name; Layout.preferredWidth: 64 }
                Repeater {
                    model: windowRow.cells
                    TextField {
                        Layout.preferredWidth: 52
                        text: modelData.value
                        selectByMouse: true
                        // albedo and conversion only apply to a closed window
                        enabled: !(windowRow.isOpen && index >= 4)
                        onEditingFinished: if (text !== modelData.value)
                                               Julia.br_edit_window(windowRow.rowNumber, modelData.key, text)
                    }
                }
                CheckBox {
                    checked: windowRow.isOpen
                    onToggled: Julia.br_window_open(windowRow.rowNumber, checked)
                    ToolTip.visible: hovered
                    ToolTip.text: "Open: a hole in the wall (burned through or uncovered)"
                }
                ToolButton {
                    text: "✕"
                    onClicked: Julia.br_remove_window(windowRow.rowNumber)
                    ToolTip.visible: hovered
                    ToolTip.text: "Remove this window"
                }
            }
        }
    }

    Label {
        text: app.targetError
        visible: text !== ""
        color: "#c42b1c"
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }
    Label {
        text: app.targetSummary
        opacity: 0.75
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
        Layout.topMargin: 8
    }
    Label {
        visible: app.targetKind === "foil"
        text: "The foil sits on the +x axis facing TCC; the black sphere at TCC stands in " +
              "for the main target (it absorbs and emits nothing)."
        opacity: 0.75
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }
}
