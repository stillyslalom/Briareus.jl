import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

StepPage {
    title: "Target"
    guidance: "Choose a hohlraum or a backlighter foil and its dimensions (mm). The viewer " +
              "colors the mesh by component; drag to rotate it. The coarse mesh keeps the " +
              "view factors on the Drive step to a few seconds; use the standard mesh for " +
              "design numbers."

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
