// Properties of the beam set: presets, beam files (beam table CSV and VisRad files),
// and how the beams are drawn.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note { text: app.beamsSummary; opacity: 0.9 }

    SectionLabel { text: "Preset"; info: "beams.preset" }
    RowLayout {
        spacing: 8
        ComboBox {
            id: presetBox
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            model: [{ text: "NIF quads (48)", value: "nif_quads" },
                    { text: "NIF beams (192)", value: "nif_beams" },
                    { text: "OMEGA hohlraum (40)", value: "omega_hohlraum" },
                    { text: "OMEGA backlighter", value: "omega_backlighter" }]
            currentIndex: Math.max(0, ["nif_quads", "nif_beams", "omega_hohlraum",
                                       "omega_backlighter"].indexOf(app.preset))
        }
        Button { text: "Load"; onClicked: Julia.br_load_preset(presetBox.currentValue) }
    }

    SectionLabel { text: "Switch" }
    RowLayout {
        Button { text: "All on"; onClicked: Julia.br_tree_on("beams", true) }
        Button { text: "All off"; onClicked: Julia.br_tree_on("beams", false) }
    }

    SectionLabel { text: "Display"; info: "beams.display" }
    RowLayout {
        spacing: 8
        Label { text: "Draw beams as"; Layout.preferredWidth: 150 }
        ComboBox {
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            model: [{ text: "Rays (axis to landing)", value: "rays" },
                    { text: "Volumes (envelope)", value: "volumes" },
                    { text: "Hidden", value: "off" }]
            currentIndex: Math.max(0, ["rays", "volumes", "off"].indexOf(app.beamStyle))
            onActivated: Julia.br_beam_style(currentValue)
        }
    }
    RowLayout {
        spacing: 8
        Label { text: "Color beams by"; Layout.preferredWidth: 150 }
        ComboBox {
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            model: [{ text: "Cone", value: "group" },
                    { text: "Clearance status", value: "clearance" },
                    { text: "TaLIS status", value: "talis" }]
            currentIndex: Math.max(0, ["group", "clearance", "talis"].indexOf(app.beamColor))
            onActivated: Julia.br_beam_color(currentValue)
        }
    }
    PropRow {
        label: "Drawn length"
        value: app.beamLength
        unit: "mm before the aim (empty: auto)"
        onEdited: (text) => Julia.br_beam_length(text)
    }
    ErrorText { text: app.beamsError }
    Note { text: app.beamsWarnings; color: "#d08a00"; opacity: 1 }
}
