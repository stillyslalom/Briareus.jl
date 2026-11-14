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

    SectionLabel { text: "Facility"; info: "beams.facility"; Layout.topMargin: 0 }
    ComboBox {
        Layout.fillWidth: true
        model: app.facilityLabels
        currentIndex: Math.max(0, app.facilityKeys.indexOf(app.facility))
        onActivated: (index) => Julia.br_set_facility(app.facilityKeys[index])
        tip: "Sets the beams, the diagnostics listed and the checks"
    }
    Note {
        visible: app.beamCount === 0
        text: "No beams: this facility has no built-in set. Import a beam table (File ▸ Import beam table)."
    }

    // the facility's presets (NIF, OMEGA)
    SectionLabel { text: "Preset"; info: "beams.preset"; visible: presetRow.visible }
    RowLayout {
        id: presetRow
        visible: presetBox.count > 0
        spacing: 8
        ComboBox {
            id: presetBox
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            readonly property var all: [{ text: "NIF quads (48)", value: "nif_quads", facility: "nif" },
                                        { text: "NIF beams (192)", value: "nif_beams", facility: "nif" },
                                        { text: "OMEGA hohlraum (40)", value: "omega_hohlraum", facility: "omega" },
                                        { text: "OMEGA backlighter", value: "omega_backlighter", facility: "omega" }]
            model: all.filter(p => p.facility === app.facility)
            currentIndex: Math.max(0, model.map(p => p.value).indexOf(app.preset))
            tip: "The facility's beam sets, aimed for this target"
        }
        Button { text: "Load"; onClicked: Julia.br_load_preset(presetBox.currentValue)
                 tip: "Replace the beams with this preset" }
    }

    SectionLabel { text: "Switch" }
    RowLayout {
        Button { text: "All on"; onClicked: Julia.br_tree_on("beams", true); tip: "Switch every beam on" }
        Button { text: "Split quads"; visible: app.hasQuads; onClicked: Julia.br_split_quads("beams")
                 tip: "Hold every quad as its four beams" }
        Button { text: "Merge quads"; visible: app.hasQuads; onClicked: Julia.br_merge_quads("beams")
                 tip: "Hold every quad as one beam" }
        Button { text: "All off"; onClicked: Julia.br_tree_on("beams", false)
                 tip: "Switch every beam off" }
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
            tip: "Rays, volumes, or hidden (drawing only)"
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
            tip: "By cone, or by status in the latest check"
        }
    }
    PropRow {
        label: "Drawn length"
        value: app.beamLength
        unit: "mm before the aim (empty: auto)"
        tip: "Length drawn before the aim point; empty: automatic"
        onEdited: (text) => Julia.br_beam_length(text)
    }
    ErrorText { text: app.beamsError }
    Note { text: app.beamsWarnings; color: "#d08a00"; opacity: 1 }
}
