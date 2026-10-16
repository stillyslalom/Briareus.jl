import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml

StepPage {
    id: page
    title: "Beams"
    guidance: "Load a facility preset or a beam table (VisRad pointing and beam-summary files " +
              "apply to the loaded table), then edit pointing per cone, quad or " +
              "beam. Aim z is along the hohlraum axis (per cone, mirrored for the lower " +
              "hemisphere), r toward the beam's own port azimuth, and the tangential offset " +
              "counterclockwise about z (mm). Positive defocus puts best focus beyond the aim " +
              "point. Power is per beam (TW). The viewer draws each beam from its port to its aim."

    RowLayout {
        spacing: 8
        ComboBox {
            id: presetBox
            Layout.preferredWidth: 190
            textRole: "text"; valueRole: "value"
            model: [{ text: "NIF quads (48)", value: "nif_quads" },
                    { text: "NIF beams (192)", value: "nif_beams" },
                    { text: "OMEGA hohlraum (40)", value: "omega_hohlraum" },
                    { text: "OMEGA backlighter", value: "omega_backlighter" }]
            currentIndex: Math.max(0, ["nif_quads", "nif_beams", "omega_hohlraum",
                                       "omega_backlighter"].indexOf(app.preset))
        }
        Button { text: "Load"; onClicked: Julia.br_load_preset(presetBox.currentValue) }
        ToolSeparator {}
        Button {
            text: "Import…"
            onClicked: importDialog.open()
            ToolTip.visible: hovered
            ToolTip.text: "Read a beam table CSV (canonical columns of write_beam_table)"
        }
        Button { text: "Export…"; onClicked: exportDialog.open() }
    }
    Flow {
        spacing: 8
        Layout.fillWidth: true
        Label { text: "VisRad files"; height: 32; verticalAlignment: Text.AlignVCenter }
        Button {
            text: "Import SST…"
            enabled: table.count > 0
            onClicked: { visradDialog.mode = "import_sst"; visradDialog.open() }
            ToolTip.visible: hovered
            ToolTip.text: "Apply a NIF SST pointing CSV (aim points and defocus) to this table; " +
                          "beams absent from the file are switched off"
        }
        Button {
            text: "Import summary…"
            enabled: table.count > 0
            onClicked: { visradDialog.mode = "import_summary"; visradDialog.open() }
            ToolTip.visible: hovered
            ToolTip.text: "Apply a VisRad Laser Beam Summary CSV (pointing, defocus, power, spot) " +
                          "to this table; beams off or absent in the file are switched off"
        }
        Button {
            text: "Export SST…"
            enabled: table.count > 0
            onClicked: { visradDialog.mode = "export_sst"; visradDialog.open() }
            ToolTip.visible: hovered
            ToolTip.text: "Write the enabled beams' pointing as a NIF SST CSV (VisRad: Import Laser Beam Data, Pointing Data)"
        }
        Button {
            text: "Export summary…"
            enabled: table.count > 0
            onClicked: { visradDialog.mode = "export_summary"; visradDialog.open() }
            ToolTip.visible: hovered
            ToolTip.text: "Write the table in VisRad's Laser Beam Summary layout"
        }
    }
    Label { text: app.beamsSummary; opacity: 0.75 }
    Label {
        text: app.beamsWarnings
        visible: text !== ""
        color: "#d08a00"
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }

    RowLayout {
        spacing: 8
        Label { text: "Rows" }
        ComboBox {
            Layout.preferredWidth: 160
            textRole: "text"; valueRole: "value"
            model: [{ text: "per cone", value: "group" }, { text: "per quad", value: "quad" },
                    { text: "per beam", value: "beam" }]
            currentIndex: Math.max(0, ["group", "quad", "beam"].indexOf(app.beamsBy))
            onActivated: Julia.br_set_grouping(currentValue)
        }
    }
    Label {
        text: app.beamsError
        visible: text !== ""
        color: "#c42b1c"
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }

    // header + rows (a ListView keeps 192 beam rows cheap)
    readonly property var columns: [
        { key: "z", label: "z" }, { key: "r", label: "r" },
        { key: "phi", label: "tang." }, { key: "defocus", label: "defocus" },
        { key: "power", label: "TW" }]
    RowLayout {
        spacing: 4
        Label { text: "On"; Layout.preferredWidth: 36 }
        Label { text: "Name"; Layout.preferredWidth: 64 }
        Label { text: "n"; Layout.preferredWidth: 22 }
        Label { text: "θ / φ (°)"; Layout.preferredWidth: 84 }
        Repeater {
            model: page.columns
            Label { text: modelData.label; Layout.preferredWidth: 60; font.pixelSize: 12 }
        }
    }
    ListView {
        id: table
        Layout.fillWidth: true
        Layout.preferredHeight: Math.min(contentHeight, 460)
        clip: true
        model: beamModel
        ScrollBar.vertical: ScrollBar {}
        delegate: RowLayout {
            id: beamRow
            readonly property int rowNumber: model.number
            // the editable cells (the inner Repeater's own `model` hides the row's)
            readonly property var cells: [
                { key: "z", value: model.z }, { key: "r", value: model.r },
                { key: "phi", value: model.phiOffset }, { key: "defocus", value: model.defocus },
                { key: "power", value: model.power }]
            spacing: 4
            CheckBox {
                Layout.preferredWidth: 36
                checked: model.on !== 0
                tristate: false
                opacity: model.on === 2 ? 0.6 : 1
                onToggled: Julia.br_set_beam_on(model.number, checked)
                ToolTip.visible: hovered && model.on === 2
                ToolTip.text: "Some beams of this row are off"
            }
            Label { text: model.key; Layout.preferredWidth: 64; elide: Text.ElideRight }
            Label { text: model.count; Layout.preferredWidth: 22; opacity: 0.75 }
            Label { text: model.theta + " / " + model.phi; Layout.preferredWidth: 84; opacity: 0.75 }
            Repeater {
                model: beamRow.cells
                TextField {
                    Layout.preferredWidth: 60
                    text: modelData.value
                    selectByMouse: true
                    onEditingFinished: if (text !== modelData.value)
                                           Julia.br_edit_beam(beamRow.rowNumber, modelData.key, text)
                }
            }
        }
    }

    FileDialog {
        id: importDialog
        title: "Import a beam table"
        nameFilters: ["Beam table (*.csv *.tsv *.txt)", "All files (*)"]
        onAccepted: Julia.br_import_beams(selectedFile.toString())
    }
    FileDialog {
        id: exportDialog
        title: "Export the beam table"
        fileMode: FileDialog.SaveFile
        defaultSuffix: "csv"
        nameFilters: ["Beam table (*.csv)"]
        onAccepted: Julia.br_export_beams(selectedFile.toString())
    }
    // one dialog for the four VisRad files
    FileDialog {
        id: visradDialog
        property string mode: "import_sst"
        readonly property bool saving: mode.startsWith("export")
        title: {
            if (mode === "import_sst") return "Apply a NIF SST pointing file"
            if (mode === "import_summary") return "Apply a VisRad Laser Beam Summary"
            if (mode === "export_sst") return "Write the pointing as a NIF SST file"
            return "Write a VisRad Laser Beam Summary"
        }
        fileMode: saving ? FileDialog.SaveFile : FileDialog.OpenFile
        defaultSuffix: "csv"
        nameFilters: ["CSV (*.csv *.txt)", "All files (*)"]
        onAccepted: {
            const url = selectedFile.toString()
            if (mode === "import_sst") Julia.br_import_sst(url)
            else if (mode === "import_summary") Julia.br_import_summary(url)
            else if (mode === "export_sst") Julia.br_export_sst(url)
            else Julia.br_export_summary(url)
        }
    }
}
