// Properties of the selected beams (a group, quad, beam, or a Ctrl-selected set): their
// 1ω, 2ω and 3ω light in the calculations and in the view,
// pointing relative to the ports or at a common aim point in their aim frame, defocus and
// power for all of them at once (only the edited field changes); their group; splitting
// and merging NIF quads.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note { text: app.bnInfo; opacity: 0.9 }
    // the beams' light by harmonic: in the calculations, drawn (the model tree's dots)
    GridLayout {
        id: matrix
        columns: 4
        columnSpacing: 6
        rowSpacing: 0
        component Box: CheckBox {
            property int at: 0                       // index into app.bnMatrix
            readonly property string st: app.bnMatrix.charAt(at)
            checked: st === "1" || st === "2"
            opacity: st === "2" ? 0.55 : 1
            padding: 2
            onToggled: Julia.br_beam_dot(app.selection, at < 3 ? "use" : "draw", at % 3 + 1, checked)
        }
        Label { text: "" }
        Label { text: "1ω"; color: "#d62728"; font.bold: true; Layout.alignment: Qt.AlignHCenter }
        Label { text: "2ω"; color: "#2ca02c"; font.bold: true; Layout.alignment: Qt.AlignHCenter }
        Label { text: "3ω"; color: "#1f5fd6"; font.bold: true; Layout.alignment: Qt.AlignHCenter }
        Label { text: "In calculations" }
        Box { at: 0; enabled: st !== "-" && app.ucLive[0]
              tip: st === "-" ? "No 1ω footprint for these beams"
                 : app.ucLive[0] ? "1ω light in the drive" : "Set the 1ω power ratio in Drive first" }
        Box { at: 1; enabled: st !== "-" && app.ucLive[1]
              tip: st === "-" ? "No 2ω footprint for these beams"
                 : app.ucLive[1] ? "2ω light in the drive" : "Set the 2ω power ratio in Drive first" }
        Box { at: 2; tip: app.bnOn === 2 ? "Some of these beams are off" : "The beams fire" }
        Label { text: "Drawn" }
        Box { at: 3; enabled: st !== "-"; tip: st === "-" ? "No 1ω footprint for these beams" : "Draw the 1ω footprint" }
        Box { at: 4; enabled: st !== "-"; tip: st === "-" ? "No 2ω footprint for these beams" : "Draw the 2ω footprint" }
        Box { at: 5; tip: "Draw the beams (drawing only)" }
    }

    SectionLabel { text: "Pointing"; info: "beams.pointing" }
    PropRow { label: "Aim z"; value: app.bn_z
              unit: app.bnMirrored ? "mm, mirrored below TCC" : "mm"
              tip: "Aim along the hohlraum axis; mirrored for the lower beams"
              onEdited: (text) => Julia.br_edit_node("z", text) }
    PropRow { label: "Aim r"; value: app.bn_r; unit: "mm toward the port azimuth"
              tip: "Radial aim offset, in the port azimuth (on NIF: the quad's)"
              onEdited: (text) => Julia.br_edit_node("r", text) }
    PropRow { label: "Tangential offset"; value: app.bn_phi; unit: "mm, counterclockwise"
              tip: "Aim offset across that azimuth"
              onEdited: (text) => Julia.br_edit_node("phi", text) }
    PropRow { label: "Defocus"; value: app.bn_defocus; unit: "mm (+: beyond the aim)"
              tip: "Best focus, this far beyond the aim point"
              onEdited: (text) => Julia.br_edit_node("defocus", text) }
    PropRow { label: "Power"; value: app.bn_power; unit: app.bnPowerUnit === "quad" ? "TW per quad" : "TW per beam"
              tip: app.bnPowerUnit === "quad" ? "3ω power of each selected quad, shared by its beams"
                                              : "3ω power of each selected beam"
              onEdited: (text) => Julia.br_edit_node("power", text) }

    SectionLabel { text: "Aim point (x, y, z)"; info: "beams.aim" }
    RowLayout {
        Label { text: "Aim frame"; Layout.preferredWidth: 110 }
        ComboBox {
            id: aimFrameBox
            objectName: "aimFrame"
            Layout.fillWidth: true
            model: app.aimChoices
            // -1: the selected beams aim in different frames
            currentIndex: app.bnAimFrame
            displayText: app.bnAimFrame < 0 ? "(several)" : currentText
            onActivated: (index) => Julia.br_aim_frame(index)
            tip: "The frame the aim is held in; the beams follow it"
        }
    }
    Button {
        text: "Aim at the frame's center"
        enabled: app.bnAimFrame >= 0
        onClicked: Julia.br_aim_at_frame()
        tip: "Point the beams at the frame's origin"
    }
    PropRow { label: "x"; value: app.bn_aim_x; unit: "mm"; tip: "In " + app.bnAimFrameLabel + "; shown: the mean aim"
              onEdited: (text) => Julia.br_edit_node("aim_x", text) }
    PropRow { label: "y"; value: app.bn_aim_y; unit: "mm"; tip: "In " + app.bnAimFrameLabel + "; shown: the mean aim"
              onEdited: (text) => Julia.br_edit_node("aim_y", text) }
    PropRow { label: "z"; value: app.bn_aim_z; unit: "mm"; tip: "In " + app.bnAimFrameLabel + "; shown: the mean aim"
              onEdited: (text) => Julia.br_edit_node("aim_z", text) }
    ErrorText { text: app.beamsError }

    SectionLabel { text: "Group"; info: "beams.group" }
    Note { text: app.bnCluster === "" ? "The selection spans several groups." : "In group " + app.bnCluster + "." }
    RowLayout {
        spacing: 8
        Layout.fillWidth: true
        ComboBox {
            id: groupBox
            Layout.fillWidth: true
            Layout.minimumWidth: 100
            editable: true
            model: app.clusterNames
            currentIndex: model ? model.indexOf(app.bnCluster) : -1
            tip: "An existing group, or a new name"
        }
        Button {
            text: "Move"
            enabled: groupBox.editText.trim() !== "" && groupBox.editText !== app.bnCluster
            onClicked: Julia.br_set_cluster(groupBox.editText)
            tip: "Move the selected beams to this group"
        }
    }
    PropRow {
        visible: app.selKind === "group" && app.selectionCount === 1
        label: "Rename group"
        value: app.selKey
        fieldWidth: 130
        tip: "A new name for this group"
        onEdited: (text) => Julia.br_rename_cluster(text)
    }

    // NIF quads: held as one beam each (a split gives the deposition per beam), or as four
    SectionLabel { text: "Quads"; info: "beams.quad"; visible: app.bnSplit > 0 || app.bnMerge > 0 }
    Flow {
        Layout.fillWidth: true
        spacing: 6
        Button {
            visible: app.bnSplit > 0
            text: app.bnSplit === 1 ? "Split into beams" : "Split " + app.bnSplit + " quads into beams"
            onClicked: Julia.br_split_quads("")
            tip: "Four beams per quad, each from its own port"
        }
        Button {
            visible: app.bnMerge > 0
            text: app.bnMerge === 1 ? "Merge into quad" : "Merge into " + app.bnMerge + " quads"
            onClicked: Julia.br_merge_quads("")
            tip: "One beam per quad again"
        }
    }

    SectionLabel { text: "View" }
    Flow {
        Layout.fillWidth: true
        spacing: 6
        Button {
            text: "Look along beam"
            onClicked: Julia.br_look_along_beam()
            tip: "View the target from the beam's port"
        }
        Button { text: "Fit selection"; onClicked: Julia.br_fit_selection(); tip: "Fit the 3D view to the selected beams" }
    }
}
