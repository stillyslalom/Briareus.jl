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
            property string tip: ""
            readonly property string st: app.bnMatrix.charAt(at)
            checked: st === "1" || st === "2"
            opacity: st === "2" ? 0.55 : 1
            padding: 2
            ToolTip.visible: hovered && tip !== ""
            ToolTip.delay: 400
            ToolTip.text: tip
            onToggled: Julia.br_beam_dot(app.selection, at < 3 ? "use" : "draw", at % 3 + 1, checked)
        }
        Label { text: "" }
        Label { text: "1ω"; color: "#d62728"; font.bold: true; Layout.alignment: Qt.AlignHCenter }
        Label { text: "2ω"; color: "#2ca02c"; font.bold: true; Layout.alignment: Qt.AlignHCenter }
        Label { text: "3ω"; color: "#1f5fd6"; font.bold: true; Layout.alignment: Qt.AlignHCenter }
        Label { text: "In calculations" }
        Box { at: 0; enabled: st !== "-" && app.ucLive[0]
              tip: st === "-" ? "No 1ω footprint for these beams (NIF quads and beams only)"
                 : app.ucLive[0] ? "The 1ω light's absorbed power adds to the drive"
                 : "Set the 1ω power ratio in the Drive's properties first" }
        Box { at: 1; enabled: st !== "-" && app.ucLive[1]
              tip: st === "-" ? "No 2ω footprint for these beams (NIF quads and beams only)"
                 : app.ucLive[1] ? "The 2ω light's absorbed power adds to the drive"
                 : "Set the 2ω power ratio in the Drive's properties first" }
        Box { at: 2; tip: app.bnOn === 2 ? "Some of these beams are off" : "Laser drive, clearance and optimization" }
        Label { text: "Drawn" }
        Box { at: 3; enabled: st !== "-"
              tip: st === "-" ? "No 1ω footprint for these beams (NIF quads and beams only)"
                              : "The 1ω footprint, and the light converging to it" }
        Box { at: 4; enabled: st !== "-"
              tip: st === "-" ? "No 2ω footprint for these beams (NIF quads and beams only)"
                              : "The 2ω footprint at the 3ω best-focus plane" }
        Box { at: 5; tip: "Hiding a beam keeps it in every calculation" }
    }

    SectionLabel { text: "Pointing"; info: "beams.pointing" }
    PropRow { label: "Aim z"; value: app.bn_z
              unit: app.bnMirrored ? "mm, mirrored below TCC" : "mm"
              tip: "Along the hohlraum axis; for beams from both hemispheres the lower ones aim at −z"
              onEdited: (text) => Julia.br_edit_node("z", text) }
    PropRow { label: "Aim r"; value: app.bn_r; unit: "mm toward each port azimuth"
              onEdited: (text) => Julia.br_edit_node("r", text) }
    PropRow { label: "Tangential offset"; value: app.bn_phi; unit: "mm, counterclockwise"
              onEdited: (text) => Julia.br_edit_node("phi", text) }
    PropRow { label: "Defocus"; value: app.bn_defocus; unit: "mm (+: beyond the aim)"
              onEdited: (text) => Julia.br_edit_node("defocus", text) }
    PropRow { label: "Power"; value: app.bn_power; unit: "TW per beam"
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
        }
    }
    Button {
        text: "Aim at the frame's center"
        enabled: app.bnAimFrame >= 0
        onClicked: Julia.br_aim_at_frame()
        ToolTip.visible: hovered
        ToolTip.text: "Point the beams at the frame's center (a part's centroid, a template frame's origin)"
    }
    Note { text: "Coordinates in " + app.bnAimFrameLabel + "; shown: the mean aim." }
    PropRow { label: "x"; value: app.bn_aim_x; unit: "mm"; onEdited: (text) => Julia.br_edit_node("aim_x", text) }
    PropRow { label: "y"; value: app.bn_aim_y; unit: "mm"; onEdited: (text) => Julia.br_edit_node("aim_y", text) }
    PropRow { label: "z"; value: app.bn_aim_z; unit: "mm"; onEdited: (text) => Julia.br_edit_node("aim_z", text) }
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
        }
        Button {
            text: "Move"
            enabled: groupBox.editText.trim() !== "" && groupBox.editText !== app.bnCluster
            onClicked: Julia.br_set_cluster(groupBox.editText)
            ToolTip.visible: hovered
            ToolTip.text: "Move the selected beams to this group (type a new name for a new group)"
        }
    }
    PropRow {
        visible: app.selKind === "group" && app.selectionCount === 1
        label: "Rename group"
        value: app.selKey
        fieldWidth: 130
        onEdited: (text) => Julia.br_rename_cluster(text)
    }

    SectionLabel { text: "Quad"; info: "beams.quad"; visible: app.bnCanSplit || app.bnCanMerge }
    Note {
        visible: app.bnCanSplit || app.bnCanMerge
        text: app.bnCanSplit ? "Quad " + app.bnQuad + " is one beam." : "Quad " + app.bnQuad + " is split into its four beams."
    }
    Flow {
        Layout.fillWidth: true
        spacing: 6
        Button { visible: app.bnCanSplit; text: "Split into beams"; onClicked: Julia.br_split_quad() }
        Button { visible: app.bnCanMerge; text: "Merge quad"; onClicked: Julia.br_merge_quad() }
    }

    SectionLabel { text: "View" }
    Flow {
        Layout.fillWidth: true
        spacing: 6
        Button {
            text: "Look along beam"
            onClicked: Julia.br_look_along_beam()
            ToolTip.visible: hovered
            ToolTip.text: "View the target from the port, down the beam axis"
        }
        Button { text: "Fit selection"; onClicked: Julia.br_fit_selection() }
    }
}
