// Properties of the selected beams (a group, quad, beam, or a Ctrl-selected set): on/off,
// pointing relative to the ports or at a common aim point, defocus and power for all of
// them at once (only the edited field changes); their group; splitting and merging
// NIF quads.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note { text: app.bnInfo; opacity: 0.9 }
    CheckBox {
        text: app.bnOn === 2 ? "On (some beams are off)" : "On"
        checked: app.bnOn !== 0
        onToggled: Julia.br_tree_on(app.selection, checked)
    }

    SectionLabel { text: "Pointing" }
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
    Note {
        text: "Values are the mean over the selection; an edit sets that field on every beam " +
              "and leaves the others. The beam direction is fixed by its port."
    }

    SectionLabel { text: "Aim point (x, y, z)" }
    Note {
        text: "Point every selected beam at one point, as VisRad points a cluster at a " +
              "reference point (e.g. a backlighter foil). Shown: the mean aim."
    }
    PropRow { label: "x"; value: app.bn_aim_x; unit: "mm"; onEdited: (text) => Julia.br_edit_node("aim_x", text) }
    PropRow { label: "y"; value: app.bn_aim_y; unit: "mm"; onEdited: (text) => Julia.br_edit_node("aim_y", text) }
    PropRow { label: "z"; value: app.bn_aim_z; unit: "mm"; onEdited: (text) => Julia.br_edit_node("aim_z", text) }
    RowLayout {
        Label { text: "Aim at"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.aimChoices
            // -1: the selected beams follow different components
            currentIndex: app.bnAimRef
            displayText: app.bnAimRef < 0 ? "(several)" : currentText
            onActivated: (index) => Julia.br_aim_at(index)
        }
    }
    Note {
        text: "Aiming at a component (its centroid: a foil, an entrance hole, a part) keeps " +
              "the beams on it when the target changes; pointing edits become offsets from it."
    }
    ErrorText { text: app.beamsError }

    SectionLabel { text: "Group" }
    Note {
        text: app.bnCluster === "" ? "The selection spans several groups."
              : "In group " + app.bnCluster + ". Groups gather beams for pointing, power and " +
                "on/off together, and are the optimizer's per-group parameters."
    }
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

    SectionLabel { text: "Quad"; visible: app.bnCanSplit || app.bnCanMerge }
    Note {
        visible: app.bnCanSplit || app.bnCanMerge
        text: app.bnCanSplit
              ? "Quad " + app.bnQuad + " is one beam. Split it into its four beams (own ports, " +
                "a quarter of the power each) to point them apart; the drive is then recomputed."
              : "Quad " + app.bnQuad + " is split into its four beams. Merge them back into one beam " +
                "(mean aim and defocus, total power)."
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
