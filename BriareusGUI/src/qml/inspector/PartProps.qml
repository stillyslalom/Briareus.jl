// Properties of a free part: name, size, placement (the frame it is placed in, center,
// axis angles, rotation, placement steps), a move along any frame's axis, material,
// sidedness, element size.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    RowLayout {
        Layout.fillWidth: true
        Note { text: app.partKind + " · " + app.compInfo; opacity: 0.9 }
        Info { topic: "part.visrad"; visible: app.partKind === "VisRad object"; Layout.alignment: Qt.AlignTop }
    }
    Note {
        visible: app.selectionCount > 1
        text: app.selectionCount + " components selected: edits apply to every selected part " +
              "that has the field (shown: the first)."
    }
    RowLayout {
        visible: app.partFile !== ""
        Layout.fillWidth: true
        Cell { text: "File: " + app.partFile; opacity: 0.72; Layout.fillWidth: true }
        Info { topic: "part.cad" }
    }
    Note {
        visible: app.partFile !== "" && app.partCopy
        text: "Read from the copy saved in the session" +
              (app.partSourceExists ? "." : " (the file is not on this computer).")
    }
    RowLayout {
        visible: app.partFile !== ""
        Button {
            text: "Relink…"
            ToolTip.visible: hovered
            ToolTip.text: "Read this file's parts from another file"
            onClicked: relinkDialog.open()
        }
        Button {
            text: "Read from the file"
            visible: app.partCopy
            enabled: app.partSourceExists
            ToolTip.visible: hovered
            ToolTip.text: "Read this file's parts from the file on disk instead of the saved copy"
            onClicked: Julia.br_part_from_source()
        }
    }
    FileDialog {
        id: relinkDialog
        title: "Read the part's CAD from another file"
        fileMode: FileDialog.OpenFile
        nameFilters: ["CAD files (*.stl *.STL *.vrw)", "All files (*)"]
        onAccepted: Julia.br_relink_part(selectedFile.toString())
    }
    PropRow {
        label: "Name"
        value: app.partName
        fieldWidth: 140
        onEdited: (text) => Julia.br_edit_part("name", text)
    }
    RowLayout {
        Label { text: "Role"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.roleLabels
            currentIndex: app.partRole
            onActivated: (index) => Julia.br_part_role(index)
        }
        Info { topic: "part.role" }
    }
    CheckLine {
        text: "Thin film: left out of the drive's cavity (a capsule tent, an LEH window)"
        checked: app.partFilm
        onToggled: Julia.br_part_film(checked)
    }
    CheckLine {
        text: "Include in calculations (laser, x-ray drive, clearance)"
        checked: app.compIncluded
        onToggled: Julia.br_node_include(app.selection, checked)
    }
    RowLayout {
        Label { text: "Placed in"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.frameLabels
            // -1 after a rejected change: show the part's frame again
            currentIndex: app.partFrame < 0 ? 0 : app.partFrame
            onActivated: (index) => Julia.br_part_frame(index)
        }
        Info { topic: "part.placement" }
    }
    // Move along [axis] of [frame] by [Δ] [−][+]: each click one edit (undone as one)
    RowLayout {
        id: moveRow
        objectName: "moveRow"
        spacing: 4
        Layout.fillWidth: true
        Label { text: "Move along"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.preferredWidth: 52
            model: ["x", "y", "z"]
            currentIndex: app.moveAxis
            onActivated: (index) => Julia.br_move_axis(index)
        }
        Label { text: "of" }
        ComboBox {
            Layout.fillWidth: true
            Layout.minimumWidth: 90
            model: app.moveFrames
            currentIndex: app.moveFrame
            onActivated: (index) => Julia.br_move_frame(index)
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "The frame whose axis the part moves along: a template's, or a part's own (its axis is z)"
        }
    }
    RowLayout {
        spacing: 4
        Layout.fillWidth: true
        Label { text: "by"; Layout.preferredWidth: 110; horizontalAlignment: Text.AlignRight }
        TextField {
            Layout.preferredWidth: 64
            text: app.moveStep
            selectByMouse: true
            onEditingFinished: if (text !== app.moveStep) Julia.br_move_step(text)
        }
        Label { text: "mm"; opacity: 0.7 }
        Button { text: "−"; Layout.preferredWidth: 34; onClicked: Julia.br_nudge_part(-1)
                 ToolTip.visible: hovered; ToolTip.text: "Move back along the axis (one edit)" }
        Button { text: "+"; Layout.preferredWidth: 34; onClicked: Julia.br_nudge_part(1)
                 ToolTip.visible: hovered; ToolTip.text: "Move forward along the axis (one edit)" }
    }
    Button {
        text: "Sweep along this axis…"
        Layout.leftMargin: 114
        onClicked: Julia.br_sweep_part()
        ToolTip.visible: hovered
        ToolTip.text: "Step the part along this axis and check each step (clearance, line of sight, capture)"
    }
    RowLayout {
        Label { text: "Position in"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: ["Cartesian (x, y, z)", "Cylindrical (r, φ, z)", "Spherical (r, θ, φ)"]
            currentIndex: app.partCoords
            onActivated: (index) => Julia.br_part_coords(index)
        }
    }
    Repeater {
        model: partFieldModel
        delegate: PropRow {
            label: model.label
            value: model.value
            unit: model.unit === "" ? "fraction" : model.unit
            onEdited: (text) => Julia.br_edit_part(model.key, text)
        }
    }
    PlacementSteps { visible: app.selectionCount <= 1 }
    CheckLine {
        text: "Two-sided (both sides radiate and absorb)"
        checked: app.partTwoSided
        enabled: !app.partAperture
        onToggled: Julia.br_part_flag("two_sided", checked)
    }
    CheckLine {
        text: "Inward (radiates from its inside, or a flat part from its back)"
        checked: app.partInward
        onToggled: Julia.br_part_flag("inward", checked)
    }
    CheckLine {
        text: "Aperture (transparent: rays pass, its flux is what escapes)"
        checked: app.partAperture
        onToggled: Julia.br_part_flag("aperture", checked)
    }
    DisplayRow {}
    ErrorText { text: app.targetError }
    Button { text: "Remove part"; onClicked: Julia.br_remove_selected() }
}
