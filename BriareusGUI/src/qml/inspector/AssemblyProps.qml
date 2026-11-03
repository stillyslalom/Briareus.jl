// Properties of an assembly: its members and frame on one line, name, the frame it is
// placed in, its origin and axis angles, placement steps, a move along any frame's axis,
// and the tree commands (ungroup, remove with its contents).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 6
    RowLayout {
        Layout.fillWidth: true
        Note { text: app.assemblyInfo; opacity: 0.9 }
        Info { topic: "assembly.intro"; Layout.alignment: Qt.AlignTop }
    }
    PropRow {
        label: "Name"
        value: app.partName
        fieldWidth: 140
        onEdited: (text) => Julia.br_edit_part("name", text)
    }
    RowLayout {
        Label { text: "Placed in"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.frameLabels
            currentIndex: app.partFrame < 0 ? 0 : app.partFrame
            onActivated: (index) => Julia.br_part_frame(index)
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "The frame its own frame is placed in (it stays where it is); its members move with it"
        }
    }
    RowLayout {
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
                 ToolTip.visible: hovered; ToolTip.text: "Move back along the axis, with its members (one edit)" }
        Button { text: "+"; Layout.preferredWidth: 34; onClicked: Julia.br_nudge_part(1)
                 ToolTip.visible: hovered; ToolTip.text: "Move forward along the axis, with its members (one edit)" }
        Info { topic: "part.triad" }
    }
    RowLayout {
        Label { text: "Origin in"; Layout.preferredWidth: 110 }
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
            unit: model.unit
            onEdited: (text) => Julia.br_edit_part(model.key, text)
        }
    }
    PlacementSteps { visible: app.selectionCount <= 1 }
    ErrorText { text: app.targetError }
    RowLayout {
        Button {
            text: "Ungroup"
            onClicked: Julia.br_ungroup_assembly()
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "Remove the assembly and keep its members, where they are, in its parent"
        }
        Button {
            text: "Remove with its contents"
            onClicked: Julia.br_remove_selected()
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "Remove the assembly and every part and assembly in it"
        }
    }
}
