// Properties of a user frame: name, where it is, the frame it is placed in, its origin
// and axis angles (a port frame: its roll alone), placement steps, a move along any
// frame's axis, its axes drawn in the 3D view.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 6
    readonly property bool port: app.framePort !== ""
    RowLayout {
        Layout.fillWidth: true
        Note { text: (port ? "Port " + app.framePort + " · " : "") + app.frameInfo; opacity: 0.9 }
        Info { topic: port ? "frames.port" : "frames.intro"; Layout.alignment: Qt.AlignTop }
    }
    PropRow {
        label: "Name"
        value: app.partName
        fieldWidth: 140
        onEdited: (text) => Julia.br_edit_part("name", text)
    }
    RowLayout {
        visible: !port
        Label { text: "Placed in"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.frameLabels
            currentIndex: app.partFrame < 0 ? 0 : app.partFrame
            onActivated: (index) => Julia.br_part_frame(index)
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "The frame it is placed in: TCC, another frame, a part's or a piece's (it stays where it is)"
        }
    }
    RowLayout {
        visible: !port
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
        visible: !port
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
    RowLayout {
        visible: !port
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
    PlacementSteps { visible: !port && app.selectionCount <= 1 }
    CheckLine {
        text: "Draw its axes in the 3D view (also while not selected)"
        checked: app.frameTriad
        onToggled: Julia.br_set_triad("", checked)
        tip: "x red, y green, z blue, with the frame's name at its origin"
    }
    RowLayout {
        Button {
            text: "Frame in it"
            onClicked: Julia.br_add_frame()
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "Add a user frame placed in this one"
        }
        Button {
            text: "Transform…"
            onClicked: Julia.br_open_transform()
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "Coordinates of a point or direction between this frame and another"
        }
    }
    ErrorText { text: app.targetError }
    Button { text: "Remove frame"; onClicked: Julia.br_remove_selected() }
}
