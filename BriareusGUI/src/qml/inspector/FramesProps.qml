// The Frames folder: add a user frame or a port frame, open the coordinate transform.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 6
    RowLayout {
        Layout.fillWidth: true
        Note { text: app.frameCount === 0 ? "No user frames" : app.frameCount + " user frame" + (app.frameCount === 1 ? "" : "s") }
        Info { topic: "frames.intro" }
    }
    RowLayout {
        Button {
            text: "Add frame"
            onClicked: Julia.br_add_frame()
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "A user frame at TCC, with the chamber's axes"
        }
        Button {
            text: "Transform…"
            onClicked: Julia.br_open_transform()
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "Coordinates of a point or direction between any two frames"
        }
    }
    RowLayout {
        Layout.fillWidth: true
        Label { text: "Port"; Layout.preferredWidth: 40 }
        ComboBox {
            id: portBox
            Layout.fillWidth: true
            model: app.portChoices
            enabled: app.portChoices.length > 0
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "The ports of the beam table (a NIF quad at the NIF port table's angles)"
        }
        Button {
            text: "Add port frame"
            enabled: app.portChoices.length > 0
            onClicked: Julia.br_add_port_frame(portBox.currentIndex)
        }
        Info { topic: "frames.port" }
    }
    ErrorText { text: app.targetError }
}
