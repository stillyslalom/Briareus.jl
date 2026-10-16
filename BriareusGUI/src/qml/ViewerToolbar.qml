// Viewer toolbar along the canvas's right edge: the cutaway of the hohlraum shell,
// reset, and saving the view as an image. Drag rotates the 3D view; the wheel zooms.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

Pane {
    id: bar
    signal saveView()
    implicitWidth: 76
    width: 76
    padding: 4
    background: Rectangle { color: bar.palette.window }

    ColumnLayout {
        anchors.fill: parent
        spacing: 2
        ToolButton {
            text: "Cutaway"
            checkable: true
            checked: app.cutaway
            Layout.fillWidth: true
            onClicked: Julia.br_set_cutaway(checked)
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "Cut away the quadrant of the hohlraum wall and caps facing the default view"
        }
        ToolSeparator { orientation: Qt.Horizontal; Layout.fillWidth: true }
        ToolButton {
            text: "Reset"
            Layout.fillWidth: true
            onClicked: Julia.br_reset_view()
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "Default viewing angle and zoom (drag rotates, the wheel zooms)"
        }
        ToolButton {
            text: "Save…"
            Layout.fillWidth: true
            onClicked: bar.saveView()
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "Save the viewer as a PNG image"
        }
        Item { Layout.fillHeight: true }
    }
}
