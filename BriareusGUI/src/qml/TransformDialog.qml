// Tools ▸ Coordinate transform (after VisRad's coordinate transform tool): a point or a
// direction given in one frame, read in another, in any coordinates; a frame can be
// made from it.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "common"

Dialog {
    id: dialog
    title: "Coordinate transform"
    modal: false
    standardButtons: Dialog.Close
    readonly property var systems: ["Cartesian (x, y, z)", "Cylindrical (r, φ, z)", "Spherical (r, θ, φ)"]
    readonly property var names: [["x", "y", "z"], ["r", "φ", "z"], ["r", "θ", "φ"]]
    readonly property var units: [["mm", "mm", "mm"], ["mm", "°", "mm"], ["mm", "°", "°"]]
    onAboutToShow: Julia.br_tf_open()

    GridLayout {
        columns: 4
        columnSpacing: 6
        rowSpacing: 4
        Label { text: "From" }
        ComboBox {
            Layout.columnSpan: 3
            Layout.preferredWidth: 300
            model: app.tfFrames
            currentIndex: app.tfFrom
            onActivated: (index) => Julia.br_tf_set("from", index)
        }
        Label { text: "To" }
        ComboBox {
            Layout.columnSpan: 3
            Layout.preferredWidth: 300
            model: app.tfFrames
            currentIndex: app.tfTo
            onActivated: (index) => Julia.br_tf_set("to", index)
        }
        Label { text: "Of a" }
        ComboBox {
            Layout.columnSpan: 2
            Layout.fillWidth: true
            model: ["Point", "Direction"]
            currentIndex: app.tfDirection ? 1 : 0
            onActivated: (index) => Julia.br_tf_set("direction", index === 1)
        }
        Info { topic: "frames.transform" }

        Label { text: "Given in" }
        ComboBox {
            Layout.columnSpan: 3
            Layout.fillWidth: true
            model: dialog.systems
            currentIndex: app.tfInput
            onActivated: (index) => Julia.br_tf_set("input", index)
        }
        Label { text: "Result in" }
        ComboBox {
            Layout.columnSpan: 3
            Layout.fillWidth: true
            model: dialog.systems
            currentIndex: app.tfOutput
            onActivated: (index) => Julia.br_tf_set("output", index)
        }
        // the three coordinates in, and out
        Repeater {
            model: 3
            delegate: RowLayout {
                Layout.columnSpan: 4
                spacing: 4
                Label { text: dialog.names[app.tfInput][index]; Layout.preferredWidth: 56 }
                TextField {
                    Layout.preferredWidth: 90
                    text: app.tfIn[index]
                    selectByMouse: true
                    onEditingFinished: if (text !== app.tfIn[index]) Julia.br_tf_set("c" + (index + 1), text)
                }
                Label { text: dialog.units[app.tfInput][index]; opacity: 0.7; Layout.preferredWidth: 28 }
                Label { text: "→  " + dialog.names[app.tfOutput][index]; Layout.preferredWidth: 40 }
                // read-only but selectable, to copy a result
                TextEdit {
                    objectName: "tfOut" + index
                    text: app.tfOut[index]
                    readOnly: true
                    selectByMouse: true
                    color: palette.windowText
                    selectionColor: palette.highlight
                    selectedTextColor: palette.highlightedText
                    font.family: "monospace"
                    Layout.preferredWidth: Math.max(90, implicitWidth)
                    horizontalAlignment: Text.AlignRight
                }
                Label { text: dialog.units[app.tfOutput][index]; opacity: 0.7 }
            }
        }
        RowLayout {
            Layout.columnSpan: 4
            Button {
                text: "Create frame"
                onClicked: Julia.br_tf_create()
                ToolTip.visible: hovered
                ToolTip.delay: 400
                ToolTip.text: app.tfDirection ? "A frame in the To frame, at its origin, with z along the direction"
                                              : "A frame in the To frame at the point, with the From frame's axes"
            }
            Item { Layout.fillWidth: true }
        }
        ErrorText { text: app.tfError; Layout.columnSpan: 4 }
    }
}
