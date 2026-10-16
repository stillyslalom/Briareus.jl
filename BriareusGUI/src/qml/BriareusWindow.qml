// The Briareus window: header, step rail, the step pages, the persistent viewer
// canvas with its toolbar, and the tick timer. All state lives in Julia: `app`
// (property map) and the item models (`stepModel`, `beamModel`, `modeModel`,
// `powerModel`, `optModeModel`, `paramModel`, `compModel`, `clearModel`,
// `worstModel`); actions call Julia.br_* functions.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Window
import QtQuick.Dialogs
import Makie
import jlqml
import "steps"

ApplicationWindow {
    id: win
    title: app.title
    width: 1480
    height: 920
    minimumWidth: 1100
    minimumHeight: 680
    visible: true

    readonly property var stepKeys: ["target", "beams", "drive", "optimize", "review"]

    // Give the canvas keyboard focus when the pointer enters it, unless a text
    // field is being edited (its typing would move to the canvas).
    function focusCanvas(canvas, hovered) {
        if (!hovered || canvas.activeFocus) return
        const item = canvas.Window.activeFocusItem
        if (item && item.hasOwnProperty("cursorPosition")) return
        canvas.forceActiveFocus(Qt.MouseFocusReason)
    }

    onClosing: Qt.quit()

    header: ToolBar {
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: 12
            anchors.rightMargin: 12
            spacing: 4
            Label { text: "Briareus"; font.weight: Font.DemiBold }
            Item { Layout.fillWidth: true }
            Label {
                text: app.status
                elide: Text.ElideRight
                Layout.maximumWidth: 640
                opacity: 0.75
            }
        }
    }

    SplitView {
        id: body
        anchors.fill: parent
        orientation: Qt.Horizontal

        Pane {
            SplitView.preferredWidth: 780
            SplitView.minimumWidth: 620
            padding: 0
            clip: true

            RowLayout {
                anchors.fill: parent
                spacing: 0

                ListView {
                    id: rail
                    Layout.fillWidth: false
                    Layout.preferredWidth: 210
                    Layout.fillHeight: true
                    Layout.topMargin: 8
                    model: stepModel
                    interactive: false
                    spacing: 2
                    delegate: StepDelegate {
                        width: rail.width
                        stepNumber: index + 1
                        current: model.key === app.step
                    }
                }

                ToolSeparator { orientation: Qt.Vertical; Layout.fillHeight: true }

                StackLayout {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    // pages scroll; their content must not widen the pane
                    Layout.preferredWidth: 0
                    Layout.minimumWidth: 0
                    currentIndex: Math.max(0, win.stepKeys.indexOf(app.step))
                    TargetPage {}
                    BeamsPage {}
                    DrivePage {}
                    OptimizePage {}
                    ReviewPage {}
                }
            }
        }

        // The canvas item lives as long as the window: QMLMakie canvases must not be
        // destroyed while the application runs.
        Item {
            SplitView.fillWidth: true
            SplitView.minimumWidth: 400
            MakieArea {
                id: canvas
                anchors.left: parent.left
                anchors.top: parent.top
                anchors.bottom: parent.bottom
                anchors.right: tools.left
                scene: app.canvas
                HoverHandler { onHoveredChanged: win.focusCanvas(canvas, hovered) }
            }
            ViewerToolbar {
                id: tools
                anchors.right: parent.right
                anchors.top: parent.top
                anchors.bottom: parent.bottom
                onSaveView: saveViewDialog.open()
            }
        }
    }

    FileDialog {
        id: saveViewDialog
        title: "Save the viewer as an image"
        fileMode: FileDialog.SaveFile
        defaultSuffix: "png"
        nameFilters: ["PNG image (*.png)"]
        onAccepted: {
            const target = selectedFile
            canvas.grabToImage(function (img) { img.saveToFile(target) })
        }
    }

    Timer {
        interval: 16
        running: true
        repeat: true
        onTriggered: {
            const r = Julia.br_tick()
            if (app.grabPath !== "") {
                const path = app.grabPath
                Julia.br_grab_started()
                body.grabToImage(function (img) { img.saveToFile(path) })
            }
            if (r === 2) {
                win.close()
            } else if (r === 1) {
                canvas.update()
            }
        }
    }
}
