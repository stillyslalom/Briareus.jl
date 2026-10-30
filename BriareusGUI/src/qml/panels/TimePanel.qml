// Time: the pulses the beams fire, the drive history they give (capsule Tr, the ticked
// synthetic diagnostics, a_l/a_0), where the energy goes, and the latest retune; a time
// cursor (click or drag the chart, or the slider) colors the 3D view by that step of the
// history.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import Makie
import jlqml
import "../common"

Item {
    function redraw() { timeCanvas.update() }

    FileDialog {
        id: pulseDialog
        title: "Pulse table (power per beam in TW against time_ns, one column per group, quad or beam)"
        fileMode: FileDialog.OpenFile
        nameFilters: ["Pulse table (*.csv)", "All files (*)"]
        onAccepted: Julia.br_time_pulses(selectedFile.toString())
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 6
        spacing: 4

        // the charts (as one segmented control), then the pulses and the run; the buttons
        // wrap to a second row on a narrow panel
        Flow {
            Layout.fillWidth: true
            spacing: 6
            Repeater {
                model: [{ key: "pulses", label: "Pulses" }, { key: "drive", label: "Drive history" },
                        { key: "diagnostics", label: "Diagnostics" }, { key: "energy", label: "Energy" },
                        { key: "retune", label: "Retune" }]
                delegate: Button {
                    text: modelData.label
                    flat: app.timeChart !== modelData.key
                    highlighted: app.timeChart === modelData.key
                    onClicked: Julia.br_time_chart(modelData.key)
                }
            }
            ToolSeparator {}
            Button {
                text: "Pulse table…"
                onClicked: pulseDialog.open()
                ToolTip.visible: hovered
                ToolTip.text: "Fire the beams with a pulse table (read_pulse_table); disabled beams stay off"
            }
            Button {
                text: "N221204 shape"
                onClicked: Julia.br_time_shape()
                ToolTip.visible: hovered
                ToolTip.text: "N221204's delivered pulse shape (Pak et al. 2024), each beam reaching its table power at the peak"
            }
            Button {
                text: "Run history"
                visible: !app.timeRunning
                enabled: app.timeHasPulses && !app.driveBuilding
                onClicked: Julia.br_time_run()
                ToolTip.visible: hovered
                ToolTip.text: "Solve the drive at each time of the pulses with the drive's view factors (drive_history)"
            }
            Button { text: "Cancel"; visible: app.timeRunning; onClicked: Julia.br_time_cancel() }
        }

        RowLayout {
            visible: app.timeRunning
            Layout.fillWidth: true
            ProgressBar { value: app.timeFraction; Layout.preferredWidth: 220 }
            Label { text: app.timeProgress; opacity: 0.75 }
        }

        Label {
            text: app.timeTitle
            font.weight: Font.DemiBold
            elide: Text.ElideRight
            Layout.fillWidth: true
        }

        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 10
            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true
                // the canvas lives as long as the window (see BriareusWindow.qml)
                MakieArea {
                    id: timeCanvas
                    objectName: "timeCanvas"
                    anchors.fill: parent
                    scene: app.timeCanvas
                    visible: app.timeEmpty === ""
                }
                Note {
                    anchors.centerIn: parent
                    width: Math.min(parent.width - 20, 520)
                    horizontalAlignment: Text.AlignHCenter
                    text: app.timeEmpty
                }
            }
            // the key: the charts' colors, as the model tree keys the 3D view's
            ColumnLayout {
                visible: app.timeEmpty === ""
                Layout.alignment: Qt.AlignTop
                Layout.preferredWidth: 190
                spacing: 3
                Repeater {
                    model: timeLegendModel
                    delegate: RowLayout {
                        id: entry
                        spacing: 6
                        readonly property string lineStyle: model.style
                        readonly property string lineColor: model.color
                        Row {
                            spacing: 2
                            Layout.preferredWidth: 22
                            Repeater {
                                model: entry.lineStyle === "solid" ? 1 : entry.lineStyle === "dash" ? 3 : 5
                                delegate: Rectangle {
                                    width: entry.lineStyle === "solid" ? 22 : entry.lineStyle === "dash" ? 6 : 2.4
                                    height: entry.lineStyle === "solid" ? 3 : 2
                                    color: entry.lineColor
                                }
                            }
                        }
                        Label { text: model.label; font.pixelSize: 12 }
                    }
                }
            }
        }

        RowLayout {
            visible: app.timeHasHistory
            Layout.fillWidth: true
            spacing: 8
            Label { text: "3D view at"; opacity: 0.75 }
            Slider {
                id: stepSlider
                Layout.fillWidth: true
                from: 1
                to: Math.max(app.timeSteps, 1)
                stepSize: 1
                snapMode: Slider.SnapAlways
                enabled: !app.timeStale
                onMoved: Julia.br_time_step(value)
            }
            // follows the cursor however it moved (a click on the chart, Now)
            Binding {
                target: stepSlider
                property: "value"
                value: app.timeStep > 0 ? app.timeStep : 1
                when: !stepSlider.pressed
            }
            Label {
                text: app.timeLabel !== "" ? app.timeLabel : "the table's powers"
                Layout.preferredWidth: 120
                font.weight: app.timeLabel !== "" ? Font.DemiBold : Font.Normal
            }
            Button { text: "Now"; enabled: app.timeLabel !== ""; onClicked: Julia.br_time_now()
                     ToolTip.visible: hovered; ToolTip.text: "Back to the drive at the beam table's powers" }
        }
        // the diagnostic definitions found (built-in, the folders', the session's): each
        // ticked one is recorded in the next history
        Flow {
            Layout.fillWidth: true
            spacing: 2
            Label { text: "Record"; opacity: 0.75; height: 28; verticalAlignment: Text.AlignVCenter }
            Repeater {
                model: diagModel
                delegate: Row {
                    required property string id
                    required property string name
                    required property bool ticked
                    required property bool usable
                    required property bool changed
                    required property string tip
                    CheckBox {
                        text: name + (changed ? " *" : "")
                        checked: ticked
                        enabled: usable
                        onToggled: Julia.br_time_diagnostic(id, checked)
                        ToolTip.visible: hovered
                        ToolTip.text: tip
                    }
                    ToolButton {
                        visible: changed
                        text: "↺"
                        onClicked: Julia.br_diag_use_path(id)
                        ToolTip.visible: hovered
                        ToolTip.text: "Use the search path's definition of " + id + " instead of the session's"
                    }
                }
            }
            ToolButton {
                text: "Folders…"
                onClicked: win.openDiagnostics()
                ToolTip.visible: hovered
                ToolTip.text: "Diagnostic definition folders, and the problems found in them"
            }
        }
        Note { text: app.timeStatus; Layout.fillWidth: true }
    }
}
