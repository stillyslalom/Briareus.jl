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
    // the diagnostic picker or the image list, opened (`request_popup`)
    function openPicker() { picker.open(); return picker }
    function openImages() { imagesPopup.open(); return imagesPopup }

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
                        { key: "irradiance", label: "Irradiance" }, { key: "retune", label: "Retune" }]
                delegate: Button {
                    text: modelData.label
                    // as wide as the name (not the style's minimum), so that the row fits
                    // a laptop-size window on one line
                    implicitWidth: implicitContentWidth + leftPadding + rightPadding
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
                text: app.timeChecking ? "Checking…" : "Run checks"
                visible: app.timeChart === "irradiance"
                enabled: app.timeHasPulses && !app.driveBuilding && !app.timeChecking
                onClicked: Julia.br_time_checks()
                ToolTip.visible: hovered
                ToolTip.text: "Irradiance histories of the unconverted 1ω and 2ω light (unconverted_history) and the specular and pre-plasma checks against the plasma threshold (Drive properties)"
            }
            Button {
                text: "Run history"
                visible: !app.timeRunning && app.timeChart !== "irradiance"
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
        // the irradiance checks: a pass or fail dot and the first crossing each (the
        // parts' peaks in the tooltip)
        Flow {
            visible: app.timeChart === "irradiance" && app.irrTexts.length > 0
            Layout.fillWidth: true
            spacing: 12
            Repeater {
                model: app.irrTexts.length
                delegate: Row {
                    spacing: 5
                    StatusDot {
                        anchors.verticalCenter: parent.verticalCenter
                        status: app.irrStale ? "todo" : app.irrPass[index] ? "ok" : "attention"
                    }
                    Label {
                        text: app.irrTexts[index]
                        HoverHandler { id: checkHover }
                        ToolTip.visible: checkHover.hovered
                        ToolTip.text: app.irrTips[index]
                    }
                }
            }
            Info { topic: "time.irradiance" }
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
        // one line: the diagnostics recorded in the next history (the ticked definitions
        // as chips, all of them in the picker), then the images the history holds of its
        // recorded imagers (each opens in the image view)
        Flow {
            Layout.fillWidth: true
            spacing: 4
            Label { text: "Record"; opacity: 0.75; height: 26; verticalAlignment: Text.AlignVCenter }
            Button {
                id: pickerButton
                text: "Diagnostics ▾"
                onClicked: picker.opened ? picker.close() : picker.open()
                ToolTip.visible: hovered && !picker.opened
                ToolTip.text: "The diagnostic definitions found (built-in, the folders', the session's); each ticked one is recorded in the next history"
            }
            Repeater {
                model: diagModel
                delegate: Rectangle {
                    required property string id
                    required property string name
                    required property bool ticked
                    required property bool usable
                    required property bool changed
                    required property string tip
                    visible: ticked
                    height: 26
                    width: chipRow.implicitWidth + 10
                    radius: 13
                    color: "transparent"
                    border.color: palette.mid
                    Row {
                        id: chipRow
                        anchors.verticalCenter: parent.verticalCenter
                        x: 8
                        Label {
                            text: name + (changed ? " *" : "")
                            anchors.verticalCenter: parent.verticalCenter
                            opacity: usable ? 1 : 0.5
                            HoverHandler { id: chipHover }
                            ToolTip.visible: chipHover.hovered
                            ToolTip.text: tip
                        }
                        ToolButton {
                            visible: changed
                            text: "↺"
                            implicitHeight: 24
                            onClicked: Julia.br_diag_use_path(id)
                            ToolTip.visible: hovered
                            ToolTip.text: "Use the search path's definition of " + id + " instead of the session's"
                        }
                        ToolButton {
                            text: "✕"
                            implicitHeight: 24
                            font.pixelSize: 10
                            onClicked: Julia.br_time_diagnostic(id, false)
                            ToolTip.visible: hovered
                            ToolTip.text: "Do not record " + name
                        }
                    }
                }
            }
            ToolButton {
                text: "Folders…"
                onClicked: win.openDiagnostics()
                ToolTip.visible: hovered
                ToolTip.text: "Diagnostic definition folders, and the problems found in them"
            }
            ToolSeparator { visible: imageRepeater.count > 0; height: 26 }
            Button {
                id: imagesButton
                visible: imageRepeater.count > 0
                text: "Images (" + imageRepeater.count + ") ▾"
                onClicked: imagesPopup.opened ? imagesPopup.close() : imagesPopup.open()
                ToolTip.visible: hovered && !imagesPopup.opened
                ToolTip.text: "The images the history holds of its recorded imagers; each opens in the image view, " +
                              "zoomed to its region with signal (Ctrl + Shift + click there: the whole field; Ctrl + click: back)"
            }
        }
        // the images of the recorded imagers, one entry each (opens above its button)
        Popup {
            id: imagesPopup
            parent: imagesButton
            x: 0
            y: -height - 2
            margins: 6
            padding: 4
            closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutsideParent
            width: imageList.implicitWidth + leftPadding + rightPadding + 14
            height: Math.min(imageList.implicitHeight + topPadding + bottomPadding,
                             (Overlay.overlay ? Overlay.overlay.height : 800) - 12)
            contentItem: ScrollView {
                clip: true
                contentWidth: imageList.implicitWidth
                Column {
                    id: imageList
                    Repeater {
                        id: imageRepeater
                        model: diagImageModel
                        delegate: ItemDelegate {
                            required property string id
                            required property string label
                            required property int frame
                            text: label
                            topPadding: 3
                            bottomPadding: 3
                            onClicked: { imagesPopup.close(); Julia.br_diag_image(id, frame) }
                            ToolTip.visible: hovered
                            ToolTip.text: frame === 0 ? "The detected radiance at the time cursor's step (W/cm²/sr)"
                                                      : "The detected fluence per steradian of the recording (J/cm²/sr)"
                        }
                    }
                }
            }
        }
        // every definition found, by kind, with a check box each (opens above its button,
        // kept inside the window)
        Popup {
            id: picker
            parent: pickerButton
            x: 0
            y: -height - 2
            margins: 6
            padding: 6
            closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutsideParent
            width: pickerList.implicitWidth + leftPadding + rightPadding + 14
            height: Math.min(pickerList.implicitHeight + topPadding + bottomPadding,
                             (Overlay.overlay ? Overlay.overlay.height : 800) - 12)
            contentItem: ScrollView {
                clip: true
                contentWidth: pickerList.implicitWidth
                Column {
                    id: pickerList
                    spacing: 0
                    Repeater {
                        model: diagModel
                        delegate: Column {
                            required property string id
                            required property string name
                            required property bool ticked
                            required property bool usable
                            required property bool changed
                            required property string tip
                            required property string group
                            required property bool head
                            required property int index
                            Label {
                                visible: head
                                text: group
                                font.weight: Font.DemiBold
                                opacity: 0.75
                                topPadding: index > 0 ? 6 : 0
                            }
                            Row {
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
                    }
                }
            }
        }
        Note { text: app.timeStatus; Layout.fillWidth: true }
    }
}
