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
    id: panel
    // the height that shows every row and the chart at its minimum (ResultsPanel grows to
    // it when the tab opens); a shorter panel scrolls
    readonly property real fitHeight: column.implicitHeight + 2 * margin
    readonly property int margin: 6
    readonly property int chartMinimum: 175
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

    // the panel's rows, the chart keeping `chartMinimum`; a scroll bar (and the wheel, but no
    // flicking, which would take the chart's drags) when the panel is shorter than `fitHeight`
    Flickable {
        id: flick
        anchors.fill: parent
        interactive: false
        clip: true
        readonly property bool scrolls: panel.fitHeight > height + 0.5
        contentWidth: width
        contentHeight: Math.max(height, panel.fitHeight)
        ScrollBar.vertical: ScrollBar { id: vbar; policy: flick.scrolls ? ScrollBar.AlwaysOn : ScrollBar.AlwaysOff }
        WheelHandler {
            enabled: flick.scrolls
            onWheel: (event) => flick.contentY = Math.max(0, Math.min(flick.contentHeight - flick.height,
                                                                      flick.contentY - event.angleDelta.y / 2))
        }

        ColumnLayout {
            id: column
            x: panel.margin
            y: panel.margin
            width: flick.width - 2 * panel.margin - (flick.scrolls ? vbar.width : 0)
            height: flick.contentHeight - 2 * panel.margin
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
                // as tall as the buttons (the style's own is taller, and sets the row's height)
                ToolSeparator { height: 26 }
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

            // one line, wrapping on a narrow panel: the chart's title; on Irradiance, the checks
            // (a pass or fail dot and the first crossing each, the parts' peaks in the
            // tooltip); the status
            Flow {
                id: titleRow
                Layout.fillWidth: true
                spacing: 12
                Label {
                    text: app.timeTitle
                    font.weight: Font.DemiBold
                    width: Math.min(implicitWidth, titleRow.width)
                    wrapMode: Text.WordWrap
                }
                Repeater {
                    model: app.timeChart === "irradiance" ? app.irrTexts.length : 0
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
                Info { topic: "time.irradiance"; visible: app.timeChart === "irradiance" && app.irrTexts.length > 0 }
                Label {
                    text: app.timeStatus
                    visible: text !== ""
                    opacity: 0.72
                    width: Math.min(implicitWidth, titleRow.width)
                    wrapMode: Text.WordWrap
                }
            }

            RowLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                // a usable plot area (the axes' labels take about 45 px of it); a taller key
                // sets the height instead
                Layout.minimumHeight: Math.max(panel.chartMinimum, timeKey.visible ? timeKey.implicitHeight : 0)
                Layout.preferredHeight: Layout.minimumHeight
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
                    id: timeKey
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

            // one line (two on a narrow panel or with many diagnostics): the time the 3D view
            // shows (the slider takes the room the rest leaves); the diagnostics recorded in
            // the next history (the ticked definitions as chips, all of them in the picker);
            // the images the history holds of its recorded imagers (each opens in the image
            // view)
            GridLayout {
                id: bottomRow
                Layout.fillWidth: true
                columnSpacing: 4
                rowSpacing: 4
                readonly property bool oneLine: !stepRow.visible || width >= 320 + recordRow.implicitWidth + 20
                columns: oneLine ? 3 : 1
                RowLayout {
                    id: stepRow
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
                ToolSeparator { visible: app.timeHasHistory && bottomRow.oneLine; Layout.preferredHeight: 26 }
                Row {
                    id: recordRow
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
                        implicitHeight: 26
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
        }
    }
}
