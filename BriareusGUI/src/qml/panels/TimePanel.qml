// Time: the pulses the beams fire, the drive history they give (the observed surface's
// Tr and symmetry), where the energy goes, the irradiance checks and the latest retune.
// The bottom row chooses what the 3D view shows: the drive at the powers set on the
// beams, or a time of the drive history (the slider, or a click or drag on the chart).
// The diagnostics a history records are in the Diagnostics tab.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import QtQml.Models
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
                    model: [{ key: "pulses", label: "Pulses", tip: "Laser power of each group against time" },
                            { key: "drive", label: "Drive history", tip: "Tr and symmetry through the pulse" },
                            { key: "energy", label: "Energy", tip: "Where the laser energy went" },
                            { key: "irradiance", label: "Irradiance", tip: "Unconverted light against the plasma threshold" },
                            { key: "retune", label: "Retune", tip: "The latest retune, before and after" }]
                    delegate: Button {
                        text: modelData.label
                        // as wide as the name (not the style's minimum), so that the row fits
                        // a laptop-size window on one line
                        implicitWidth: implicitContentWidth + leftPadding + rightPadding
                        flat: app.timeChart !== modelData.key
                        highlighted: app.timeChart === modelData.key
                        onClicked: Julia.br_time_chart(modelData.key)
                        tip: modelData.tip
                    }
                }
                // as tall as the buttons (the style's own is taller, and sets the row's height)
                ToolSeparator { height: 26 }
                // the pulses the beams fire: a pulse table, or a published shape
                Button {
                    id: pulsesButton
                    text: "Pulses ▾"
                    onClicked: pulsesMenu.popup(pulsesButton, 0, pulsesButton.height)
                    tip: app.pulsesLabel !== "" ? app.pulsesLabel : "Load a pulse table or a published shape"
                    AutoMenu {
                        id: pulsesMenu
                        MenuItem { text: "Pulse table (CSV)…"; onTriggered: pulseDialog.open() }
                        MenuSeparator {}
                        MenuItem { text: "Published shapes"; enabled: false }
                        Instantiator {
                            model: app.pulseShapeKeys
                            delegate: MenuItem {
                                required property string modelData
                                required property int index
                                text: app.pulseShapeLabels[index]
                                onTriggered: Julia.br_time_shape(modelData)
                            }
                            onObjectAdded: (index, object) => pulsesMenu.insertItem(3 + index, object)
                            onObjectRemoved: (index, object) => pulsesMenu.removeItem(object)
                        }
                    }
                }
                // a pulse table's cone or quad columns: the power of each quad (shared by
                // its beams, held as one beam or as four) or of each beam
                ComboBox {
                    visible: app.pulseUnit !== ""
                    width: 150
                    height: 26
                    textRole: "text"; valueRole: "value"
                    model: [{ text: "TW per quad", value: "quad" }, { text: "TW per beam", value: "beam" }]
                    currentIndex: app.pulseUnit === "beam" ? 1 : 0
                    onActivated: Julia.br_time_unit(currentValue)
                    tip: "What the table's cone and quad columns give"
                }
                Button {
                    text: app.timeChecking ? "Checking…" : "Run checks"
                    visible: app.timeChart === "irradiance"
                    enabled: app.timeHasPulses && !app.driveBuilding && !app.timeChecking
                    onClicked: Julia.br_time_checks()
                    tip: app.timeHasPulses ? "Check the unconverted light against the plasma threshold" : "Load pulses first"
                }
                Button {
                    text: "Run history"
                    visible: !app.timeRunning && app.timeChart !== "irradiance"
                    enabled: app.timeHasPulses && !app.driveBuilding
                    onClicked: Julia.br_time_run()
                    tip: app.timeHasPulses ? "Compute the drive at each time of the pulses" : "Load pulses first"
                }
                Button { text: "Cancel"; visible: app.timeRunning; onClicked: Julia.br_time_cancel(); tip: "Stop the history" }
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
                Info { topic: "time." + app.timeChart }
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
                            Tip { visible: checkHover.hovered; text: app.irrTips[index]; maxWidth: 420 }
                        }
                    }
                }
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
                        anchors.top: parent.top
                        anchors.left: parent.left
                        width: Math.min(parent.width - 20, 520)
                        text: app.timeEmpty
                    }
                }
                // the key: the charts' colors, as the model tree keys the 3D view's
                ColumnLayout {
                    id: timeKey
                    visible: app.timeEmpty === ""
                    Layout.alignment: Qt.AlignTop
                    Layout.maximumWidth: 260
                    spacing: 3
                    Repeater {
                        model: timeLegendModel
                        delegate: RowLayout {
                            id: entry
                            Layout.fillWidth: true
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
                            Cell { text: model.label; font.pixelSize: 12; Layout.fillWidth: true }
                        }
                    }
                }
            }

            // what the 3D view shows: the drive at the powers set on the beams, or a time of
            // the drive history (the slider, enabled in that mode)
            RowLayout {
                visible: app.timeHasHistory
                Layout.fillWidth: true
                spacing: 8
                Label { text: "3D view shows"; opacity: 0.75 }
                Button {
                    text: "Beam powers"
                    implicitWidth: implicitContentWidth + leftPadding + rightPadding
                    flat: app.timeMode !== "beams"
                    highlighted: app.timeMode === "beams"
                    onClicked: Julia.br_time_mode("beams")
                    tip: "The drive at the powers set on the beams"
                }
                Button {
                    text: "Drive history"
                    implicitWidth: implicitContentWidth + leftPadding + rightPadding
                    flat: app.timeMode !== "history"
                    highlighted: app.timeMode === "history"
                    enabled: !app.timeStale
                    onClicked: Julia.br_time_mode("history")
                    tip: app.timeStale ? "The history is out of date: run it again" : "The drive at a time of the history"
                }
                Slider {
                    id: stepSlider
                    Layout.fillWidth: true
                    from: 1
                    to: Math.max(app.timeSteps, 1)
                    stepSize: 1
                    snapMode: Slider.SnapAlways
                    enabled: !app.timeStale
                    opacity: app.timeMode === "history" ? 1 : 0.5
                    onMoved: Julia.br_time_step(value)
                }
                // follows the cursor however it moved (a click on the chart, the mode)
                Binding {
                    target: stepSlider
                    property: "value"
                    value: app.timeStep > 0 ? app.timeStep : 1
                    when: !stepSlider.pressed
                }
                // as wide as the widest label, so that the slider keeps its length
                Label {
                    text: app.timeMode === "history" ? app.timeLabel : ""
                    Layout.preferredWidth: 88
                    font.weight: Font.DemiBold
                }
                Info { topic: "time.cursor" }
            }
        }
    }
}
