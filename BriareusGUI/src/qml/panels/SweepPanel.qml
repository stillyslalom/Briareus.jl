// Placement sweep: a row per step (least clearance and where, failing and review pairs,
// the line of sight, the captured power, the objective) beside a chart of the least
// clearance and the blocked line of sight against the shift; a step can be applied.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import Makie
import jlqml
import "../common"

Item {
    function redraw() { sweepCanvas.update() }
    function statusColor(fails, reviews, los) {
        if (fails > 0 || los === "blocked") return "#c42b1c"
        if (reviews > 0) return "#d08a00"
        return palette.windowText
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 6
        spacing: 4
        RowLayout {
            Layout.fillWidth: true
            spacing: 8
            Cell { text: app.swTitle; font.weight: Font.DemiBold; Layout.fillWidth: true; Layout.maximumWidth: implicitWidth }
            Info { topic: "sweep.results" }
            Item { Layout.fillWidth: true }
            Label { text: "inputs changed since"; visible: app.swStale; color: "#d08a00" }
            Button {
                text: "Apply step"
                enabled: app.swHasResult && app.swSelected > 0 && !app.swStale && !app.swRunning
                onClicked: Julia.br_sw_apply(app.swSelected)
                ToolTip.visible: hovered
                ToolTip.text: "Move the part to the selected step (one edit, undone as one)"
            }
        }
        RowLayout {
            visible: app.swRunning
            Layout.fillWidth: true
            ProgressBar { value: app.swTotal > 0 ? app.swDone / app.swTotal : 0; Layout.preferredWidth: 220 }
            Label { text: "step " + app.swDone + " of " + app.swTotal; opacity: 0.75 }
        }
        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 10
            // the table
            ColumnLayout {
                Layout.fillHeight: true
                Layout.preferredWidth: head.totalWidth + 12
                Layout.maximumWidth: head.totalWidth + 12
                Layout.minimumWidth: 200
                spacing: 2
                HScroll {
                    visible: app.swHasResult
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    minimumContentWidth: head.totalWidth
                    spacing: 2
                    TableHeader {
                        id: head
                        visible: app.swHasResult
                        columns: [{ t: "Shift", w: 58 }, { t: "Clearance", w: 70 }, { t: "Worst", w: 110 },
                                  { t: "Fail", w: 32 }, { t: "Rev.", w: 32 }, { t: "LOS", w: 56 },
                                  { t: "Blocked", w: 58 }, { t: "Captured", w: 66 }]
                    }
                    ListView {
                        id: rowsView
                        visible: app.swHasResult
                        Layout.fillWidth: true
                        Layout.fillHeight: true
                        Layout.minimumHeight: 60
                        clip: true
                        model: sweepModel
                        boundsBehavior: Flickable.StopAtBounds
                        ScrollBar.vertical: ScrollBar {}
                        delegate: Rectangle {
                            width: rowsView.width
                            height: line.implicitHeight + 2
                            color: model.number === app.swSelected ? palette.highlight : "transparent"
                            readonly property color ink: model.number === app.swSelected ? palette.highlightedText
                                                         : statusColor(model.fails, model.reviews, model.los)
                            RowLayout {
                                id: line
                                spacing: 4
                                readonly property var w: head.widths
                                Cell { text: model.value; Layout.preferredWidth: line.w[0]; color: parent.parent.ink }
                                Cell { text: model.clearance; Layout.preferredWidth: line.w[1]; color: parent.parent.ink }
                                Cell { text: model.worst; Layout.preferredWidth: line.w[2]; color: parent.parent.ink }
                                Cell { text: model.fails; Layout.preferredWidth: line.w[3]; color: parent.parent.ink }
                                Cell { text: model.reviews; Layout.preferredWidth: line.w[4]; color: parent.parent.ink }
                                Cell { text: model.los === "" ? "–" : model.los; Layout.preferredWidth: line.w[5]; color: parent.parent.ink }
                                Cell { text: model.blocked; Layout.preferredWidth: line.w[6]; color: parent.parent.ink }
                                Cell { text: model.captured; Layout.preferredWidth: line.w[7]; color: parent.parent.ink }
                            }
                            MouseArea {
                                anchors.fill: parent
                                onClicked: Julia.br_sw_select(model.number)
                                onDoubleClicked: Julia.br_sw_apply(model.number)
                            }
                        }
                    }
                }
                Note {
                    visible: !app.swHasResult
                    text: app.swRunning ? "" : "Set up the sweep in Properties (select it in the model tree, or " +
                                               "Sweep along this axis… on a part) and run it."
                }
            }
            // the chart and its key
            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.minimumWidth: 200
                // the canvas lives as long as the window (see BriareusWindow.qml)
                MakieArea {
                    id: sweepCanvas
                    objectName: "sweepCanvas"
                    anchors.fill: parent
                    scene: app.sweepCanvas
                    visible: app.swEmpty === ""
                }
                Note {
                    anchors.centerIn: parent
                    width: Math.min(parent.width - 20, 420)
                    horizontalAlignment: Text.AlignHCenter
                    text: app.swEmpty
                }
            }
            ColumnLayout {
                visible: app.swEmpty === ""
                Layout.alignment: Qt.AlignTop
                Layout.preferredWidth: 130
                spacing: 3
                Repeater {
                    model: sweepLegendModel
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
    }
}
