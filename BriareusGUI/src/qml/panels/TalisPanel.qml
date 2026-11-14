// TaLIS report: the notes (short lines, wrapped to the panel), the per-rule summary,
// every check colored by status, CSV export.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 6
    function statusColor(s) {
        if (s === "fail") return "#c42b1c"
        if (s === "review") return "#d08a00"
        return "#2e9e44"
    }
    RowLayout {
        visible: !app.tlHasReport
        spacing: 6
        Note {
            Layout.fillWidth: false
            text: app.tlRunning ? "checking… " + app.tlElapsed : "No report yet: run TaLIS rules"
        }
        Info { topic: "talis.intro" }
    }
    RowLayout {
        visible: app.tlHasReport
        spacing: 8
        Cell { text: app.tlBrief; font.weight: Font.DemiBold; Layout.fillWidth: true; Layout.maximumWidth: implicitWidth }
        Info { topic: "talis.report" }
        Item { Layout.fillWidth: true }
        Label { text: "inputs changed since"; visible: app.tlStale; color: "#d08a00" }
        Button { text: "Export CSV…"; enabled: !app.tlRunning; onClicked: win.openTalisExport()
                 tip: "Write every check to a CSV file" }
    }
    // notes and the rule summary beside every check, split by a draggable handle
    SplitView {
        Layout.fillWidth: true
        Layout.fillHeight: true
        orientation: Qt.Horizontal
        ScrollView {
            id: notesView
            SplitView.preferredWidth: 330
            SplitView.minimumWidth: 200
            SplitView.fillWidth: !app.tlHasReport
            contentWidth: availableWidth
            contentHeight: notesColumn.implicitHeight
            clip: true
            ColumnLayout {
                id: notesColumn
                width: notesView.availableWidth - 12          // clear of the split handle
                spacing: 4
                // the tolerance note first, then the rule summary, then the other notes
                Repeater {
                    model: noteModel
                    delegate: RowLayout {
                        visible: model.kind === "tolerance"
                        Layout.fillWidth: true
                        spacing: 6
                        // (the ⓘ follows the text; a long note wraps to the panel)
                        Label { text: model.text; opacity: 0.72; wrapMode: Text.WordWrap
                                Layout.fillWidth: true; Layout.maximumWidth: implicitWidth }
                        Info { topic: "talis.tolerance"; Layout.alignment: Qt.AlignTop }
                        Item { Layout.fillWidth: true }
                    }
                }
                SectionLabel { text: "Rules"; visible: app.tlHasReport }
                // one rule in two lines: its counts, then its worst check
                Repeater {
                    model: ruleModel
                    delegate: ColumnLayout {
                        spacing: 0
                        Layout.fillWidth: true
                        Tip {
                            visible: ruleHover.hovered
                            text: model.requirement + "\nworst " + model.worst + " · " + model.where_
                        }
                        HoverHandler { id: ruleHover }
                        Label {
                            text: model.label + "  <font color='" + palette.windowText + "'>" +
                                  (model.fail > 0 ? "<font color='#c42b1c'>" + model.fail + " fail</font> · " : "") +
                                  (model.review > 0 ? "<font color='#d08a00'>" + model.review + " review</font> · " : "") +
                                  model.pass + " pass</font>"
                            textFormat: Text.StyledText
                            color: statusColor(model.status)
                            font.weight: Font.DemiBold
                            wrapMode: Text.WordWrap
                            Layout.fillWidth: true
                        }
                        Label {
                            text: "worst " + model.worst + " · " + model.where_
                            opacity: 0.7
                            wrapMode: Text.WordWrap
                            Layout.fillWidth: true
                        }
                    }
                }
                SectionLabel { text: "Notes"; visible: app.tlHasReport }
                Repeater {
                    model: noteModel
                    delegate: RowLayout {
                        visible: model.kind !== "tolerance"
                        Layout.fillWidth: true
                        spacing: 6
                        Label {
                            text: "· " + model.text
                            opacity: model.kind === "scope" || model.kind === "source" ? 0.7 : 1
                            wrapMode: Text.WordWrap
                            Layout.fillWidth: true
                            Layout.maximumWidth: implicitWidth
                        }
                        Info { topic: "talis.scope"; visible: model.kind === "scope"; Layout.alignment: Qt.AlignTop }
                        Item { Layout.fillWidth: true }
                    }
                }
            }
        }
        // every check
        ColumnLayout {
            visible: app.tlHasReport
            SplitView.fillWidth: app.tlHasReport
            SplitView.minimumWidth: 300
            spacing: 4
            HScroll {
                Layout.fillWidth: true
                Layout.fillHeight: true
                minimumContentWidth: head.totalWidth + 10
                TableHeader {
                    id: head
                    Layout.leftMargin: 10
                    columns: [{ t: "Beam", w: 60 }, { t: "Rule", w: 110 }, { t: "Component", w: 96 },
                              { t: "Status", w: 50 }, { t: "Margin", w: 80 }, { t: "Worst", w: 76 },
                              { t: "Required", w: 70 }]
                }
                ListView {
                    Layout.leftMargin: 10
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    clip: true
                    model: talisModel
                    boundsBehavior: Flickable.StopAtBounds
                    ScrollBar.vertical: ScrollBar {}
                    delegate: RowLayout {
                        spacing: 4
                        Tip { visible: rowHover.hovered && model.note !== ""; text: model.note }
                        HoverHandler { id: rowHover }
                        Cell { text: model.beam; Layout.preferredWidth: head.widths[0] }
                        Cell { text: model.rule; Layout.preferredWidth: head.widths[1] }
                        Cell { text: model.component; Layout.preferredWidth: head.widths[2] }
                        Cell { text: model.status; Layout.preferredWidth: head.widths[3]; color: statusColor(model.status) }
                        Cell { text: model.margin; Layout.preferredWidth: head.widths[4] }
                        Cell { text: model.worst; Layout.preferredWidth: head.widths[5]; opacity: 0.7 }
                        Cell { text: model.required; Layout.preferredWidth: head.widths[6]; opacity: 0.7 }
                    }
                }
            }
        }
    }
}
