// TaLIS report: notes (the ASSUMED tolerance interpretation first), per-rule summary,
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
    Note {
        visible: !app.tlHasReport
        text: app.tlRunning ? "checking… " + app.tlElapsed
              : "Set up the TaLIS checks in Properties (select them in the model tree) and run them."
        opacity: 1
    }
    RowLayout {
        visible: app.tlHasReport
        spacing: 8
        Label { text: app.tlBrief; font.weight: Font.DemiBold; Layout.fillWidth: true; elide: Text.ElideRight }
        Label { text: "beams changed since"; visible: app.tlStale; color: "#d08a00" }
        Button { text: "Export CSV…"; enabled: !app.tlRunning; onClicked: win.openTalisExport() }
    }
    // notes and the rule summary beside every check, split by a draggable handle; the
    // notes wrap at the width of their side
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
                Repeater {
                    model: noteModel
                    delegate: Label {
                        text: (model.kind === "assumed" ? "⚠ " : "· ") + model.text
                        color: model.kind === "assumed" ? "#d08a00" : palette.windowText
                        font.weight: model.kind === "assumed" ? Font.DemiBold : Font.Normal
                        opacity: model.kind === "scope" || model.kind === "source" ? 0.7 : 1
                        wrapMode: Text.WordWrap
                        Layout.fillWidth: true
                    }
                }
                SectionLabel { text: "Rules"; visible: app.tlHasReport }
                // one rule in two lines: its counts, then its worst check
                Repeater {
                    model: ruleModel
                    delegate: ColumnLayout {
                        spacing: 0
                        Layout.fillWidth: true
                        ToolTip.visible: ruleHover.hovered
                        ToolTip.delay: 400
                        ToolTip.text: model.requirement
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
            }
        }
        // every check
        ColumnLayout {
            visible: app.tlHasReport
            SplitView.fillWidth: app.tlHasReport
            SplitView.minimumWidth: 300
            spacing: 4
            RowLayout {
                spacing: 4
                Layout.leftMargin: 10
                Repeater {
                    model: [{ t: "Beam", w: 60 }, { t: "Rule", w: 110 }, { t: "Component", w: 96 },
                            { t: "Status", w: 50 }, { t: "Margin", w: 80 }, { t: "Worst", w: 76 },
                            { t: "Required", w: 70 }]
                    Label { text: modelData.t; Layout.preferredWidth: modelData.w; font.weight: Font.DemiBold }
                }
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
                    ToolTip.visible: rowHover.hovered && model.note !== ""
                    ToolTip.text: model.note
                    HoverHandler { id: rowHover }
                    Label { text: model.beam; Layout.preferredWidth: 60; elide: Text.ElideRight }
                    Label { text: model.rule; Layout.preferredWidth: 110; elide: Text.ElideRight }
                    Label { text: model.component; Layout.preferredWidth: 96; elide: Text.ElideRight }
                    Label { text: model.status; Layout.preferredWidth: 50; color: statusColor(model.status) }
                    Label { text: model.margin; Layout.preferredWidth: 80 }
                    Label { text: model.worst; Layout.preferredWidth: 76; opacity: 0.7 }
                    Label { text: model.required; Layout.preferredWidth: 70; opacity: 0.7 }
                }
            }
        }
    }
}
