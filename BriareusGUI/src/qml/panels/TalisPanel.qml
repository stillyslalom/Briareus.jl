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
    RowLayout {
        Layout.fillWidth: true
        Layout.fillHeight: true
        spacing: 24
        // notes and the rule summary: narrower than 420 when the panel is, so that the
        // checks' table (and the Export button above) stay inside it
        ScrollView {
            Layout.fillWidth: true
            Layout.minimumWidth: 260
            Layout.maximumWidth: 420
            Layout.preferredWidth: 420
            Layout.fillHeight: true
            contentWidth: availableWidth
            clip: true
            ColumnLayout {
                width: parent.width
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
                Repeater {
                    model: ruleModel
                    delegate: RowLayout {
                        spacing: 6
                        Label {
                            text: model.label
                            Layout.preferredWidth: 130
                            color: statusColor(model.status)
                            ToolTip.visible: ruleHover.hovered
                            ToolTip.text: model.requirement
                            HoverHandler { id: ruleHover }
                        }
                        Label { text: model.pass + " pass"; Layout.preferredWidth: 56 }
                        Label { text: model.review + " review"; Layout.preferredWidth: 64; color: model.review > 0 ? "#d08a00" : palette.windowText }
                        Label { text: model.fail + " fail"; Layout.preferredWidth: 46; color: model.fail > 0 ? "#c42b1c" : palette.windowText }
                        Label { text: model.worst; Layout.preferredWidth: 90 }
                        Label { text: model.where_; elide: Text.ElideRight; opacity: 0.7; Layout.fillWidth: true }
                    }
                }
            }
        }
        // every check
        ColumnLayout {
            visible: app.tlHasReport
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 4
            RowLayout {
                spacing: 4
                Repeater {
                    model: [{ t: "Beam", w: 60 }, { t: "Rule", w: 110 }, { t: "Component", w: 90 },
                            { t: "Status", w: 50 }, { t: "Margin", w: 80 }, { t: "Worst", w: 76 },
                            { t: "Required", w: 70 }]
                    Label { text: modelData.t; Layout.preferredWidth: modelData.w; font.weight: Font.DemiBold }
                }
            }
            ListView {
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
                    Label { text: model.component; Layout.preferredWidth: 90; elide: Text.ElideRight }
                    Label { text: model.status; Layout.preferredWidth: 50; color: statusColor(model.status) }
                    Label { text: model.margin; Layout.preferredWidth: 80 }
                    Label { text: model.worst; Layout.preferredWidth: 76; opacity: 0.7 }
                    Label { text: model.required; Layout.preferredWidth: 70; opacity: 0.7 }
                }
            }
        }
    }
}
