// The Review step's TaLIS tab: capture and checked components, the rule tolerance and
// its interpretation, the check (on a worker), the per-rule summary, every check
// colored by status, and the report's notes (the ASSUMED tolerance interpretation
// always first and visible).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml

ColumnLayout {
    id: section
    spacing: 12

    function statusColor(s) {
        if (s === "fail") return "#c42b1c"
        if (s === "review") return "#d08a00"
        return "#2e9e44"
    }

    Label { text: "Components"; font.weight: Font.DemiBold }
    ColumnLayout {
        spacing: 2
        enabled: !app.tlRunning
        Repeater {
            model: talisCompModel
            delegate: RowLayout {
                id: compRow
                readonly property string compName: model.name
                readonly property string compRole: model.role
                spacing: 16
                Label { text: compRow.compName; Layout.preferredWidth: 120 }
                ComboBox {
                    Layout.preferredWidth: 260
                    textRole: "text"; valueRole: "value"
                    model: [{ text: "capture: beams land on it", value: "capture" },
                            { text: "check 3ω and 1ω clearance", value: "hardware" },
                            { text: "ignore", value: "none" }]
                    currentIndex: Math.max(0, ["capture", "hardware", "none"].indexOf(compRow.compRole))
                    onActivated: Julia.br_tl_role(compRow.compName, currentValue)
                }
            }
        }
    }

    Label { text: "Rule tolerance"; font.weight: Font.DemiBold; Layout.topMargin: 8 }
    Label {
        text: "The angular tolerance becomes a lateral allowance r·tan(tolerance) at each " +
              "hardware point; choose the point r is measured from (a point: e.g. a " +
              "component's mount)."
        opacity: 0.75
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }
    GridLayout {
        columns: 3
        columnSpacing: 12
        rowSpacing: 6
        enabled: !app.tlRunning
        Label { text: "Tolerance"; Layout.preferredWidth: 160 }
        TextField {
            Layout.preferredWidth: 90
            text: app.tlTolerance
            onEditingFinished: if (text !== app.tlTolerance) Julia.br_tl_edit("tolerance", text)
        }
        Label { text: "degrees (0: none)"; opacity: 0.75 }
        Label { text: "Distance r from" }
        ComboBox {
            id: pivotBox
            Layout.preferredWidth: 260
            Layout.columnSpan: 2
            textRole: "text"; valueRole: "value"
            model: [{ text: "TCC", value: "tcc" },
                    { text: "the beam's aim point", value: "aim" },
                    { text: "the farther of TCC and the aim", value: "max" },
                    { text: "a point (x, y, z mm)", value: "point" }]
            currentIndex: Math.max(0, ["tcc", "aim", "max", "point"].indexOf(app.tlPivot))
            onActivated: Julia.br_tl_pivot(currentValue)
        }
        Label { text: "Pivot point"; visible: app.tlPivot === "point" }
        TextField {
            visible: app.tlPivot === "point"
            Layout.preferredWidth: 160
            text: app.tlPivotPoint
            onEditingFinished: if (text !== app.tlPivotPoint) Julia.br_tl_edit("pivot_point", text)
        }
        Label { text: "x, y, z (mm)"; opacity: 0.75; visible: app.tlPivot === "point" }
        Label { text: "1ω rules" }
        Switch { checked: app.tlOmega1; onToggled: Julia.br_tl_omega1(checked) }
        Label { text: "NIF beams only"; opacity: 0.75 }
    }
    Label {
        text: "Error set (from the Clearance tab): pointing " + app.rv_pointing + " mm, positioning " +
              app.rv_target_offset + " mm, hardware " + app.rv_hardware_tolerance + " mm, focus ± " +
              app.rv_focus + " mm; " + app.rvQuality + " sampling."
        opacity: 0.75
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }
    Label {
        text: app.tlError
        visible: text !== ""
        color: "#c42b1c"
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }

    RowLayout {
        spacing: 8
        Layout.topMargin: 8
        Button {
            text: "Run TaLIS checks"
            highlighted: true
            enabled: !app.tlRunning
            onClicked: Julia.br_tl_start()
        }
        Button { text: "Cancel"; enabled: app.tlRunning; onClicked: Julia.br_tl_cancel() }
        Button {
            text: "Export CSV…"
            enabled: app.tlHasReport && !app.tlRunning
            onClicked: talisExport.open()
        }
        BusyIndicator {
            running: app.tlRunning
            visible: running
            Layout.preferredHeight: 32
            Layout.preferredWidth: 32
        }
    }
    Label {
        text: app.tlRunning ? "checking… " + app.tlElapsed :
              app.tlStale ? "The beams or the target changed since this check. " + app.tlBrief : app.tlStatus
        visible: text !== ""
        color: app.tlStale && !app.tlRunning ? "#d08a00" : palette.windowText
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }

    // notes: the assumed tolerance interpretation first, always shown
    Label { text: "Notes"; font.weight: Font.DemiBold; Layout.topMargin: 8 }
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

    Label { text: "Rules"; font.weight: Font.DemiBold; visible: app.tlHasReport; Layout.topMargin: 8 }
    RowLayout {
        visible: app.tlHasReport
        spacing: 4
        Repeater {
            model: [{ t: "Rule", w: 118 }, { t: "pass", w: 36 }, { t: "review", w: 46 },
                    { t: "fail", w: 32 }, { t: "Worst margin", w: 92 }, { t: "Beam / component", w: 130 }]
            Label { text: modelData.t; Layout.preferredWidth: modelData.w; font.weight: Font.DemiBold }
        }
    }
    Repeater {
        model: ruleModel
        delegate: RowLayout {
            spacing: 4
            Label {
                text: model.label
                Layout.preferredWidth: 118
                color: statusColor(model.status)
                ToolTip.visible: ruleHover.hovered
                ToolTip.text: model.requirement
                HoverHandler { id: ruleHover }
            }
            Label { text: model.pass; Layout.preferredWidth: 36 }
            Label { text: model.review; Layout.preferredWidth: 46; color: model.review > 0 ? "#d08a00" : palette.windowText }
            Label { text: model.fail; Layout.preferredWidth: 32; color: model.fail > 0 ? "#c42b1c" : palette.windowText }
            Label { text: model.worst; Layout.preferredWidth: 92 }
            Label { text: model.where_; Layout.preferredWidth: 130; elide: Text.ElideRight; opacity: 0.75 }
        }
    }

    Label { text: "Every check"; font.weight: Font.DemiBold; visible: app.tlHasReport; Layout.topMargin: 8 }
    Label {
        text: "Margin: worst case less the requirement (positive passes; n/a: nothing to measure, " +
              "e.g. a capture surface without a free edge or hardware downstream of the landing). " +
              "Worst, nominal and required in the rule's unit. Hover a row for its note."
        visible: app.tlHasReport
        opacity: 0.75
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }
    RowLayout {
        visible: app.tlHasReport
        spacing: 4
        Repeater {
            model: [{ t: "Beam", w: 56 }, { t: "Rule", w: 96 }, { t: "Component", w: 76 },
                    { t: "Status", w: 44 }, { t: "Margin", w: 72 }, { t: "Worst", w: 66 },
                    { t: "Req.", w: 60 }]
            Label { text: modelData.t; Layout.preferredWidth: modelData.w; font.weight: Font.DemiBold }
        }
    }
    ListView {
        visible: app.tlHasReport
        Layout.fillWidth: true
        Layout.preferredHeight: Math.min(contentHeight, 420)
        clip: true
        model: talisModel
        ScrollBar.vertical: ScrollBar {}
        delegate: RowLayout {
            spacing: 4
            ToolTip.visible: rowHover.hovered && model.note !== ""
            ToolTip.text: model.note
            HoverHandler { id: rowHover }
            Label { text: model.beam; Layout.preferredWidth: 56; elide: Text.ElideRight }
            Label { text: model.rule; Layout.preferredWidth: 96; elide: Text.ElideRight }
            Label { text: model.component; Layout.preferredWidth: 76; elide: Text.ElideRight }
            Label { text: model.status; Layout.preferredWidth: 44; color: statusColor(model.status) }
            Label { text: model.margin; Layout.preferredWidth: 72 }
            Label { text: model.worst; Layout.preferredWidth: 66; opacity: 0.75 }
            Label { text: model.required; Layout.preferredWidth: 60; opacity: 0.75 }
        }
    }

    FileDialog {
        id: talisExport
        title: "Export the TaLIS report"
        fileMode: FileDialog.SaveFile
        defaultSuffix: "csv"
        nameFilters: ["TaLIS report (*.csv)"]
        onAccepted: Julia.br_tl_export(selectedFile.toString())
    }
}
