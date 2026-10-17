// Clearance report: the closest beam per component and every beam–component row,
// colored by status, with CSV export.
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
        visible: !app.rvHasReport
        text: app.rvRunning ? "checked " + app.rvDone + " of " + app.rvTotal + " beams"
              : "Set up the clearance check in Properties (select it in the model tree) and run it."
        opacity: 1
    }
    RowLayout {
        visible: app.rvHasReport
        spacing: 8
        Label { text: app.rvBrief; font.weight: Font.DemiBold; Layout.fillWidth: true; elide: Text.ElideRight }
        Label { text: "beams changed since"; visible: app.rvStale; color: "#d08a00" }
        Button { text: "Export CSV…"; enabled: !app.rvRunning; onClicked: win.openClearanceExport() }
    }
    Repeater {
        model: worstModel
        delegate: Label {
            text: model.component + ": " + model.beam + " · " +
                  (model.worst === "–" ? "no hardware point before the landing"
                                       : "clearance " + model.worst + " mm (nominal " + model.nominal + " mm)") +
                  " · " + model.power + " TW · " + model.status
            color: statusColor(model.status)
            elide: Text.ElideRight
            Layout.fillWidth: true
        }
    }
    RowLayout {
        visible: app.rvHasReport
        spacing: 4
        Repeater {
            model: [{ t: "Beam", w: 64 }, { t: "Component", w: 110 }, { t: "Status", w: 56 },
                    { t: "Worst (mm)", w: 80 }, { t: "Nominal", w: 72 },
                    { t: "Power (TW)", w: 80 }, { t: "Peak (W/cm²)", w: 96 }]
            Label { text: modelData.t; Layout.preferredWidth: modelData.w; font.weight: Font.DemiBold }
        }
    }
    ListView {
        visible: app.rvHasReport
        Layout.fillWidth: true
        Layout.fillHeight: true
        Layout.minimumHeight: 60
        clip: true
        model: clearModel
        boundsBehavior: Flickable.StopAtBounds
        ScrollBar.vertical: ScrollBar {}
        delegate: RowLayout {
            spacing: 4
            Label { text: model.beam; Layout.preferredWidth: 64 }
            Label { text: model.component; Layout.preferredWidth: 110; elide: Text.ElideRight }
            Label { text: model.status; Layout.preferredWidth: 56; color: statusColor(model.status) }
            Label { text: model.worst; Layout.preferredWidth: 80 }
            Label { text: model.nominal; Layout.preferredWidth: 72; opacity: 0.7 }
            Label { text: model.power; Layout.preferredWidth: 80 }
            Label { text: model.intensity; Layout.preferredWidth: 96 }
        }
    }
}
