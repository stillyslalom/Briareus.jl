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
        Cell { text: app.rvBrief; font.weight: Font.DemiBold; Layout.fillWidth: true }
        Label { text: "inputs changed since"; visible: app.rvStale; color: "#d08a00" }
        Button { text: "Export CSV…"; enabled: !app.rvRunning; onClicked: win.openClearanceExport() }
    }
    Repeater {
        model: worstModel
        delegate: Cell {
            text: model.component + ": " + model.beam + " · " +
                  (model.worst === "–" ? "no hardware point before the landing"
                                       : "clearance " + model.worst + " mm (nominal " + model.nominal + " mm)") +
                  " · " + model.power + " TW · " + model.status
            color: statusColor(model.status)
            Layout.fillWidth: true
        }
    }
    HScroll {
        visible: app.rvHasReport
        Layout.fillWidth: true
        Layout.fillHeight: true
        Layout.minimumHeight: 80
        minimumContentWidth: head.totalWidth
        TableHeader {
            id: head
            columns: [{ t: "Beam", w: 64 }, { t: "Component", w: 110 }, { t: "Status", w: 56 },
                      { t: "Worst (mm)", w: 80 }, { t: "Nominal", w: 72 },
                      { t: "Power (TW)", w: 80 }, { t: "Peak (W/cm²)", w: 96 }]
        }
        ListView {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Layout.minimumHeight: 60
            clip: true
            model: clearModel
            boundsBehavior: Flickable.StopAtBounds
            ScrollBar.vertical: ScrollBar {}
            delegate: RowLayout {
                spacing: 4
                Cell { text: model.beam; Layout.preferredWidth: head.widths[0] }
                Cell { text: model.component; Layout.preferredWidth: head.widths[1] }
                Cell { text: model.status; Layout.preferredWidth: head.widths[2]; color: statusColor(model.status) }
                Cell { text: model.worst; Layout.preferredWidth: head.widths[3] }
                Cell { text: model.nominal; Layout.preferredWidth: head.widths[4]; opacity: 0.7 }
                Cell { text: model.power; Layout.preferredWidth: head.widths[5] }
                Cell { text: model.intensity; Layout.preferredWidth: head.widths[6] }
            }
        }
    }
}
