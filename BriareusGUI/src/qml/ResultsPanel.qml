// The bottom panel: the beam table and each analysis's results, one tab each (Time:
// pulses, drive histories and the retune over time; Sweep: a placement sweep's steps).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "panels"

Pane {
    id: results
    padding: 0
    property alias currentIndex: tabs.currentIndex
    readonly property var keys: ["beams", "drive", "time", "optimize", "clearance", "sweep", "talis"]
    function show(key) { const i = keys.indexOf(key); if (i >= 0) tabs.currentIndex = i }
    function redrawTime() { timePanel.redraw() }
    function redrawSweep() { sweepPanel.redraw() }
    // leaving NIF hides the TaLIS tab: show the drive rather than its stale page
    readonly property bool talisShown: app.facility === "nif"
    onTalisShownChanged: if (!talisShown && tabs.currentIndex === keys.indexOf("talis")) show("drive")

    ColumnLayout {
        anchors.fill: parent
        spacing: 0
        TabBar {
            id: tabs
            Layout.fillWidth: true
            // the tabs' side padding, less when the panel is narrow, so that every tab shows
            // its whole name rather than being scrolled out of sight (from the names' widths,
            // which do not depend on the padding)
            readonly property var names: ["Beam table", "Drive", "Time", "Optimization", "Clearance", "Sweep", "TaLIS"]
            FontMetrics { id: tabFont; font: tabs.font }
            readonly property real pad: {
                const shown = app.facility === "nif" ? names : names.slice(0, -1)   // (TaLIS: NIF beams)
                const w = shown.reduce((s, t) => s + tabFont.advanceWidth(t), 0)
                // less the spacing between tabs and a margin for the bar's own padding
                const room = results.width - w - spacing * (shown.length - 1) - leftPadding - rightPadding - 24
                return Math.max(4, Math.min(18, Math.floor(room / (2 * shown.length))))
            }
            component Tab: TabButton {
                leftPadding: tabs.pad
                rightPadding: tabs.pad
                width: Math.ceil(tabFont.advanceWidth(text)) + 2 * tabs.pad + 2
            }
            Tab { text: "Beam table" }
            Tab { text: "Drive" }
            Tab { text: "Time" }
            Tab { text: "Optimization" }
            Tab { text: "Clearance" }
            Tab { text: "Sweep" }
            // NIF's rules: the tab only for NIF beams
            Tab { text: "TaLIS"; visible: app.facility === "nif"; width: visible ? Math.ceil(tabFont.advanceWidth(text)) + 2 * tabs.pad + 2 : 0 }
        }
        StackLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            currentIndex: tabs.currentIndex
            Item {
                BeamTablePanel { anchors.fill: parent; anchors.margins: 10 }
            }
            DrivePanel {}
            TimePanel { id: timePanel }
            OptimizePanel {}
            Item {
                ClearancePanel { anchors.fill: parent; anchors.margins: 10 }
            }
            SweepPanel { id: sweepPanel }
            Item {
                TalisPanel { anchors.fill: parent; anchors.margins: 10 }
            }
        }
    }
}
