// The bottom panel shows the beam table and the results of each analysis, one tab each.
// The Time tab shows pulses, drive histories and the retune over time, the Diagnostics tab
// the synthetic diagnostics a history records, and the Sweep tab the steps of a placement
// sweep.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "panels"

Pane {
    id: results
    padding: 0
    property alias currentIndex: tabs.currentIndex
    readonly property var keys: ["beams", "drive", "time", "diagnostics", "optimize", "clearance", "sweep", "talis"]
    function show(key) { const i = keys.indexOf(key); if (i >= 0) tabs.currentIndex = i }
    function redrawTime() { timePanel.redraw() }
    function redrawDiagnostics() { diagnosticsPanel.redraw() }
    function openDiagnosticPicker() { return diagnosticsPanel.openPicker() }
    function redrawSweep() { sweepPanel.redraw() }
    // The height the current tab needs, which the window grows the panel to. The Time and
    // Diagnostics tabs need their rows and the chart at its minimum height; the others need
    // 0, since they fit any height.
    readonly property real fitHeight: tabs.currentIndex === keys.indexOf("time") ? tabs.height + timePanel.fitHeight
                                    : tabs.currentIndex === keys.indexOf("diagnostics") ? tabs.height + diagnosticsPanel.fitHeight : 0
    // Leaving NIF hides the TaLIS tab, so the panel shows the drive instead of the stale
    // TaLIS page.
    readonly property bool talisShown: app.facility === "nif"
    onTalisShownChanged: if (!talisShown && tabs.currentIndex === keys.indexOf("talis")) show("drive")

    ColumnLayout {
        anchors.fill: parent
        spacing: 0
        TabBar {
            id: tabs
            Layout.fillWidth: true
            // The side padding of the tabs, less when the panel is narrow, so that every tab
            // shows its whole name instead of being scrolled out of sight. It is computed
            // from the widths of the names, which do not depend on the padding.
            readonly property var names: ["Beam table", "Drive", "Time", "Diagnostics", "Optimization", "Clearance", "Sweep", "TaLIS"]
            FontMetrics { id: tabFont; font: tabs.font }
            readonly property real pad: {
                const shown = app.facility === "nif" ? names : names.slice(0, -1)   // TaLIS needs NIF beams
                const w = shown.reduce((s, t) => s + tabFont.advanceWidth(t), 0)
                // less the spacing between tabs and a margin for the bar padding
                const room = results.width - w - spacing * (shown.length - 1) - leftPadding - rightPadding - 24
                return Math.max(4, Math.min(18, Math.floor(room / (2 * shown.length))))
            }
            topPadding: 0
            bottomPadding: 0
            implicitHeight: 32
            component Tab: TabButton {
                leftPadding: tabs.pad
                rightPadding: tabs.pad
                topPadding: 4
                bottomPadding: 6
                implicitHeight: 30
                width: Math.ceil(tabFont.advanceWidth(text)) + 2 * tabs.pad + 2
            }
            Tab { text: "Beam table" }
            Tab { text: "Drive" }
            Tab { text: "Time" }
            Tab { text: "Diagnostics" }
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
                BeamTablePanel { anchors.fill: parent; anchors.margins: 6 }
            }
            DrivePanel {}
            TimePanel { id: timePanel }
            DiagnosticsPanel { id: diagnosticsPanel }
            OptimizePanel {}
            Item {
                ClearancePanel { anchors.fill: parent; anchors.margins: 6 }
            }
            SweepPanel { id: sweepPanel }
            Item {
                TalisPanel { anchors.fill: parent; anchors.margins: 6 }
            }
        }
    }
}
