// The bottom panel: the beam table and each analysis's results, one tab each (Time:
// pulses, drive histories and the retune over time).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "panels"

Pane {
    id: results
    padding: 0
    property alias currentIndex: tabs.currentIndex
    readonly property var keys: ["beams", "drive", "time", "optimize", "clearance", "talis"]
    function show(key) { const i = keys.indexOf(key); if (i >= 0) tabs.currentIndex = i }
    function redrawTime() { timePanel.redraw() }
    // leaving NIF hides the TaLIS tab: show the drive rather than its stale page
    readonly property bool talisShown: app.facility === "nif"
    onTalisShownChanged: if (!talisShown && tabs.currentIndex === keys.indexOf("talis")) show("drive")

    ColumnLayout {
        anchors.fill: parent
        spacing: 0
        TabBar {
            id: tabs
            Layout.fillWidth: true
            TabButton { text: "Beam table"; width: implicitWidth }
            TabButton { text: "Drive"; width: implicitWidth }
            TabButton { text: "Time"; width: implicitWidth }
            TabButton { text: "Optimization"; width: implicitWidth }
            TabButton { text: "Clearance"; width: implicitWidth }
            // NIF's rules: the tab only for NIF beams
            TabButton { text: "TaLIS"; width: app.facility === "nif" ? implicitWidth : 0; visible: app.facility === "nif" }
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
            Item {
                TalisPanel { anchors.fill: parent; anchors.margins: 10 }
            }
        }
    }
}
