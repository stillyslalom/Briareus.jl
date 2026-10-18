// The bottom panel: the beam table and each analysis's results, one tab each.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "panels"

Pane {
    id: results
    padding: 0
    property alias currentIndex: tabs.currentIndex
    readonly property var keys: ["beams", "drive", "optimize", "clearance", "talis"]
    function show(key) { const i = keys.indexOf(key); if (i >= 0) tabs.currentIndex = i }

    ColumnLayout {
        anchors.fill: parent
        spacing: 0
        TabBar {
            id: tabs
            Layout.fillWidth: true
            TabButton { text: "Beam table"; width: implicitWidth }
            TabButton { text: "Drive"; width: implicitWidth }
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
