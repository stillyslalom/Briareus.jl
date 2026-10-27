// Drive results: capsule modes, foil spot, laser power accounting, VisRad benchmark kit.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ScrollPage {
    Note {
        visible: !app.driveHasResult
        text: app.driveBuilding ? "Computing the view factors… " + app.driveElapsed
              : "Compute the drive (Compute ▸ Drive, or color the target by a drive field) to see " +
                (app.hasFoil ? "the foil spot" : app.hasCapsule ? "the capsule modes" : "the drive") +
                " and the laser power accounting."
    }
    Label { text: app.driveSummary; font.weight: Font.DemiBold; visible: app.driveHasResult; wrapMode: Text.WordWrap; Layout.fillWidth: true }
    RowLayout {
        visible: app.driveHasResult
        spacing: 40
        Layout.fillWidth: true
        ColumnLayout {
            visible: app.driveMean !== ""
            Layout.alignment: Qt.AlignTop
            spacing: 4
            SectionLabel { text: "Capsule drive"; Layout.topMargin: 0 }
            Label { text: "Mean incident flux " + app.driveMean }
            GridLayout {
                columns: 3
                columnSpacing: 24
                rowSpacing: 2
                Label { text: "l"; font.weight: Font.DemiBold }
                Label { text: "ε_l"; font.weight: Font.DemiBold }
                Label { text: "P_l / P0"; font.weight: Font.DemiBold }
                Repeater {
                    model: modeModel
                    delegate: Label { Layout.row: model.l; Layout.column: 0; text: model.l }
                }
                Repeater {
                    model: modeModel
                    delegate: Label { Layout.row: model.l; Layout.column: 1; text: model.epsilon }
                }
                Repeater {
                    model: modeModel
                    delegate: Label { Layout.row: model.l; Layout.column: 2; text: model.legendre }
                }
            }
            Note { text: "ε_l: rotation-invariant amplitude over all m (= |P_l/P0| when axisymmetric)." }
        }
        ColumnLayout {
            visible: app.spotText !== ""
            Layout.alignment: Qt.AlignTop
            spacing: 4
            SectionLabel { text: "Foil spot"; Layout.topMargin: 0 }
            Label { text: app.spotText; wrapMode: Text.WordWrap; Layout.maximumWidth: 420 }
        }
        ColumnLayout {
            Layout.alignment: Qt.AlignTop
            spacing: 4
            SectionLabel { text: "Laser power"; Layout.topMargin: 0 }
            Label { text: app.delivered; opacity: 0.75 }
            GridLayout {
                columns: 3
                columnSpacing: 24
                rowSpacing: 2
                Repeater {
                    model: powerModel
                    delegate: Label { Layout.row: index; Layout.column: 0; text: model.name }
                }
                Repeater {
                    model: powerModel
                    delegate: Label { Layout.row: index; Layout.column: 1; text: model.power + " TW" }
                }
                Repeater {
                    model: powerModel
                    delegate: Label { Layout.row: index; Layout.column: 2; text: model.percent }
                }
            }
        }
        Item { Layout.fillWidth: true }
    }
    RowLayout {
        spacing: 8
        Layout.topMargin: 6
        Button {
            text: "Write VisRad benchmark kit…"
            enabled: app.kitReady && !app.kitRunning
            onClicked: win.openKitDialog()
            ToolTip.visible: hovered
            ToolTip.text: "Write the pointing (SST), the beam summary, the beam table, a README with " +
                          "the geometry and the VisRad steps, and these results to compare against"
        }
        BusyIndicator { running: app.kitRunning; visible: running; Layout.preferredHeight: 28; Layout.preferredWidth: 28 }
        Note { text: app.kitStatus }
    }
}
