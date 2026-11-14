// Drive results: capsule modes, foil spot, laser power accounting.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ScrollPage {
    RowLayout {
        visible: !app.driveHasResult
        spacing: 6
        Note {
            Layout.fillWidth: false
            text: app.driveBuilding ? "Computing the view factors… " + app.driveElapsed : "No drive yet: Compute ▸ Drive"
        }
        Info { topic: "drive.intro" }
    }
    RowLayout {
        visible: app.driveHasResult
        spacing: 40
        Layout.fillWidth: true
        // the observed surface's drive, and under it a hohlraum part's wall
        ColumnLayout {
            visible: app.driveMean !== "" || app.cavityNames.length > 0
            Layout.alignment: Qt.AlignTop
            spacing: 4
            ColumnLayout {
                visible: app.driveMean !== ""
                spacing: 4
                SectionLabel {
                    text: app.driveFace ? "Observed face" : "Capsule drive"
                    info: app.driveFace ? "" : "drive.modes"
                    Layout.topMargin: 0
                }
                Label { text: "Mean incident flux " + app.driveMean }
                GridLayout {
                    visible: !app.driveFace
                    columns: 3
                    columnSpacing: 24
                    rowSpacing: 2
                    Label { text: "ℓ"; font.weight: Font.DemiBold }
                    Label { text: "εₗ"; font.weight: Font.DemiBold }
                    Label { text: "Pₗ/P₀"; font.weight: Font.DemiBold }
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
            }
            ColumnLayout {
                visible: app.cavityNames.length > 0
                spacing: 4
                SectionLabel { text: "Hohlraum wall"; info: "drive.walls"; Layout.topMargin: app.driveMean !== "" ? 12 : 0 }
                ComboBox {
                    id: cavityBox
                    visible: app.cavityNames.length > 1
                    Layout.preferredWidth: Math.ceil(cavityBox.widestEntry) + 48
                    model: app.cavityNames
                    currentIndex: app.cavityIndex
                    onActivated: (index) => Julia.br_drive_cavity(index)
                    Tip {
                        visible: parent.hovered
                        text: "The hohlraum whose wall temperature is shown"
                    }
                }
                Label { text: "Mean incident flux " + app.cavityTr }
            }
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
}
