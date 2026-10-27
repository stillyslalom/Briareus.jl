// Optimization results: objective, modes before and after, clipping and capsule
// clearance, the parameter table, and applying the result.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ScrollPage {
    Note {
        visible: !app.optHasResult
        text: app.optRunning ? app.optProgress : "Set up the optimization in Properties (select it in the model tree) and run it."
        opacity: 1
    }
    Label { text: app.optValues; visible: app.optHasResult; wrapMode: Text.WordWrap; Layout.fillWidth: true; font.weight: Font.DemiBold }
    Note { text: app.optClipped; visible: app.optHasResult }
    Note { text: app.optSpot }
    RowLayout {
        visible: app.optHasResult
        spacing: 40
        Layout.fillWidth: true
        GridLayout {
            Layout.alignment: Qt.AlignTop
            columns: 3
            columnSpacing: 24
            rowSpacing: 2
            visible: optModeRepeater.count > 0
            Label { text: "l"; font.weight: Font.DemiBold }
            Label { text: "εₗ before"; font.weight: Font.DemiBold }
            Label { text: "εₗ after"; font.weight: Font.DemiBold }
            Repeater {
                id: optModeRepeater
                model: optModeModel
                delegate: Label { Layout.row: model.l; Layout.column: 0; text: model.l }
            }
            Repeater {
                model: optModeModel
                delegate: Label { Layout.row: model.l; Layout.column: 1; text: model.before }
            }
            Repeater {
                model: optModeModel
                delegate: Label { Layout.row: model.l; Layout.column: 2; text: model.after }
            }
        }
        GridLayout {
            Layout.alignment: Qt.AlignTop
            columns: 5
            columnSpacing: 16
            rowSpacing: 2
            Label { text: "Parameter"; font.weight: Font.DemiBold }
            Label { text: "Start"; font.weight: Font.DemiBold }
            Label { text: "Optimized"; font.weight: Font.DemiBold }
            Label { text: "Lower"; font.weight: Font.DemiBold }
            Label { text: "Upper"; font.weight: Font.DemiBold }
            Repeater {
                model: paramModel
                delegate: Label { Layout.row: index + 1; Layout.column: 0; text: model.name }
            }
            Repeater {
                model: paramModel
                delegate: Label { Layout.row: index + 1; Layout.column: 1; text: model.start }
            }
            Repeater {
                model: paramModel
                delegate: Label {
                    Layout.row: index + 1; Layout.column: 2
                    text: model.optimized + (model.atBound ? "  (at bound)" : "")
                    color: model.atBound ? "#d08a00" : palette.windowText
                }
            }
            Repeater {
                model: paramModel
                delegate: Label { Layout.row: index + 1; Layout.column: 3; text: model.lower; opacity: 0.7 }
            }
            Repeater {
                model: paramModel
                delegate: Label { Layout.row: index + 1; Layout.column: 4; text: model.upper; opacity: 0.7 }
            }
        }
        Item { Layout.fillWidth: true }
    }
    Button {
        visible: app.optHasResult
        text: app.optApplied ? "Applied to the beam table" : "Apply to the beam table"
        enabled: !app.optApplied && !app.optRunning
        onClicked: Julia.br_opt_accept()
    }
}
