// A part's placement steps, applied in order after its placement fields: per step a
// heading (number, kind, what it fixes; up, down, remove) and its fields (the frame it
// refers to, numbers, an axis or side, a flip), then a menu that adds a step and what
// stays free. One `stepModel` row per heading or field (`kind`: head, number, choice,
// check; `choices` joined by "|", `pick` the chosen index or the check state).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 4
    Layout.fillWidth: true
    RowLayout {
        Layout.fillWidth: true
        SectionLabel { text: "Placement steps"; Layout.fillWidth: true }
        Info { topic: "part.steps" }
    }
    Repeater {
        model: stepModel
        delegate: ColumnLayout {
            id: stepRow
            property string choicesText: model.choices
            property string rowKind: model.kind
            property int stepNumber: model.step
            property string rowKey: model.key
            property string rowValue: model.value
            property int rowPick: model.pick
            spacing: 0
            Layout.fillWidth: true
            Layout.topMargin: rowKind === "head" ? 6 : 0
            RowLayout {
                visible: stepRow.rowKind === "head"
                spacing: 4
                Layout.fillWidth: true
                Label { text: model.label; font.bold: true }
                Label {
                    text: stepRow.rowValue
                    opacity: 0.72
                    wrapMode: Text.WordWrap
                    Layout.fillWidth: true
                }
                Button {
                    text: "Sweep…"
                    onClicked: Julia.br_step_sweep(stepRow.stepNumber)
                    ToolTip.visible: hovered
                    ToolTip.text: "Step a number of this placement step through a range and check each value"
                }
                Button {
                    text: "↑"; Layout.preferredWidth: 30; enabled: stepRow.stepNumber > 1
                    onClicked: Julia.br_step_move(stepRow.stepNumber, -1)
                    ToolTip.visible: hovered; ToolTip.text: "Apply this step earlier"
                }
                Button {
                    text: "↓"; Layout.preferredWidth: 30
                    onClicked: Julia.br_step_move(stepRow.stepNumber, 1)
                    ToolTip.visible: hovered; ToolTip.text: "Apply this step later"
                }
                Button {
                    text: "✕"; Layout.preferredWidth: 30
                    onClicked: Julia.br_step_remove(stepRow.stepNumber)
                    ToolTip.visible: hovered; ToolTip.text: "Remove this step"
                }
            }
            PropRow {
                visible: stepRow.rowKind === "number"
                label: model.label
                value: stepRow.rowValue
                unit: model.unit
                labelWidth: 128
                Layout.leftMargin: 12
                onEdited: (text) => Julia.br_step_edit(stepRow.stepNumber, stepRow.rowKey, text)
            }
            RowLayout {
                visible: stepRow.rowKind === "choice"
                Layout.fillWidth: true
                Layout.leftMargin: 12
                Label { text: model.label; Layout.preferredWidth: 128 }
                ComboBox {
                    Layout.fillWidth: true
                    model: stepRow.choicesText.split("|")
                    currentIndex: stepRow.rowPick
                    onActivated: (index) => Julia.br_step_choice(stepRow.stepNumber, stepRow.rowKey, index)
                }
            }
            CheckLine {
                visible: stepRow.rowKind === "check"
                Layout.leftMargin: 12
                text: model.label
                checked: stepRow.rowPick === 1
                onToggled: Julia.br_step_choice(stepRow.stepNumber, stepRow.rowKey, checked ? 1 : 0)
            }
        }
    }
    ComboBox {
        Layout.fillWidth: true
        displayText: "Add a step…"
        model: app.stepKindLabels
        onActivated: (index) => Julia.br_step_add(index)
    }
    Note { text: app.placementNote; visible: app.placementNote !== "" }
}
