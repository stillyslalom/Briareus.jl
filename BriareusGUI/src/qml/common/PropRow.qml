// One editable property: label, text field (committed on Enter or when focus leaves),
// unit. `value` is the current text from Julia; `edited(text)` fires on a change. The
// unit stands beside the field when it fits there, otherwise on a second, smaller line
// under the field, wrapped to the panel. With `check` set the label is a check box
// (`checked`, `toggled(on)`), and `fieldEnabled` can disable the field alone.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import "."

ColumnLayout {
    id: row
    property string label: ""
    property string value: ""
    property string unit: ""
    property string tip: ""
    property int labelWidth: 140
    property int fieldWidth: 90
    property bool check: false
    property bool checked: false
    property bool fieldEnabled: true
    signal edited(string text)
    signal toggled(bool on)
    spacing: 0
    Layout.fillWidth: true

    // the room beside the field, and whether the unit fits in it
    readonly property real unitRoom: width - labelWidth - fieldWidth - 2 * line.spacing
    readonly property bool unitBeside: unitMetrics.advanceWidth <= unitRoom
    TextMetrics { id: unitMetrics; font: unitInline.font; text: row.unit }

    RowLayout {
        id: line
        spacing: 8
        Layout.fillWidth: true
        Label {
            visible: !row.check
            text: row.label
            Layout.preferredWidth: row.labelWidth
            elide: Text.ElideRight
            // the tip, or the whole label when it is cut short
            ToolTip.visible: (row.tip !== "" || truncated) && labelHover.hovered
            ToolTip.delay: 400
            ToolTip.text: row.tip !== "" ? (truncated ? row.label + ": " : "") + row.tip : row.label
            HoverHandler { id: labelHover }
        }
        CheckBox {
            visible: row.check
            text: row.label
            checked: row.checked
            Layout.preferredWidth: row.labelWidth
            onToggled: row.toggled(checked)
            ToolTip.visible: row.tip !== "" && hovered
            ToolTip.delay: 400
            ToolTip.text: row.tip
        }
        TextField {
            id: field
            Layout.preferredWidth: row.fieldWidth
            Layout.minimumWidth: 60
            enabled: row.fieldEnabled
            text: row.value
            selectByMouse: true
            onEditingFinished: if (text !== row.value) row.edited(text)
        }
        Label {
            id: unitInline
            text: row.unitBeside ? row.unit : ""
            opacity: 0.7
            Layout.fillWidth: true
            elide: Text.ElideRight
        }
    }
    Label {
        visible: !row.unitBeside && row.unit !== ""
        text: row.unit
        opacity: 0.7
        font.pixelSize: 11
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
        Layout.leftMargin: row.labelWidth + line.spacing
    }
}
