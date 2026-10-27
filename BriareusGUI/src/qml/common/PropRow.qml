// One editable property: label, text field (committed on Enter or when focus leaves),
// unit. `value` is the current text from Julia; `edited(text)` fires on a change.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

RowLayout {
    id: row
    property string label: ""
    property string value: ""
    property string unit: ""
    property string tip: ""
    property int labelWidth: 140
    property int fieldWidth: 90
    signal edited(string text)
    spacing: 8
    Label {
        text: row.label
        Layout.preferredWidth: row.labelWidth
        elide: Text.ElideRight
        // the tip, or the whole label when it is cut short
        ToolTip.visible: (row.tip !== "" || truncated) && labelHover.hovered
        ToolTip.delay: 400
        ToolTip.text: row.tip !== "" ? (truncated ? row.label + ": " : "") + row.tip : row.label
        HoverHandler { id: labelHover }
    }
    TextField {
        Layout.preferredWidth: row.fieldWidth
        Layout.minimumWidth: 60
        text: row.value
        selectByMouse: true
        onEditingFinished: if (text !== row.value) row.edited(text)
    }
    Label { text: row.unit; opacity: 0.7; Layout.fillWidth: true; elide: Text.ElideRight }
}
