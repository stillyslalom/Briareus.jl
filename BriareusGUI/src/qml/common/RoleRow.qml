// A component's role in a check (the clearance check's, the TaLIS checks'): its name,
// cut short with the whole name on hover, and a choice of `roles` ({ text, value, tip }
// each); `role` is the current value and `chosen(value)` fires on a change. The box's
// tooltip is the current role's `tip`.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import "."

RowLayout {
    id: row
    property string label: ""
    property string role: ""
    property var roles: []
    property int labelWidth: 120
    signal chosen(string value)
    spacing: 8
    Layout.fillWidth: true
    readonly property int current: Math.max(0, roles.map(r => r.value).indexOf(role))
    Cell { text: row.label; Layout.preferredWidth: row.labelWidth }
    ComboBox {
        Layout.fillWidth: true
        textRole: "text"; valueRole: "value"
        model: row.roles
        currentIndex: row.current
        onActivated: row.chosen(currentValue)
        tip: row.roles.length > 0 ? row.roles[row.current].tip : ""
    }
}
