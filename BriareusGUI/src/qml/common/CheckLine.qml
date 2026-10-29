// A check box whose label wraps to the panel width (a CheckBox keeps its text on one
// line). A click on the box or on the label toggles it; `toggled()` then fires, with
// `checked` already changed, as for a CheckBox. `tip` shows on hover.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import "."

RowLayout {
    id: line
    property alias text: label.text
    property alias checked: box.checked
    property string tip: ""
    readonly property bool hovered: lineHover.hovered
    signal toggled()
    spacing: 0
    Layout.fillWidth: true
    HoverHandler { id: lineHover }
    ToolTip.visible: tip !== "" && hovered
    ToolTip.delay: 400
    ToolTip.text: tip
    CheckBox {
        id: box
        Layout.alignment: Qt.AlignTop
        rightPadding: 0
        spacing: 0
        // as wide as its box: without a text, a style may center the box in a wider control
        Layout.preferredWidth: leftPadding + (indicator ? indicator.width : 20) + 8
        onToggled: line.toggled()
    }
    Label {
        id: label
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
        Layout.alignment: Qt.AlignVCenter
        TapHandler { onTapped: { box.toggle(); line.toggled() } }
    }
}
