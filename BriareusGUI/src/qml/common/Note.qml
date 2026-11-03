// A short status line (a summary, a state, an outcome) that wraps to the panel width;
// explanations go in help.js, shown by an ⓘ (Info.qml), or in `tip`, shown on hover.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Label {
    property string tip: ""
    opacity: 0.72
    wrapMode: Text.WordWrap
    Layout.fillWidth: true
    visible: text !== ""
    ToolTip.visible: tip !== "" && noteHover.hovered
    ToolTip.delay: 400
    ToolTip.text: tip
    HoverHandler { id: noteHover }
}
